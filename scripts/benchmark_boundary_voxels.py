#!/usr/bin/env python3
"""Paired Octree vs sorted-voxel flat benchmark; identical trees and query sets."""
import argparse
import hashlib
import json
from pathlib import Path
import statistics

from benchmark_boundary_tree import REPO, execute, fixtures, save
from benchmark_octree_phase4 import write_xyz

MODES=("octree","voxels")

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--build-dir",type=Path,default=REPO/"build")
    parser.add_argument("--axis",choices=("x","y","z","best-axis"),default="y")
    parser.add_argument("--repeats",type=int,default=7)
    parser.add_argument("--warmups",type=int,default=1)
    parser.add_argument("--query-repeats",type=int,default=3)
    parser.add_argument("--cases",nargs="+")
    args=parser.parse_args()
    if args.repeats<1 or args.warmups<0 or args.query_repeats<1: parser.error("Invalid trial count")
    cases=fixtures()+[("aligned_gap",[(.02,.02,.02),(.10,.02,.02)])]
    if args.cases:
        if set(args.cases)-{n for n,_ in cases}: parser.error("Unknown fixture")
        cases=[(n,p) for n,p in cases if n in args.cases]
    output=args.output.resolve(); output.mkdir(parents=True,exist_ok=False)
    binary=args.build_dir.resolve()/"build_halfspace_tree"
    manifest={"binary_sha256":hashlib.sha256(binary.read_bytes()).hexdigest(),"axis":args.axis,"voxel_size":.04,
              "repeats":args.repeats,"warmups":args.warmups,"query_repeats":args.query_repeats,
              "order":"serial alternating pairs","modes":MODES,"inputs":[],
              "validation":"same fine-voxel reference tree, 1000 seeded random probes and complete sign diagnostics",
              "memory":"owned occupancy vector capacities + objects, excluding allocator overhead, partition, temporaries and RSS",
              "core_excludes":"certificates, probes, validation reference tree, input/output I/O and query benchmark"}
    rows=[]
    for name,points in cases:
        print(name,flush=True)
        folder=output/name; folder.mkdir(); cloud=folder/"input.xyz"; write_xyz(cloud,points)
        manifest["inputs"].append({"case":name,"raw_points":len(points),"sha256":hashlib.sha256(cloud.read_bytes()).hexdigest()})
        save(output/"manifest.json",manifest)
        def command(mode,out,query=False):
            cmd=[str(binary),str(cloud),str(out),"--pipeline","octree","--tree-source","boundary",
                 "--voxel-size",".04","--boundary-sweep",args.axis,"--boundary-expression","flat",
                 "--boundary-occupancy",mode,"--boundary-validation-reference","voxels",
                 "--boundary-diagnostics","true","--sign-propagation","scalar",
                 "--validation-samples","0" if query else "1000"]
            if query: cmd += ["--query-benchmark","true"]
            return cmd
        def read(out):
            m=json.loads((out/"benchmark.json").read_text())[0]
            v=json.loads((out/"validation_report.json").read_text())
            assert v["passed"] and m["validation_reference"]=="fine_voxels"
            assert m["false_free"]==m["missed_free"]==m["scalar_mismatches"]==m["raw_points_classified_free"]==0
            if m["occupancy_backend"]=="voxels":
                assert bool(m["octree_built"])==bool(m["fallback"]),"Unexpected Octree in normal voxel path"
                if not m["fallback"]: assert m["octree_nodes"]==m["baseline_partition_ms"]==m["lazy_fallback_ms"]==0
            return m
        def check_pair(a,b):
            ma,mb=read(a),read(b)
            keys=("fallback","selected_axis","tree_nodes","raw_tree_nodes","plane_references","occupied_voxels",
                  "free_prisms_after_merge","baseline_tree_nodes","tested_points","candidate_point_checks",
                  "sign_certificate_applied","strict_sign_certified","artificial_zero_faces","artificial_zero_edges",
                  "artificial_zero_vertices","raw_points_classified_free")
            for k in keys: assert ma[k]==mb[k],(name,k,ma[k],mb[k])
            entry=json.loads((a/"strict_free_space.json").read_text())
            for f in (entry["tree"],entry["raw_tree"],"boundary.off","boundary_supports.json","free_prisms.json"):
                assert (a/f).read_bytes()==(b/f).read_bytes(),(name,f)
            if (a/"query_samples.json").exists():
                assert (a/"query_samples.json").read_bytes()==(b/"query_samples.json").read_bytes()
                assert ma["query_zero_guard_fraction"]==mb["query_zero_guard_fraction"]
            return hashlib.sha256((a/entry["tree"]).read_bytes()).hexdigest()
        records={m:[] for m in MODES}; queries={m:[] for m in MODES}
        for mode in MODES: (folder/mode).mkdir()
        for trial in range(-args.warmups,args.repeats):
            order=MODES if trial%2==0 else MODES[::-1]
            for mode in order:
                out=folder/mode; wall=execute(command(mode,out),out/f"trial_{trial}.log"); m=read(out)
                if trial>=0: records[mode].append({"metrics":m,"process_wall_ms":wall})
            digest=check_pair(folder/"octree",folder/"voxels")
        for trial in range(args.query_repeats):
            for mode in (MODES if trial%2==0 else MODES[::-1]):
                out=folder/"query"/mode; out.mkdir(parents=True,exist_ok=True)
                execute(command(mode,out,True),out/f"query_{trial}.log"); queries[mode].append(read(out))
            check_pair(folder/"query"/"octree",folder/"query"/"voxels")
        for mode in MODES:
            save(folder/mode/"trials.json",records[mode]); save(folder/mode/"query_trials.json",queries[mode])
            m=records[mode][-1]["metrics"]
            med=lambda key: statistics.median(r["metrics"][key] for r in records[mode])
            qmed=lambda key: statistics.median(q[key] for q in queries[mode])
            row={"case":name,"occupancy":mode,"raw_points":len(points),"occupied_voxels":m["occupied_voxels"],
                 "tree_nodes":m["tree_nodes"],"plane_references":m["plane_references"],"tree_sha256":digest,
                 "core_ms":med("core_runtime_ms"),"core_min_ms":min(r["metrics"]["core_runtime_ms"] for r in records[mode]),
                 "core_max_ms":max(r["metrics"]["core_runtime_ms"] for r in records[mode]),
                 "occupancy_ms":med("occupancy_ms"),"rect_partition_ms":med("baseline_partition_ms"),
                 "lazy_fallback_ms":med("lazy_fallback_ms"),"build_validation_ms":med("total_ms"),
                 "validation_ms":med("validation_ms"),"occupancy_storage_bytes":m["occupancy_storage_bytes"],
                 "octree_nodes":m["octree_nodes"],"fallback":bool(m["fallback"]),
                 "strict_query_us":qmed("strict_free_query_us"),"occupancy_query_us":qmed("occupancy_query_us"),
                 "boundary_query_us":qmed("boundary_strict_free_query_us"),
                 "boundary_occupancy_query_us":qmed("boundary_occupancy_query_us"),
                 "seam_query_us":qmed("artificial_seam_strict_query_us"),
                 "seam_query_count":queries[mode][-1]["artificial_seam_query_count"],
                 "zero_guard_fraction":qmed("query_zero_guard_fraction"),
                 "sign_certificate_applied":bool(m["sign_certificate_applied"]),
                 "artificial_zeros_fev":[m["artificial_zero_faces"],m["artificial_zero_edges"],m["artificial_zero_vertices"]],
                 "byte_identical_geometry":True,"byte_identical_query_sets":True,"passed":True}
            rows.append(row)
        save(output/"summary.json",rows)
    report=["# Flat: Octree removal comparison","",
            f"axis={args.axis}; h=0.04 m; {args.warmups} warmup pairs, {args.repeats} measured pairs; serial alternating order.",
            f"Query: {args.query_repeats} paired processes, each median of 7 x 5 rounds; exported query corpora compared byte-for-byte.","",
            "Both modes use the same fine-voxel validation tree, seeded probes and sign diagnostics. Core excludes all validation and I/O.",
            "Raw/simplified trees, support registry, free prisms and boundary mesh must be byte-identical. Artificial zeros are unchanged.","",
            "| Case | Occupancy | Octree nodes | Tree nodes | Core ms | Occupancy ms | Rect ms | Build+validation ms | Storage KiB | Query us | Seam us | Fallback |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |"]
    for r in rows:
        report.append(f'| {r["case"]} | {r["occupancy"]} | {r["octree_nodes"]} | {r["tree_nodes"]} | '
                      f'{r["core_ms"]:.3f} | {r["occupancy_ms"]:.3f} | {r["rect_partition_ms"]:.3f} | '
                      f'{r["build_validation_ms"]:.3f} | {r["occupancy_storage_bytes"]/1024:.2f} | '
                      f'{r["strict_query_us"]:.3f} | {r["seam_query_us"]:.3f} | {r["fallback"]} |')
    report += ["","Storage is estimated owned occupancy payload (vector capacities + objects), NOT peak memory or process RSS.",
               "Voxel fallback includes BOTH sparse occupancy and lazily built legacy Octree/rect; lazy cost is included in core and retained in summary.",
               "The --pipeline octree CLI routing name and tree JSON pipeline tag are kept for compatibility; occupancy_backend identifies the actual implementation.",
               "Default remains octree. Only opt-in flat voxels mode removes hierarchy construction; no CBF/strict-zero-level-set claim."]
    (output/"report.md").write_text("\n".join(report)+"\n")
    print(output/"report.md",flush=True)

if __name__=="__main__": main()
