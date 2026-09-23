#!/usr/bin/env python3
"""Flat voxel sweep vs boundary-direct; common geometry and timed query corpora."""
import argparse
import hashlib
import json
from pathlib import Path
import statistics
from benchmark_boundary_tree import REPO, execute, fixtures, save
from benchmark_octree_phase4 import write_xyz

MODES=("flat","direct")
DEFAULT_CASES=("solid_cuboid","l_voxels","u_voxels","staircase","closed_cavity","two_disconnected",
               "aligned_gap","edge_touch","vertex_touch","camera")

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("output",type=Path); p.add_argument("--build-dir",type=Path,default=REPO/"build")
    p.add_argument("--repeats",type=int,default=5); p.add_argument("--warmups",type=int,default=1)
    p.add_argument("--query-repeats",type=int,default=1); p.add_argument("--cases",nargs="+",default=DEFAULT_CASES)
    args=p.parse_args()
    if args.repeats<1 or args.warmups<0 or args.query_repeats<1: p.error("Invalid trials")
    cases=fixtures()+[("aligned_gap",[(.02,.02,.02),(.10,.02,.02)]),
                      ("edge_touch",[(.02,.02,.02),(.06,.06,.02)]),
                      ("vertex_touch",[(.02,.02,.02),(.06,.06,.06)])]
    if set(args.cases)-{n for n,_ in cases}: p.error("Unknown case")
    cases=[(n,points) for n,points in cases if n in args.cases]
    output=args.output.resolve(); output.mkdir(parents=True,exist_ok=False)
    binary=args.build_dir.resolve()/"build_halfspace_tree"
    manifest={"binary_sha256":hashlib.sha256(binary.read_bytes()).hexdigest(),"voxel_size":.04,
              "flat_sweep":"y","order":"serial alternating pairs","repeats":args.repeats,"warmups":args.warmups,
              "query_repeats":args.query_repeats,"query_rounds":"median of 7 x 5, original recursive evaluator, no memoization",
              "core_excludes":"I/O, certificates, probes and query benchmark",
              "validation":"Both certify geometry. Direct requires all-strata strict signs; flat retains occupancy zero guard. Probe/reference workloads differ.",
              "tree_nodes":"unique serialized nodes; expanded references also reported (not unique planes)","inputs":[]}
    rows=[]
    for name,points in cases:
        print(name,flush=True); folder=output/name; folder.mkdir(); cloud=folder/"input.xyz"; write_xyz(cloud,points)
        manifest["inputs"].append({"case":name,"raw_points":len(points),"sha256":hashlib.sha256(cloud.read_bytes()).hexdigest()})
        save(output/"manifest.json",manifest)
        def command(mode,out,query=False):
            cmd=[str(binary),str(cloud),str(out),"--pipeline","octree","--tree-source","boundary",
                 "--boundary-expression",mode,"--boundary-occupancy","voxels","--voxel-size",".04",
                 "--boundary-diagnostics","true","--sign-propagation","scalar",
                 "--validation-samples","0" if query else "1000"]
            if mode=="flat": cmd += ["--boundary-sweep","y","--boundary-validation-reference","voxels"]
            else: cmd += ["--boundary-query","recursive","--execution-sharing","none"] # preserve historical evaluator/validator
            if query: cmd += ["--query-benchmark","true"]
            return cmd
        def read(out,mode):
            m=json.loads((out/"benchmark.json").read_text())[0]
            assert json.loads((out/"validation_report.json").read_text())["passed"]
            assert m["passed"] and not m["fallback"] and m["raw_points_classified_free"]==0
            if mode=="direct":
                assert m["strict_sign_certified"] and m["world_sign_errors"]==0
                assert m["artificial_zero_faces"]==m["artificial_zero_edges"]==m["artificial_zero_vertices"]==0
                assert not (out/"free_prisms.json").exists()
            return m
        def pair(a,b,query=False):
            assert (a/"boundary.off").read_bytes()==(b/"boundary.off").read_bytes()
            if query:
                qa,qb=(json.loads((out/"query_samples.json").read_text()) for out in (a,b))
                assert qa["random"]==qb["random"] and qa["boundary"]==qb["boundary"]
        records={m:[] for m in MODES}; queries={m:[] for m in MODES}
        for mode in MODES: (folder/mode).mkdir()
        for trial in range(-args.warmups,args.repeats):
            for mode in (MODES if trial%2==0 else MODES[::-1]):
                out=folder/mode; wall=execute(command(mode,out),out/f"trial_{trial}.log"); m=read(out,mode)
                if trial>=0: records[mode].append({"metrics":m,"process_wall_ms":wall})
            pair(folder/"flat",folder/"direct")
        for trial in range(args.query_repeats):
            for mode in (MODES if trial%2==0 else MODES[::-1]):
                out=folder/"query"/mode; out.mkdir(parents=True,exist_ok=True)
                execute(command(mode,out,True),out/f"query_{trial}.log"); queries[mode].append(read(out,mode))
            pair(folder/"query"/"flat",folder/"query"/"direct",True)
        for mode in MODES:
            save(folder/mode/"trials.json",records[mode]); save(folder/mode/"query_trials.json",queries[mode])
            med=lambda k: statistics.median(v["metrics"][k] for v in records[mode])
            qmed=lambda k: statistics.median(v[k] for v in queries[mode])
            m=records[mode][-1]["metrics"]
            rows.append({"case":name,"mode":mode,"raw_points":len(points),"occupied_voxels":m["occupied_voxels"],
                         "tree_nodes":m["tree_nodes"],"plane_references":m["plane_references"],
                         "expanded_nodes":m["expanded_node_references"],"core_ms":med("core_runtime_ms"),
                         "core_min_ms":min(v["metrics"]["core_runtime_ms"] for v in records[mode]),
                         "core_max_ms":max(v["metrics"]["core_runtime_ms"] for v in records[mode]),
                         "validation_ms":med("validation_ms"),"build_validation_ms":med("total_ms"),
                         "strict_query_us":qmed("strict_free_query_us"),"boundary_query_us":qmed("boundary_strict_free_query_us"),
                         "strict_sign_certified":bool(m["strict_sign_certified"]),"strata_count":m.get("strata_count"),
                         "artificial_zeros_fev":[m["artificial_zero_faces"],m["artificial_zero_edges"],m["artificial_zero_vertices"]],
                         "same_boundary_mesh":True,"same_timed_queries":True,"passed":True})
        save(output/"summary.json",rows)
    report=["# Boundary direct vs flat","",f"h=0.04 m; flat sweep Y; {args.warmups} warmup pairs + {args.repeats} measured pairs.",
            f"Query: {args.query_repeats} paired processes, each median of 7 x 5 rounds, same random/boundary corpora.",
            "Core excludes I/O, validation and query timing. Total includes full validation; validation workloads differ.",
            "Tree nodes count shared nodes once. Plane refs / expanded nodes reflect original recursive evaluator work, WITHOUT memoization.","",
            "| Case | Mode | Nodes | Plane refs | Expanded nodes | Core ms | Build+validation ms | Query us | Boundary us | Artificial zeros F/E/V | Strict signs |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |"]
    for r in rows:
        report.append(f'| {r["case"]} | {r["mode"]} | {r["tree_nodes"]:.0f} | {r["plane_references"]:.0f} | {r["expanded_nodes"]:.0f} | '
                      f'{r["core_ms"]:.4f} | {r["build_validation_ms"]:.3f} | {r["strict_query_us"]:.3f} | {r["boundary_query_us"]:.3f} | '
                      f'{"/".join(str(int(v)) for v in r["artificial_zeros_fev"])} | {r["strict_sign_certified"]} |')
    report += ["","Direct uses certified root<0 alone. Flat strict-query timings include its required zero occupancy guard.",
               "Same occupied voxel boundary, different scalar fields; no cross-backend scalar equality is claimed.",
               "This prototype does not establish a distance field, differentiable CBF or realtime deadline guarantee."]
    (output/"report.md").write_text("\n".join(report)+"\n"); print(output/"report.md",flush=True)

if __name__=="__main__": main()
