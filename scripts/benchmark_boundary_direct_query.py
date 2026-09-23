#!/usr/bin/env python3
"""Paired recursive/DAG direct evaluators plus unchanged flat reference."""
import argparse
import hashlib
import json
from pathlib import Path
import statistics
from benchmark_boundary_tree import REPO, execute, fixtures, save
from benchmark_boundary_direct import DEFAULT_CASES
from benchmark_octree_phase4 import write_xyz

MODES=("recursive","dag","flat")

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path); parser.add_argument("--build-dir",type=Path,default=REPO/"build")
    parser.add_argument("--repeats",type=int,default=5); parser.add_argument("--warmups",type=int,default=1)
    parser.add_argument("--query-repeats",type=int,default=3); parser.add_argument("--cases",nargs="+",default=DEFAULT_CASES)
    args=parser.parse_args()
    if args.repeats<1 or args.warmups<0 or args.query_repeats<1: parser.error("Invalid trial count")
    cases=fixtures()+[("aligned_gap",[(.02,.02,.02),(.10,.02,.02)]),
                      ("edge_touch",[(.02,.02,.02),(.06,.06,.02)]),
                      ("vertex_touch",[(.02,.02,.02),(.06,.06,.06)])]
    if set(args.cases)-{n for n,_ in cases}: parser.error("Unknown fixture")
    cases=[(n,p) for n,p in cases if n in args.cases]
    output=args.output.resolve(); output.mkdir(parents=True,exist_ok=False)
    binary=args.build_dir.resolve()/"build_halfspace_tree"
    manifest={"binary_sha256":hashlib.sha256(binary.read_bytes()).hexdigest(),"voxel_size":.04,"modes":MODES,
              "order":"serial rotating triples","repeats":args.repeats,"warmups":args.warmups,"query_repeats":args.query_repeats,
              "query":"7 x 5 median per process, then median across processes; same random/boundary corpora",
              "core":"geometry plus evaluator compilation/workspace, excluding I/O, validation and query benchmark",
              "cross_check":"separate untimed audit: recursive/DAG scalar + signed-zero equality at every validation probe",
              "signs":"mandatory full all-strata certificate for direct; identical-pointer raw/simplified evaluated once in BOTH modes",
              "validation":"same direct probes and certificates; only chosen scalar evaluator changes. Flat retains its existing reference/guard.",
              "inputs":[]}
    rows=[]
    for name,points in cases:
        print(name,flush=True); folder=output/name; folder.mkdir(); cloud=folder/"input.xyz"; write_xyz(cloud,points)
        manifest["inputs"].append({"case":name,"raw_points":len(points),"sha256":hashlib.sha256(cloud.read_bytes()).hexdigest()})
        save(output/"manifest.json",manifest)
        def command(mode,out,query=False,audit=False):
            cmd=[str(binary),str(cloud),str(out),"--pipeline","octree","--tree-source","boundary",
                 "--boundary-expression","flat" if mode=="flat" else "direct","--boundary-occupancy","voxels",
                 "--voxel-size",".04","--boundary-diagnostics","true","--sign-propagation","scalar",
                 "--validation-samples","0" if query else "1000"]
            if mode=="flat": cmd += ["--boundary-sweep","y","--boundary-validation-reference","voxels"]
            else: cmd += ["--boundary-query",mode,"--boundary-query-cross-check","true" if audit else "false",
                          "--execution-sharing","none"] # freeze the original pointer-DAG experiment
            if query: cmd += ["--query-benchmark","true"]
            return cmd
        def read(out,mode):
            m=json.loads((out/"benchmark.json").read_text())[0]
            assert json.loads((out/"validation_report.json").read_text())["passed"] and m["passed"] and not m["fallback"]
            assert m["raw_points_classified_free"]==0
            if mode!="flat":
                assert m["query_backend"]==mode and m["world_sign_errors"]==m["scalar_cross_check_errors"]==0
                assert m["strict_sign_certified"] and m["artificial_zero_faces"]==m["artificial_zero_edges"]==m["artificial_zero_vertices"]==0
                assert m["query_node_visits"]==m["tree_nodes" if mode=="dag" else "expanded_node_references"]
            return m
        def compare(base,query=False):
            a,b=base/"recursive",base/"dag"
            for f in ("free_direct_tree_raw.json","free_direct_tree_simplified.json","boundary_supports.json",
                      "direct_splits.json","boundary.off","boundary_sign_diagnostics.json","occupancy_reference.json"):
                assert (a/f).read_bytes()==(b/f).read_bytes(),(name,f)
            ma,mb=read(a,"recursive"),read(b,"dag")
            for k in ("tested_points","strata_count","tree_nodes","plane_references","occupied_voxels"):
                assert ma[k]==mb[k],(name,k)
            assert (a/"boundary.off").read_bytes()==(base/"flat"/"boundary.off").read_bytes()
            if query:
                q=[json.loads((base/m/"query_samples.json").read_text()) for m in MODES]
                assert all(x["random"]==q[0]["random"] and x["boundary"]==q[0]["boundary"] for x in q)
        # Audit costs are kept OUT of measured records. This is additional
        # recursive validation, not a substitute for geometry/sign certificates.
        audit=folder/"audit"; audit.mkdir(); execute(command("dag",audit,audit=True),audit/"audit.log")
        am=read(audit,"dag"); assert am["scalar_cross_checks"]==am["tested_points"]>0
        records={m:[] for m in MODES}; queries={m:[] for m in MODES}
        for mode in MODES: (folder/mode).mkdir()
        for trial in range(-args.warmups,args.repeats):
            offset=trial%3
            for mode in MODES[offset:]+MODES[:offset]:
                out=folder/mode; wall=execute(command(mode,out),out/f"trial_{trial}.log"); m=read(out,mode)
                if trial>=0: records[mode].append({"metrics":m,"process_wall_ms":wall})
            compare(folder)
        for trial in range(args.query_repeats):
            offset=trial%3
            for mode in MODES[offset:]+MODES[:offset]:
                out=folder/"query"/mode; out.mkdir(parents=True,exist_ok=True)
                execute(command(mode,out,query=True),out/f"query_{trial}.log"); queries[mode].append(read(out,mode))
            compare(folder/"query",True)
        for mode in MODES:
            save(folder/mode/"trials.json",records[mode]); save(folder/mode/"query_trials.json",queries[mode])
            m=records[mode][-1]["metrics"]
            med=lambda k: statistics.median(v["metrics"].get(k,0) for v in records[mode])
            qmed=lambda k: statistics.median(v[k] for v in queries[mode])
            rows.append({"case":name,"mode":mode,"raw_points":len(points),"occupied_voxels":m["occupied_voxels"],
                "tree_nodes":m["tree_nodes"],"expanded_nodes":m["expanded_node_references"],"plane_references":m["plane_references"],
                "query_node_visits":m.get("query_node_visits",m["expanded_node_references"]),
                "core_ms":med("core_runtime_ms"),"compile_ms":med("evaluator_compile_ms"),
                "core_min_ms":min(v["metrics"]["core_runtime_ms"] for v in records[mode]),
                "core_max_ms":max(v["metrics"]["core_runtime_ms"] for v in records[mode]),
                "certificate_ms":med("certificate_ms"),"probe_ms":med("probe_validation_ms"),
                "build_validation_ms":med("total_ms"),"query_us":qmed("strict_free_query_us"),
                "boundary_query_us":qmed("boundary_strict_free_query_us"),"program_bytes":m.get("compiled_program_bytes",0),
                "workspace_bytes":m.get("query_workspace_bytes",0),"strict_sign_certified":bool(m["strict_sign_certified"]),
                "artificial_zeros_fev":[m["artificial_zero_faces"],m["artificial_zero_edges"],m["artificial_zero_vertices"]],
                "audit_scalar_checks":am["scalar_cross_checks"] if mode!="flat" else None,
                "audit_scalar_errors":0 if mode!="flat" else None,"same_trees_recursive_dag":True,"same_timed_queries":True,"passed":True})
        save(output/"summary.json",rows)
    report=["# Direct query optimization: unchanged tree, recursive vs compiled DAG","",
            f"h=.04 m; {args.warmups} warmup + {args.repeats} measured rotating triples; {args.query_repeats} query processes (7 x 5 rounds each).",
            "Compilation and reusable workspace preparation are INCLUDED in core. No per-query allocation or cross-point cache.",
            "Both direct modes run the same all-strata certificate and geometry probes; separate audit verifies full scalar and signed-zero equality.",
            "Direct raw/simplified trees, splits, coefficients, certificates and occupancy are byte-identical across evaluators.","",
            "| Case | Mode | Nodes | Visited nodes/query | Core ms | Compile ms | Build+validation ms | Query us | Boundary us | Zeros F/E/V |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |"]
    for r in rows:
        report.append(f'| {r["case"]} | {r["mode"]} | {r["tree_nodes"]:.0f} | {r["query_node_visits"]:.0f} | '
                      f'{r["core_ms"]:.4f} | {r["compile_ms"]:.4f} | {r["build_validation_ms"]:.3f} | '
                      f'{r["query_us"]:.3f} | {r["boundary_query_us"]:.3f} | '
                      f'{"/".join(str(int(v)) for v in r["artificial_zeros_fev"])} |')
    report += ["","Flat geometry is identical but its scalar field differs and strict querying retains the original occupancy guard.",
               "Expanded references are unchanged; DAG visits each stored node once. This is execution optimization, not tree minimization.",
               "Both current direct modes reuse the raw/simplified symbolic result only when their expression pointers are identical.",
               "Core and build+validation exclude input/output I/O and optional query benchmark. No realtime deadline or CBF guarantee."]
    (output/"report.md").write_text("\n".join(report)+"\n"); print(output/"report.md",flush=True)

if __name__=="__main__": main()
