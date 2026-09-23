#!/usr/bin/env python3
"""Fixed-axis volume/flat/local comparison; incomplete strict-sign results retained."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import statistics

from benchmark_boundary_tree import REPO, execute, fixtures, save
from benchmark_octree_phase4 import write_xyz

MODES=("volume","flat","local")
def cases():
    selected={"single","adjacent_cuboids","l_voxels","u_voxels","staircase","closed_cavity",
              "two_disconnected","camera","sparse_disconnected"}
    out=[(n,p) for n,p in fixtures() if n in selected]
    cavity=[]
    for x,y,z in itertools.product(range(-2,2),range(-1,3),range(-1,2)):
        if not(z==0 and ((x==-1 and y==0) or (x==0 and y in (0,1)))):
            cavity.append(tuple(.04*(a+.5) for a in (x,y,z)))
    out += [("l_free_cavity",cavity),("aligned_gap",[(.02,.02,.02),(.10,.02,.02)])]
    return out

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--build-dir",type=Path,default=REPO/"build")
    parser.add_argument("--axis",choices=("x","y","z"),default="y")
    parser.add_argument("--repeats",type=int,default=5)
    parser.add_argument("--warmups",type=int,default=1)
    parser.add_argument("--cases",nargs="+")
    args=parser.parse_args()
    if args.repeats<1 or args.warmups<0: parser.error("Invalid trial count")
    inputs=cases()
    if args.cases:
        if set(args.cases)-{n for n,_ in inputs}: parser.error("Unknown case")
        inputs=[(n,p) for n,p in inputs if n in args.cases]
    output=args.output.resolve(); output.mkdir(parents=True,exist_ok=False)
    binary=args.build_dir.resolve()/"build_halfspace_tree"
    manifest={"binary_sha256":hashlib.sha256(binary.read_bytes()).hexdigest(),"axis":args.axis,"voxel_size":.04,
              "warmups":args.warmups,"repeats":args.repeats,"modes":MODES,"inputs":[],
              "core_excludes":"sign certificates, probe validation, input/output I/O, query benchmark"}
    rows=[]
    for name,points in inputs:
        print(name,flush=True)
        folder=output/name; folder.mkdir(); cloud=folder/"input.xyz"; write_xyz(cloud,points)
        manifest["inputs"].append({"name":name,"sha256":hashlib.sha256(cloud.read_bytes()).hexdigest()})
        save(output/"manifest.json",manifest)
        def command(mode,out,query=False):
            cmd=[str(binary),str(cloud),str(out),"--pipeline","octree","--voxel-size",".04",
                 "--tree-source","volume" if mode=="volume" else "boundary","--validation-samples","0" if query else "1000"]
            if mode!="volume":
                cmd += ["--boundary-sweep",args.axis,"--boundary-expression",mode,"--boundary-diagnostics","true",
                        "--sign-propagation","scalar"]
            if query: cmd += ["--query-benchmark","true"]
            return cmd
        records={m:[] for m in MODES}
        for mode in MODES: (folder/mode).mkdir()
        for trial in range(-args.warmups,args.repeats):
            k=trial%3
            for mode in MODES[k:]+MODES[:k]:
                out=folder/mode; wall=execute(command(mode,out),out/f"trial_{trial}.log")
                metrics=json.loads((out/"benchmark.json").read_text())[0]
                assert json.loads((out/"validation_report.json").read_text())["passed"]
                if mode!="volume" and not metrics["fallback"]:
                    assert metrics["sign_certificate_applied"] and metrics["unsafe_free_count"]==metrics["true_boundary_sign_errors"]==0
                if trial>=0: records[mode].append({"metrics":metrics,"process_wall_ms":wall})
        # The flat reference exported by local must be identical to standalone flat.
        if not records["local"][-1]["metrics"]["fallback"]:
            for suffix in ("raw","simplified"):
                f=f"free_closure_tree_{suffix}.json"
                assert (folder/"flat"/f).read_bytes()==(folder/"local"/f).read_bytes()
        query={}
        for mode in ("flat","local"):
            out=folder/"query"/mode; out.mkdir(parents=True)
            execute(command(mode,out,True),out/"query.log")
            query[mode]=json.loads((out/"benchmark.json").read_text())[0]
        for mode,runs in records.items():
            save(folder/mode/"trials.json",runs)
            m=runs[-1]["metrics"]
            median=lambda key: statistics.median(r["metrics"][key] for r in runs)
            row={"case":name,"mode":mode,"tree_nodes":m["tree_nodes"],"constructed_nodes":m.get("constructed_tree_nodes"),
                 "expanded_node_references":m.get("expanded_node_references"),"plane_references":m.get("plane_references",m.get("planes_after_reduce")),
                 "core_ms":median("core_runtime_ms"),"build_validation_ms":median("total_ms"),
                 "strict_query_us":query[mode]["strict_free_query_us"] if mode!="volume" else
                    statistics.median(q["baseline_strict_free_query_us"] for q in query.values()),
                 "strict_sign_certified":bool(m.get("strict_sign_certified",False)) if mode!="volume" else None,
                 "artificial_zero_faces":m.get("artificial_zero_faces"),"artificial_zero_edges":m.get("artificial_zero_edges"),
                 "artificial_zero_vertices":m.get("artificial_zero_vertices"),"sign_certificate_applied":bool(m.get("sign_certificate_applied",False)),
                 "logical_merges":m.get("logical_merge_count",0),"adjacency_faces":m.get("adjacency_face_count",0),
                 "factoring_failed":m.get("factoring_failed_count",0),"local_simplify_ms":median("local_simplify_ms") if mode=="local" else None,
                 "fallback":bool(m.get("fallback",False)),"classification_passed":True}
            rows.append(row)
        save(output/"summary.json",rows)
    report=["# Fixed-axis boundary local prototype","",f"axis={args.axis}; h=0.04 m; {args.warmups} warmups, {args.repeats} measured trials; serial rotated order.","",
            "Strict-sign success and guarded classification success are DIFFERENT. Incomplete candidates retain the occupancy guard.",
            "Construction already determines the new scalar; subsequent simplification must preserve it and does not repair seams.","",
            "| Case | Mode | Nodes | Expanded refs | Core ms | Build+validation ms | Strict query us | Free zero F/E/V | Strict sign |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | --- | --- |"]
    for r in rows:
        zeros="-" if not r["sign_certificate_applied"] else f'{r["artificial_zero_faces"]}/{r["artificial_zero_edges"]}/{r["artificial_zero_vertices"]}'
        status="fallback" if r["fallback"] else "not assessed" if r["mode"]=="volume" else str(r["strict_sign_certified"])
        report.append(f'| {r["case"]} | {r["mode"]} | {r["tree_nodes"]} | {r["expanded_node_references"] or "-"} | '
                      f'{r["core_ms"]:.3f} | {r["build_validation_ms"]:.3f} | {r["strict_query_us"]:.3f} | {zeros} | {status} |')
    report += ["","Core includes flat reference construction, local adjacency/templates, and existing structural simplification when local is enabled.",
               "Full sign certificates and all probes are charged to validation. These timings exclude file I/O and separate query benchmarking.",
               "Zeros count relative-open strata in the compressed plane arrangement, not original patches. Unbounded intervals are included.",
               "Constructed unique nodes, simplified unique nodes and recursively expanded references can differ substantially.",
               "Query timings use the existing recursive evaluator, including any required zero guard; paired baseline estimates are retained.",
               "No BSP, R-function, smoothing, Nef or global Boolean minimization was introduced."]
    (output/"report.md").write_text("\n".join(report)+"\n")
    print(output/"report.md",flush=True)

if __name__=="__main__": main()
