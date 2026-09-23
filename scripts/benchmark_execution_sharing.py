#!/usr/bin/env python3
"""Paired execution-only CSE or binary LSE kernel benchmark; unchanged function."""
import argparse
import hashlib
import json
from pathlib import Path
import statistics

from benchmark_boundary_tree import REPO, execute, fixtures, save
from benchmark_octree_phase4 import write_xyz

MODES=("none","structural")
UNCHANGED=("free_direct_tree_raw.json","free_direct_tree_simplified.json","boundary.off",
    "boundary_patches.json","boundary_supports.json","direct_splits.json","occupancy_reference.json",
    "boundary_sign_diagnostics.json","strict_free_space.json","validation_report.json",
    "smooth_distance.json","distance_samples.csv")

def compare(a,b,ignore_work=True):
    for name in UNCHANGED:
        assert (a/name).read_bytes()==(b/name).read_bytes(),(name,a,b)
    left=json.loads((a/"distance_validation.json").read_text())
    right=json.loads((b/"distance_validation.json").read_text())
    for audit in (left,right):
        for key in (("audit_ms","work_used") if ignore_work else ("audit_ms",)): audit["metrics"].pop(key)
    assert left==right,(a,b,"distance audit changed")

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--build-dir",type=Path,default=REPO/"build")
    parser.add_argument("--comparison",choices=("sharing","kernel"),default="sharing")
    parser.add_argument("--repeats",type=int,default=3)
    parser.add_argument("--warmups",type=int,default=1)
    parser.add_argument("--query-repeats",type=int,default=3)
    parser.add_argument("--errors",type=float,nargs="+",default=[.001])
    parser.add_argument("--cases",nargs="+",default=["solid_cuboid","l_voxels","u_voxels","staircase",
        "closed_cavity","two_disconnected","thin_wall","narrow_corridor","camera"])
    args=parser.parse_args(); cases=fixtures()
    sharing=args.comparison=="sharing"; modes=MODES if sharing else ("generic","binary")
    if args.repeats<1 or args.warmups<0 or args.query_repeats<1: parser.error("Invalid repeat count")
    if any(not 0<e<1e6 for e in args.errors) or len(set(args.errors))!=len(args.errors): parser.error("Invalid/duplicate errors")
    if set(args.cases)-{n for n,_ in cases}: parser.error("Unknown case")
    output=args.output.resolve(); output.mkdir(parents=True,exist_ok=False)
    binary=args.build_dir.resolve()/"build_halfspace_tree"
    manifest={"binary_sha256":hashlib.sha256(binary.read_bytes()).hexdigest(),"voxel_size_m":.04,
        "repeats":args.repeats,"warmups":args.warmups,"query_repeats":args.query_repeats,
        "errors_m":args.errors,"modes":modes,"comparison_kind":args.comparison,"order":"serial rotating pairs",
        "fixed_execution_sharing":None if sharing else "structural","fixed_lse_kernel":"generic" if sharing else None,
        "query":"separate processes; median across processes of 7 x 5 medians, same corpora",
        "timing":"core includes sharing compilation; full validation unchanged; I/O and query benchmark excluded",
        "comparison":"original tree, geometry, function parameters, CSV values/gradients byte-identical; audit unchanged except time/work",
        "inputs":[]}
    rows=[]
    for name,points in cases:
        if name not in args.cases: continue
        folder=output/name; folder.mkdir(); cloud=folder/"input.xyz"; write_xyz(cloud,points)
        manifest["inputs"].append({"case":name,"sha256":hashlib.sha256(cloud.read_bytes()).hexdigest(),"raw_points":len(points)})
        save(output/"manifest.json",manifest)
        for error in args.errors:
            print(name,error,flush=True); group=folder/("error_"+str(error)); group.mkdir()
            def command(mode,out,query=False):
                cmd=[str(binary),str(cloud),str(out),"--pipeline","octree","--tree-source","boundary",
                    "--boundary-expression","direct","--boundary-query","dag","--voxel-size",".04",
                    "--distance-field","lse","--lse-error",str(error),"--execution-sharing",mode if sharing else "structural",
                    "--lse-kernel","generic" if sharing else mode,"--validation-cache","none",
                    "--sign-propagation","scalar"]
                if query: cmd += ["--query-benchmark","true"]
                return cmd
            records={m:[] for m in modes}; queries={m:[] for m in modes}
            for mode in modes: (group/mode).mkdir()
            def read(out):
                m=json.loads((out/"benchmark.json").read_text())[0]
                assert m["passed"] and m["strict_sign_certified"] and m["distance_false_free_count"]==0
                assert m["distance_upper_order_errors"]==m["distance_bound_errors"]==0
                assert m["distance_gradient_checks"]>0
                return m
            for trial in range(-args.warmups,args.repeats):
                for mode in modes[trial%2:]+modes[:trial%2]:
                    out=group/mode; execute(command(mode,out),out/f"trial_{trial}.log")
                    m=read(out)
                    if trial>=0: records[mode].append(m)
                compare(group/modes[0],group/modes[1],sharing)
            for mode in modes:
                save(group/mode/"trials.json",records[mode]); (group/mode/"query").mkdir()
            for trial in range(args.query_repeats):
                for mode in modes[trial%2:]+modes[:trial%2]:
                    out=group/mode/"query"; execute(command(mode,out,True),out/f"query_{trial}.log")
                    queries[mode].append(read(out))
                compare(group/modes[0]/"query",group/modes[1]/"query",sharing)
                assert (group/modes[0]/"query/query_samples.json").read_bytes()==(group/modes[1]/"query/query_samples.json").read_bytes()
            for mode in modes:
                save(group/mode/"query_trials.json",queries[mode])
                m=records[mode][0]; median=lambda key: statistics.median(r[key] for r in records[mode])
                qmedian=lambda key: statistics.median(r[key] for r in queries[mode])
                rows.append({"case":name,"error_budget_m":error,"mode":mode,"source_nodes":m["tree_nodes"],
                    "comparison_kind":args.comparison,"binary_nodes":m["smooth_binary_nodes"],
                    "generic_nodes":m["smooth_generic_nodes"],"binary_kernel":m["smooth_binary_kernel"],
                    "execution_nodes":m["smooth_execution_nodes"],"shared_nodes":m["smooth_shared_nodes"],
                    "exp_calls":m["smooth_exp_calls_per_query"],"log_calls":m["smooth_log_calls_per_query"],
                    "core_ms":median("core_runtime_ms"),"hard_compile_ms":median("evaluator_compile_ms"),
                    "smooth_compile_ms":median("smooth_compile_ms"),"total_ms":median("total_ms"),
                    "smooth_audit_ms":median("sdf_validation_ms"),"hard_query_us":qmedian("strict_free_query_us"),
                    "value_query_us":qmedian("smooth_value_query_us"),"gradient_query_us":qmedian("smooth_gradient_query_us"),
                    "program_bytes":m["smooth_program_bytes"],"gradient_workspace_bytes":m["smooth_gradient_workspace_bytes"],
                    "same_tree_and_geometry":True,"same_values_and_gradients":True,"same_timed_queries":True,"passed":True})
            save(output/"summary.json",rows)
    report=["# "+("Execution-only structural sharing" if sharing else "Binary LSE kernel, structural sharing fixed"),"",
        "No child removal/reordering, no Boolean simplification, no geometry/LSE formula change.",
        "Paired runs, full validation retained, compilation cost included. Performance is local measurement, not a guarantee.","",
        "| Scene | Budget mm | Mode | Source / execution nodes | Binary nodes | exp / log calls | Core ms | Full build+validation ms | Hard us | LSE value us | Value+gradient us |",
        "| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for r in rows:
        report.append(f'| {r["case"]} | {r["error_budget_m"]*1000:g} | {r["mode"]} | {r["source_nodes"]:g} / {r["execution_nodes"]:g} | '
            f'{r["binary_nodes"]:g} | {r["exp_calls"]:g} / {r["log_calls"]:g} | {r["core_ms"]:.3f} | {r["total_ms"]:.3f} | '
            f'{r["hard_query_us"]:.3f} | {r["value_query_us"]:.3f} | {r["gradient_query_us"]:.3f} |')
    report += ["","All original bundles and sampled value/gradient CSVs are byte-identical between modes.",
        "Numerical distance/gradient/ray audits match except runtime/work estimates. Exact all-strata validation still runs.",
        ("Execution buffer capacities are retained; workspace shrinks when CSE merges, but program capacity need not shrink."
         if sharing else "Both kernels have identical nodes, edges, exp/log counts, work estimates and workspace sizes."),
        "This optimization does not improve Euclidean SDF accuracy or change the smoothed zero set."]
    (output/"report.md").write_text("\n".join(report)+"\n")
    print(output/"report.md",flush=True)

if __name__=="__main__": main()
