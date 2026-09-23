#!/usr/bin/env python3
"""Measure smoothing vs geometric SDF error separately, using the unchanged direct DAG."""
import argparse
import hashlib
import json
from pathlib import Path
import statistics
from benchmark_boundary_tree import REPO, execute, fixtures, save
from benchmark_octree_phase4 import write_xyz

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("output",type=Path); p.add_argument("--build-dir",type=Path,default=REPO/"build")
    p.add_argument("--errors",nargs="+",type=float,default=[.001,.005,.01])
    p.add_argument("--repeats",type=int,default=3); p.add_argument("--warmups",type=int,default=1)
    p.add_argument("--cases",nargs="+",default=["solid_cuboid","l_voxels","u_voxels","staircase","closed_cavity","thin_wall","narrow_corridor","two_disconnected","camera"])
    args=p.parse_args()
    if args.repeats<1 or args.warmups<0 or any(not 0<e<1e6 for e in args.errors) or len(set(args.errors))!=len(args.errors): p.error("Invalid trials/errors")
    cases=fixtures()
    if set(args.cases)-{n for n,_ in cases}: p.error("Unknown scene")
    output=args.output.resolve(); output.mkdir(parents=True,exist_ok=False); binary=args.build_dir.resolve()/"build_halfspace_tree"
    manifest={"binary_sha256":hashlib.sha256(binary.read_bytes()).hexdigest(),"voxel_size":.04,"scalar_error_budgets_m":args.errors,
        "warmups":args.warmups,"repeats":args.repeats,"query_rounds":"one separate process per mode; median of 7 x 5 rounds",
        "errors":"seeded sampling, finite rectangular boundary reference; no global SDF or Hausdorff guarantee",
        "timing":"core includes smooth compiler/workspaces; full validation and reference queries reported separately; I/O excluded", "inputs":[]}
    rows=[]
    for name,points in cases:
        if name not in args.cases: continue
        print(name,flush=True); folder=output/name; folder.mkdir(); cloud=folder/"input.xyz"; write_xyz(cloud,points)
        manifest["inputs"].append({"case":name,"raw_points":len(points),"sha256":hashlib.sha256(cloud.read_bytes()).hexdigest()})
        save(output/"manifest.json",manifest)
        def command(out,error=None,query=False):
            cmd=[str(binary),str(cloud),str(out),"--pipeline","octree","--tree-source","boundary","--boundary-expression","direct",
                 "--voxel-size",".04","--validation-samples","1000","--boundary-query","dag",
                 "--execution-sharing","none","--sign-propagation","scalar"] # freeze the original smoothing benchmark
            if error is not None: cmd += ["--distance-field","lse","--lse-error",str(error),"--lse-kernel","generic",
                                           "--validation-cache","none"]
            if query: cmd += ["--query-benchmark","true"]
            return cmd
        exact=folder/"exact"; exact.mkdir(); execute(command(exact,query=True),exact/"run.log")
        exact_m=json.loads((exact/"benchmark.json").read_text())[0]
        for error in args.errors:
            out=folder/("error_"+str(error)); out.mkdir(); records=[]
            for trial in range(-args.warmups,args.repeats):
                execute(command(out,error),out/f"trial_{trial}.log"); m=json.loads((out/"benchmark.json").read_text())[0]
                if trial>=0: records.append(m)
            save(out/"trials.json",records)
            for f in ("free_direct_tree_raw.json","free_direct_tree_simplified.json","boundary.off","strict_free_space.json","boundary_sign_diagnostics.json"):
                assert (out/f).read_bytes()==(exact/f).read_bytes(),(name,error,f)
            validation=json.loads((out/"distance_validation.json").read_text()); metrics=validation["metrics"]
            assert validation["passed"] and metrics["false_free_count"]==metrics["upper_order_errors"]==metrics["bound_errors"]==0
            assert metrics["gradient_checks"]>0
            query=out/"query"; query.mkdir(); execute(command(query,error,True),query/"run.log")
            q=json.loads((query/"benchmark.json").read_text())[0]
            assert json.loads((query/"query_samples.json").read_text())==json.loads((exact/"query_samples.json").read_text())
            med=lambda k: statistics.median(r[k] for r in records)
            rows.append({"case":name,"raw_points":len(points),"error_budget_m":error,"beta":m["lse_beta"],"tree_nodes":m["tree_nodes"],
                "core_ms":med("core_runtime_ms"),"smooth_compile_ms":med("smooth_compile_ms"),"build_validation_ms":med("total_ms"),
                "smooth_validation_ms":med("sdf_validation_ms"),"hard_query_us":exact_m["strict_free_query_us"],
                "smooth_value_query_us":q["smooth_value_query_us"],"smooth_gradient_query_us":q["smooth_gradient_query_us"],
                "reference_query_us":q["reference_sdf_query_us"],"same_geometry":True,"same_timed_queries":True,"passed":True,**metrics})
            save(output/"summary.json",rows)
    report=["# Conservative LSE distance proxy: smoothing and SDF errors","",f"h=.04 m; scalar error budgets {args.errors}; {args.warmups} warmup + {args.repeats} measured runs.",
        "Separate query process: median of 7 x 5 rounds, identical query corpora. These are local machine samples, not realtime guarantees.",
        "The reference is Euclidean distance to FINITE exposed patches, signed by original closed voxel occupancy.",
        "LSE error bound controls |smooth-hard|, NOT |smooth-true SDF|. Strict exact-tree bundle remains unchanged.","",
        "| Scene | Budget mm | Beta 1/m | Smooth error max mm | Hard SDF error max mm | Smooth SDF error max mm | Smooth SDF RMSE mm | Zero ray offset max mm | Value us | Value+grad us | Core ms | Build+validation ms |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for r in rows:
        report.append(f'| {r["case"]} | {1000*r["error_budget_m"]:g} | {r["beta"]:.1f} | {1000*r["smoothing_error_max_m"]:.3f} | '
            f'{1000*r["hard_sdf_error_max_m"]:.3f} | {1000*r["sdf_error_max_m"]:.3f} | {1000*r["sdf_error_rmse_m"]:.3f} | '
            f'{1000*r["zero_ray_max_offset_m"]:.3f} | {r["smooth_value_query_us"]:.3f} | {r["smooth_gradient_query_us"]:.3f} | '
            f'{r["core_ms"]:.3f} | {r["build_validation_ms"]:.3f} |')
    report += ["","passed means numerical bounds/gradients/false-free checks passed, NOT certified SDF accuracy.",
        "Zero-ray offsets are sampled outward brackets within 2 voxels, not global Hausdorff distance; blocked/unbracketed rays are retained in JSON.",
        "Conservativeness is proved for the real arithmetic formula; no interval-certified libm, distance reinitialization or CBF controller is claimed."]
    (output/"report.md").write_text("\n".join(report)+"\n"); print(output/"report.md",flush=True)

if __name__=="__main__": main()
