#!/usr/bin/env python3
"""Paired probe-cache and packed symbolic-validation benchmark."""
import argparse
import hashlib
import json
from pathlib import Path
import statistics

from benchmark_boundary_tree import REPO, execute, fixtures, save
from benchmark_octree_phase4 import write_xyz

MODES=("scalar_uncached","scalar_cached","packed_cached")
UNCHANGED=("free_direct_tree_raw.json","free_direct_tree_simplified.json","boundary.off",
    "boundary_patches.json","boundary_supports.json","direct_splits.json","occupancy_reference.json",
    "boundary_sign_diagnostics.json","strict_free_space.json","validation_report.json",
    "smooth_distance.json","distance_samples.csv")
CACHE_METRICS=("audit_ms","probe_cache_enabled","probe_cache_entries","probe_cache_hits","probe_cache_misses",
    "probe_evaluations","diagnostic_cache_hits","diagnostic_evaluations","hard_cache_hits","hard_evaluations")

def compare(left,right):
    for name in UNCHANGED:
        assert (left/name).read_bytes()==(right/name).read_bytes(),(name,left,right)
    a=json.loads((left/"distance_validation.json").read_text())
    b=json.loads((right/"distance_validation.json").read_text())
    for audit in (a,b):
        for key in CACHE_METRICS: audit["metrics"].pop(key)
    assert a==b,(left,right,"validation result changed")

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--build-dir",type=Path,default=REPO/"build")
    parser.add_argument("--repeats",type=int,default=5)
    parser.add_argument("--warmups",type=int,default=1)
    parser.add_argument("--validation-samples",type=int,default=1000)
    parser.add_argument("--cases",nargs="+",default=["solid_cuboid","l_voxels","u_voxels","staircase",
        "closed_cavity","two_disconnected","thin_wall","narrow_corridor","camera"])
    args=parser.parse_args(); cases=fixtures()
    if args.repeats<1 or args.warmups<0 or args.validation_samples<0: parser.error("Invalid trial/sample count")
    if set(args.cases)-{name for name,_ in cases}: parser.error("Unknown case")
    output=args.output.resolve(); output.mkdir(parents=True,exist_ok=False)
    binary=args.build_dir.resolve()/"build_halfspace_tree"
    manifest={"binary_sha256":hashlib.sha256(binary.read_bytes()).hexdigest(),"voxel_size_m":.04,
        "repeats":args.repeats,"warmups":args.warmups,"validation_samples":args.validation_samples,
        "modes":MODES,"order":"serial rotating triples","scalar_error_budget_m":.001,
        "work_limits":"original per-occurrence conservative accounting in every mode",
        "coverage":"all volume/face/edge/vertex strata, all original world and distance probe occurrences",
        "comparison":"tree, geometry, probe rows, gradients, rays and validation results unchanged; cache/runtime counters excluded",
        "inputs":[]}
    rows=[]
    for name,points in cases:
        if name not in args.cases: continue
        print(name,flush=True); folder=output/name; folder.mkdir(); cloud=folder/"input.xyz"; write_xyz(cloud,points)
        manifest["inputs"].append({"case":name,"raw_points":len(points),"sha256":hashlib.sha256(cloud.read_bytes()).hexdigest()})
        save(output/"manifest.json",manifest)
        def command(mode,out):
            packed=mode=="packed_cached"; cached=mode!="scalar_uncached"
            return [str(binary),str(cloud),str(out),"--pipeline","octree","--tree-source","boundary",
                "--boundary-expression","direct","--boundary-query","dag","--voxel-size",".04",
                "--validation-samples",str(args.validation_samples),"--distance-field","lse","--lse-error",".001",
                "--execution-sharing","structural","--lse-kernel","binary",
                "--validation-cache","exact" if cached else "none",
                "--sign-propagation","packed" if packed else "scalar"]
        records={mode:[] for mode in MODES}
        for mode in MODES: (folder/mode).mkdir()
        for trial in range(-args.warmups,args.repeats):
            offset=trial%len(MODES)
            for mode in MODES[offset:]+MODES[:offset]:
                out=folder/mode; execute(command(mode,out),out/f"trial_{trial}.log")
                metrics=json.loads((out/"benchmark.json").read_text())[0]
                assert metrics["passed"] and metrics["strict_sign_certified"]
                assert metrics["distance_false_free_count"]==metrics["distance_upper_order_errors"]==metrics["distance_bound_errors"]==0
                assert metrics["distance_work_used"]>0 and metrics["distance_gradient_checks"]>0
                if trial>=0: records[mode].append(metrics)
            for mode in MODES[1:]: compare(folder/MODES[0],folder/mode)
        for mode in MODES:
            save(folder/mode/"trials.json",records[mode]); sample=records[mode][0]
            median=lambda key: statistics.median(record[key] for record in records[mode])
            rows.append({"case":name,"mode":mode,"strata":sample["strata_count"],
                "sign_batches":sample["sign_propagation_batches"],"distance_probes":sample["distance_distance_probes"],
                "probe_cache_hits":sample["distance_probe_cache_hits"],"probe_evaluations":sample["distance_probe_evaluations"],
                "hard_cache_hits":sample["distance_hard_cache_hits"],"hard_evaluations":sample["distance_hard_evaluations"],
                "diagnostic_cache_hits":sample["distance_diagnostic_cache_hits"],
                "certificate_ms":median("certificate_ms"),"world_probe_ms":median("probe_validation_ms"),
                "distance_audit_ms":median("distance_audit_ms"),"full_validation_ms":median("validation_ms"),
                "total_ms":median("total_ms"),"same_outputs":True,"passed":True})
        save(output/"summary.json",rows)
    report=["# Validation acceleration","",
        f"h=.04 m; 1 mm LSE budget; {args.warmups} warmup + {args.repeats} measured rotating triples.",
        "All modes retain the original work limits and every volume, face, edge, vertex and numerical probe occurrence.",
        "Timing is a local measurement; output equality is checked after every triple.","",
        "| Scene | Mode | Strata / batches | Distance probes / evaluated | Probe / hard / diagnostic hits | Certificate ms | Distance audit ms | Full validation ms | Total ms |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in rows:
        report.append(f'| {row["case"]} | {row["mode"]} | {row["strata"]:.0f} / {row["sign_batches"]:.0f} | '
            f'{row["distance_probes"]:.0f} / {row["probe_evaluations"]:.0f} | {row["probe_cache_hits"]:.0f} / '
            f'{row["hard_cache_hits"]:.0f} / {row["diagnostic_cache_hits"]:.0f} | {row["certificate_ms"]:.3f} | '
            f'{row["distance_audit_ms"]:.3f} | {row["full_validation_ms"]:.3f} | {row["total_ms"]:.3f} |')
    report += ["","Cached modes preserve duplicate samples in statistics and exported CSV; only exact-coordinate computations are reused.",
        "Packed propagation evaluates the same ordered MIN/MAX DAG over the same strata, 64 stratum signs per machine word.",
        "Tree, geometry, sampled values/gradients, rays and validation outcomes are identical across all three modes."]
    (output/"report.md").write_text("\n".join(report)+"\n")
    print(output/"report.md",flush=True)

if __name__=="__main__": main()
