#!/usr/bin/env python3
"""Same-process boundary-direct frame snapshot reuse benchmark."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

REPO=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--build-dir",type=Path,default=REPO/"build")
    parser.add_argument("--repeats",type=int,default=21)
    parser.add_argument("--cases",nargs="+",choices=("l_voxels","camera"),default=("l_voxels","camera"))
    args=parser.parse_args()
    if args.repeats<1: parser.error("Repeats must be positive")
    output=args.output.resolve(); output.mkdir(parents=True,exist_ok=False)
    binary=args.build_dir.resolve()/"boundary_direct_frame_benchmark"
    cases={"l_voxels":(REPO/"tests/data/convex_l.xyz",.01),
           "camera":(REPO/"data/realistic_sparse_camera.xyz",.04)}
    rows=[]
    for name in args.cases:
        cloud,h=cases[name]
        run=subprocess.run([str(binary),str(cloud),str(h),str(args.repeats)],check=True,
            capture_output=True,text=True,timeout=120)
        row=json.loads(run.stdout); row["case"]=name; row["passed"]=True; rows.append(row)
    (output/"summary.json").write_text(json.dumps(rows,indent=2)+"\n")
    manifest={"binary_sha256":hashlib.sha256(binary.read_bytes()).hexdigest(),"repeats":args.repeats,
        "process_model":"each case is one process; cold and hit trials execute within that process",
        "frame_change":"raw points move to their occupied voxel centers; sorted occupied keys remain identical",
        "cache_key":"exact voxel-size bits, sorted occupied keys, build/certificate/execution/LSE configuration",
        "coverage":"every hit rebuilds occupied keys and checks every raw point against the closed cached occupancy"}
    (output/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    report=["# Multi-frame geometry and program reuse","",
        f"Release same-process measurements, {args.repeats} cold frames and {args.repeats} warmed cache hits per case.",
        "Cold trials create and certify a new snapshot. Hit trials still voxelize the new raw frame and check all raw points.","",
        "| Scene | Raw points | Voxels | Nodes | Cold prepare ms | Snapshot build ms | Hit prepare ms | Reduction |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in rows:
        report.append(f'| {row["case"]} | {row["raw_points"]} | {row["occupied_voxels"]} | {row["execution_nodes"]} | '
            f'{row["cold_prepare_ms"]:.6f} | {row["cold_snapshot_build_ms"]:.6f} | {row["hit_prepare_ms"]:.6f} | '
            f'{100*row["reduction_fraction"]:.1f}% |')
    report += ["","All hit trials reused the same certified geometry, immutable compiled program, LSE object and workspaces.",
        "Function values and occupied keys were checked in the benchmark; timing excludes process startup and file I/O."]
    (output/"report.md").write_text("\n".join(report)+"\n")
    print(output/"report.md")

if __name__=="__main__": main()
