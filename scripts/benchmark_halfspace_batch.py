#!/usr/bin/env python3
"""Same-process scalar and parallel halfspace batch-query benchmark."""
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
    parser.add_argument("--repeats",type=int,default=5)
    parser.add_argument("--sizes",type=int,nargs="+",default=(256,4096,16384))
    parser.add_argument("--workers",type=int,nargs="+",default=(2,4,8))
    args=parser.parse_args()
    if args.repeats<1 or any(v<1 for v in (*args.sizes,*args.workers)): parser.error("Counts must be positive")
    output=args.output.resolve(); output.mkdir(parents=True,exist_ok=False)
    binary=args.build_dir.resolve()/"halfspace_batch_benchmark"
    camera=REPO/"data/realistic_sparse_camera.xyz"; rows=[]
    combinations=[(size,max(args.workers)) for size in args.sizes]
    combinations += [(max(args.sizes),workers) for workers in args.workers if workers!=max(args.workers)]
    for size,workers in combinations:
        run=subprocess.run([str(binary),str(camera),".04",str(size),str(args.repeats),str(workers)],
            check=True,capture_output=True,text=True,timeout=120)
        row=json.loads(run.stdout); row["passed"]=True; rows.append(row)
        (output/f"points_{size}_workers_{workers}.json").write_text(json.dumps(row,indent=2)+"\n")
    (output/"summary.json").write_text(json.dumps(rows,indent=2)+"\n")
    manifest={"binary_sha256":hashlib.sha256(binary.read_bytes()).hexdigest(),
        "input_sha256":hashlib.sha256(camera.read_bytes()).hexdigest(),"voxel_size_m":.04,
        "lse_error_budget_m":.001,"repeats":args.repeats,"sizes":args.sizes,"workers":args.workers,
        "timing":"same process; rotating scalar/batch-1/parallel triples; median microseconds per point; one untimed warmup per mode; reusable workspace thread creation excluded",
        "comparison":"original scalar loop, one-worker batch, and N-worker batch on identical ordered points",
        "correctness":"hard and LSE doubles plus all value-gradient components compared bit-for-bit before timing and after every triple",
        "trial_arrays":"hard/smooth/gradient_trials_us each contains scalar, batch-1, parallel raw per-point times"}
    (output/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    report=["# Halfspace batch and parallel evaluation","",
        "Camera tree, h=.04 m, 1 mm LSE budget, structural sharing and binary LSE kernel.",
        "Times are median microseconds per point in one process. Reusable workspace thread construction is excluded.","",
        "| Points | Workers | Hard scalar / batch-1 / parallel us | Hard x | LSE scalar / batch-1 / parallel us | LSE x | Gradient scalar / batch-1 / parallel us | Gradient x |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in rows:
        report.append(f'| {row["points"]} | {row["workers"]} | {row["hard_scalar_us"]:.3f} / {row["hard_batch_one_us"]:.3f} / {row["hard_parallel_us"]:.3f} | {row["hard_speedup"]:.2f} | '
            f'{row["smooth_scalar_us"]:.3f} / {row["smooth_batch_one_us"]:.3f} / {row["smooth_parallel_us"]:.3f} | {row["smooth_speedup"]:.2f} | '
            f'{row["gradient_scalar_us"]:.3f} / {row["gradient_batch_one_us"]:.3f} / {row["gradient_parallel_us"]:.3f} | {row["gradient_speedup"]:.2f} |')
    report += ["","Batch-1 measures API overhead without parallelism. Parallel output retains input order.",
        "Every hard/LSE value and gradient component matched the scalar path bit-for-bit before timing and after every triple.",
        "Per-trial timings are saved in the JSON records; validation is outside the measured intervals."]
    (output/"report.md").write_text("\n".join(report)+"\n"); print(output/"report.md")

if __name__=="__main__": main()
