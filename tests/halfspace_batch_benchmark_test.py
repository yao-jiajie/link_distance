"""Smoke-test batch benchmark correctness and metrics."""
import json
from pathlib import Path
import subprocess
import statistics
import sys

binary=Path(sys.argv[1]).resolve(); repo=Path(__file__).resolve().parents[1]
run=subprocess.run([str(binary),str(repo/"data/realistic_sparse_camera.xyz"),".04","512","1","2"],
    check=True,capture_output=True,text=True,timeout=120)
row=json.loads(run.stdout)
assert row["points"]==512 and row["workers"]==2 and row["execution_nodes"]==429
assert row["exact_hard"] and row["exact_smooth"] and row["exact_gradient"]
for kind in ("hard","smooth","gradient"):
    assert row[kind+"_scalar_us"]>0 and row[kind+"_batch_one_us"]>0 and row[kind+"_parallel_us"]>0
    for mode,trials in zip(("scalar","batch_one","parallel"),row[kind+"_trials_us"]):
        assert len(trials)==row["repeats"] and statistics.median(trials)==row[kind+"_"+mode+"_us"]
assert row["gradient_batch_workspace_bytes"]>row["smooth_batch_workspace_bytes"]>=row["hard_batch_workspace_bytes"]
print("Halfspace batch benchmark: ordered hard/LSE/gradient results match scalar execution")
