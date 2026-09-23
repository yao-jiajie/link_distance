import json
from pathlib import Path
import subprocess
import sys
import tempfile

repo=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix="smooth_distance_benchmark_") as tmp:
    out=Path(tmp)/"results"
    subprocess.run([sys.executable,str(repo/"scripts/benchmark_smooth_distance.py"),str(out),"--build-dir",sys.argv[1],
        "--errors",".001",".01","--warmups","0","--repeats","1","--cases","solid_cuboid","l_voxels"],check=True,timeout=100)
    rows=json.loads((out/"summary.json").read_text()); assert len(rows)==4
    for r in rows:
        assert r["passed"] and r["same_geometry"] and r["same_timed_queries"] and r["false_free_count"]==0
        assert r["smoothing_error_max_m"]<=r["error_budget_m"]+1e-12
        assert r["smooth_gradient_query_us"]>0 and r["gradient_checks"]>0
        assert r["core_ms"]<r["build_validation_ms"] and r["hard_sdf_error_max_m"]>0
print("Smooth benchmark: separate smoothing/SDF errors, intact exact tree, matching queries and conservative probes passed")
