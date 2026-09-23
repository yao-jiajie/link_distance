"""Smoke-test the same-process multi-frame benchmark contract."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

repo=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix="multi_frame_reuse_") as temporary:
    out=Path(temporary)/"results"
    subprocess.run([sys.executable,str(repo/"scripts/benchmark_multi_frame_reuse.py"),str(out),
        "--build-dir",sys.argv[1],"--cases","l_voxels","camera","--repeats","3"],check=True,timeout=120)
    rows=json.loads((out/"summary.json").read_text())
    assert len(rows)==2
    for row in rows:
        assert row["passed"] and row["same_occupancy"] and row["same_function"] and row["same_snapshot"]
        assert row["workspaces_reused"] and row["cache_hits"]==3 and row["cache_misses"]==1
        assert row["cold_prepare_ms"]>0 and row["cold_snapshot_build_ms"]>0 and row["hit_prepare_ms"]>0
print("Multi-frame reuse: same-process geometry/program/workspace hits preserve frame checks and values")
