import json
from pathlib import Path
import subprocess
import sys
import tempfile

repo=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix="direct_benchmark_") as directory:
    out=Path(directory)/"results"
    subprocess.run([sys.executable,str(repo/"scripts/benchmark_boundary_direct.py"),str(out),"--build-dir",sys.argv[1],
                    "--repeats","1","--warmups","0","--query-repeats","1","--cases","solid_cuboid","l_voxels","aligned_gap"],
                   check=True,timeout=100)
    rows=json.loads((out/"summary.json").read_text()); assert len(rows)==6
    for r in rows:
        assert r["passed"] and r["same_boundary_mesh"] and r["same_timed_queries"]
        assert 0<r["core_ms"]<r["build_validation_ms"] and r["strict_query_us"]>0
        if r["mode"]=="direct": assert r["strict_sign_certified"] and sum(r["artificial_zeros_fev"])==0
    assert any(sum(r["artificial_zeros_fev"]) for r in rows if r["mode"]=="flat")
    assert (out/"report.md").exists()
print("Direct paired benchmark: certified signs, common geometry and identical query corpora passed")
