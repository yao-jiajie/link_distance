import json
from pathlib import Path
import subprocess
import sys
import tempfile

repo=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix="direct_query_benchmark_") as directory:
    out=Path(directory)/"results"
    subprocess.run([sys.executable,str(repo/"scripts/benchmark_boundary_direct_query.py"),str(out),"--build-dir",sys.argv[1],
                    "--repeats","1","--warmups","0","--query-repeats","1","--cases","l_voxels","two_disconnected","aligned_gap"],
                   check=True,timeout=100)
    rows=json.loads((out/"summary.json").read_text()); assert len(rows)==9
    for r in rows:
        assert r["passed"] and r["same_trees_recursive_dag"] and r["same_timed_queries"]
        assert 0<r["core_ms"]<r["build_validation_ms"] and r["query_us"]>0
        if r["mode"]!="flat": assert r["strict_sign_certified"] and r["audit_scalar_checks"]>0 and r["audit_scalar_errors"]==0
        if r["mode"]=="dag": assert r["query_node_visits"]==r["tree_nodes"] and r["workspace_bytes"]==8*r["tree_nodes"]
    assert any(r["query_node_visits"]<r["expanded_nodes"] for r in rows if r["mode"]=="dag")
    assert (out/"report.md").exists()
print("Direct query paired benchmark: identical trees/corpora, scalar audit and certified signs passed")
