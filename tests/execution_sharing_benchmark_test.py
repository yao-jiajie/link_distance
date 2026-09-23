"""End-to-end CSE comparison, including the scene with real structural duplicates."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

repo=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix="execution_sharing_") as tmp:
    out=Path(tmp)/"results"
    subprocess.run([sys.executable,str(repo/"scripts/benchmark_execution_sharing.py"),str(out),
        "--build-dir",sys.argv[1],"--cases","l_voxels","camera","--warmups","0",
        "--repeats","1","--query-repeats","1"],check=True,timeout=100)
    rows=json.loads((out/"summary.json").read_text()); assert len(rows)==4
    for r in rows:
        assert r["passed"] and r["same_tree_and_geometry"] and r["same_values_and_gradients"] and r["same_timed_queries"]
        assert r["execution_nodes"]+r["shared_nodes"]==r["source_nodes"]
        assert r["gradient_workspace_bytes"]==32*r["execution_nodes"]
        assert r["core_ms"]>r["smooth_compile_ms"]>0 and r["gradient_query_us"]>0
        if r["mode"]=="none": assert r["shared_nodes"]==0
        if r["case"]=="camera" and r["mode"]=="structural":
            assert r["source_nodes"]==524 and r["execution_nodes"]==429
            assert r["exp_calls"]==439 and r["log_calls"]==368
    # Independent signature count from the EXPORTED, unmodified DAG.
    tree=json.loads((out/"camera/error_0.001/none/free_direct_tree_raw.json").read_text())
    mapped={}; signatures={}
    for n in tree["nodes"]:
        key=("leaf",n["leaf_id"]) if n["type"]=="leaf" else (n["type"],tuple(mapped[c] for c in n["children"]))
        if key not in signatures: signatures[key]=len(signatures)
        mapped[n["id"]]=signatures[key]
    assert len(signatures)==429
print("Execution CSE: original exports/gradients/bounds unchanged, 524->429 independent signature count, paired benchmarks passed")
