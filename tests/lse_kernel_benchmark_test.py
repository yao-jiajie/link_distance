"""Whole-scene binary/generic comparison with the SAME structural CSE program."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

repo=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix="lse_kernel_") as tmp:
    out=Path(tmp)/"results"
    subprocess.run([sys.executable,str(repo/"scripts/benchmark_execution_sharing.py"),str(out),
        "--comparison","kernel","--build-dir",sys.argv[1],"--cases","solid_cuboid","camera",
        "--warmups","0","--repeats","1","--query-repeats","1"],check=True,timeout=100)
    rows=json.loads((out/"summary.json").read_text()); assert len(rows)==4
    for case in ("solid_cuboid","camera"):
        a,b=[r for r in rows if r["case"]==case]
        assert a["mode"]=="generic" and b["mode"]=="binary"
        for r in (a,b):
            assert r["passed"] and r["same_tree_and_geometry"] and r["same_values_and_gradients"] and r["same_timed_queries"]
            assert r["comparison_kind"]=="kernel" and r["value_query_us"]>0 and r["gradient_query_us"]>0
        assert a["binary_nodes"]==a["binary_kernel"]==0 and b["binary_nodes"]>0 and b["binary_kernel"]==1
        for k in ("source_nodes","execution_nodes","shared_nodes","exp_calls","log_calls","gradient_workspace_bytes"):
            assert a[k]==b[k],k
        assert b["binary_nodes"]+b["generic_nodes"]==b["log_calls"]
        if case=="camera": assert b["execution_nodes"]==429 and b["binary_nodes"]==297 and b["generic_nodes"]==71
print("Binary LSE: paired output/gradient/diagnostic identity, same execution program, query corpora and work counts passed")
