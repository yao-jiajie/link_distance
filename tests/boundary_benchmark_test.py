"""Exercise five-way benchmarking, including best-axis costs and work fallback."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

repo = Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix="boundary_xyz_benchmark_") as temporary:
    output = Path(temporary)/"results"
    subprocess.run([sys.executable,str(repo/"scripts/benchmark_boundary_tree.py"),str(output),
                    "--build-dir",sys.argv[1],"--warmups","0","--repeats","1",
                    "--cases","single","u_voxels","sparse_disconnected"],check=True,timeout=100)
    rows = json.loads((output/"summary.json").read_text())
    modes = {"volume","sweep-x","sweep-y","sweep-z","best-axis"}
    assert len(rows) == 15 and all(r["passed"] for r in rows)
    for name in {r["case"] for r in rows}:
        group = [r for r in rows if r["case"]==name]
        assert {r["mode"] for r in group} == modes
        for row in group:
            assert row["core_ms"] <= row["build_validation_ms"]
            if row["mode"]=="volume":
                continue
            assert row["candidate_count"] == (3 if row["mode"]=="best-axis" else 1)
            assert row["core_ms"] >= row["all_candidates_ms"]+row["selection_ms"]
            assert row["all_candidates_ms"] >= row["selected_candidate_build_ms"]
            if name=="sparse_disconnected":
                assert row["fallback"] and row["fallback_reason"]=="compressed_cell_limit"
    assert len(json.loads((output/"manifest.json").read_text())["inputs"]) == 3
print("Five-way benchmark, selection and complete-cost accounting passed")
