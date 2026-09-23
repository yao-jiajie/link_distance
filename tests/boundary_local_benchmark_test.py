"""Smoke test for the fixed-axis local benchmark, retaining negative results."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

repo=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix="boundary_local_benchmark_") as directory:
    output=Path(directory)/"results"
    subprocess.run([sys.executable,str(repo/"scripts/benchmark_boundary_local.py"),str(output),
                    "--build-dir",sys.argv[1],"--repeats","1","--warmups","0",
                    "--cases","single","l_free_cavity","aligned_gap"],check=True,timeout=100)
    rows=json.loads((output/"summary.json").read_text())
    assert len(rows)==9 and all(r["classification_passed"] for r in rows)
    for name in ("single","l_free_cavity","aligned_gap"):
        group={r["mode"]:r for r in rows if r["case"]==name}
        assert set(group)=={"volume","flat","local"}
        for r in group.values():
            assert r["core_ms"]<=r["build_validation_ms"] and r["strict_query_us"]>0
        assert group["local"]["sign_certificate_applied"]
        assert group["local"]["strict_sign_certified"]==(name!="aligned_gap")
    gap=next(r for r in rows if r["case"]=="aligned_gap" and r["mode"]=="local")
    assert gap["adjacency_faces"]==4 and gap["logical_merges"]==2
    assert gap["artificial_zero_faces"]==2 and gap["artificial_zero_edges"]==4
    cavity=next(r for r in rows if r["case"]=="l_free_cavity" and r["mode"]=="local")
    assert cavity["artificial_zero_faces"]==cavity["artificial_zero_edges"]==cavity["artificial_zero_vertices"]==0
    assert (output/"report.md").is_file()
print("Fixed-axis local benchmark preserves successful and incomplete candidates")
