"""Exercise the reproducible frontend benchmark and its raw timing records."""
import json
from pathlib import Path
import statistics
import subprocess
import sys
import tempfile

repo = Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix="boundary_frontend_") as temporary:
    output = Path(temporary)/"results"
    subprocess.run([sys.executable, str(repo/"scripts/benchmark_boundary_frontend.py"), str(output),
        "--build-dir", sys.argv[1], "--cases", "closed_cavity", "camera", "--repeats", "2",
        "--baseline-binary", str(Path(sys.argv[1])/"build_halfspace_tree"), "--pipeline-repeats", "1"], check=True, timeout=120)
    rows = json.loads((output/"summary.json").read_text())
    assert [row["case"] for row in rows] == ["closed_cavity", "camera"]
    for row in rows:
        assert row["passed"] and row["exact_boundary"] and row["exact_registry"]
        assert 0 < row["supports"] <= row["patches"] <= row["faces"]
        for mode in ("legacy", "contiguous"):
            for metric in ("extraction_ms", "merge_ms", "registry_ms", "frontend_ms"):
                trials = row[mode][metric+"_trials"]
                assert len(trials) == row["repeats"] == 2 and all(t > 0 for t in trials)
                assert row[mode][metric] == statistics.median(trials)
        for metric in ("merge_ms", "registry_ms", "frontend_ms"):
            assert row[metric+"_reduction"] == 1-row["contiguous"][metric]/row["legacy"][metric]
    assert rows[0]["faces"] == 60 and rows[0]["patches"] == 12
    manifest = json.loads((output/"manifest.json").read_text())
    assert len(manifest["inputs"]) == 2 and manifest["warmups"] == 2
    assert "camera" in (output/"report.md").read_text()
    comparisons = json.loads((output/"end_to_end_summary.json").read_text())
    assert len(comparisons) == 2 and all(r["same_outputs"] and r["same_non_timing_metrics"] for r in comparisons)
print("Frontend benchmark: exact legacy comparison and paired timing records passed")
