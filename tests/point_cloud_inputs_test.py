"""Run all diverse XYZ inputs through the current hard/LSE pipeline and report."""
import json
from pathlib import Path
import statistics
import subprocess
import sys
import tempfile

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO/"scripts"))
from generate_point_cloud_cases import write_cases


def load(path):
    return json.loads(path.read_text())


with tempfile.TemporaryDirectory(prefix="point_cloud_inputs_") as temporary:
    root = Path(temporary); inputs = root/"inputs"; out = root/"results"
    manifest = write_cases(inputs)
    subprocess.run([sys.executable,str(REPO/"scripts/benchmark_point_cloud_inputs.py"),str(out),
        "--input-dir",str(inputs),"--build-dir",sys.argv[1],"--repeats","1","--warmups","0",
        "--validation-samples","100"],check=True,timeout=160)
    rows = load(out/"summary.json")
    assert [r["case"] for r in rows]==[c["name"] for c in manifest["cases"]] and len(rows)==20
    for r in rows:
        assert r["passed"] and r["deterministic_outputs"] and r["independent_sample_checks"]==24
        assert r["gradient_checks"]>0 and r["zero_rays_tested"]>0 and r["strata_count"]>0
        assert r["function_ms"]==statistics.median(r["function_trials_ms"])>0
        assert r["function_min_ms"]<=r["function_ms"]<=r["function_p95_ms"]
        record = load(out/r["case"]/"trials.json")[0]
        assert r["function_ms"]==record["core_runtime_ms"]
        assert record["core_runtime_ms"]>=record["geometry_core_ms"] and record["core_runtime_ms"]<record["total_ms"]
    keys = {r["case"]: set(map(tuple,load(out/r["case"]/"occupancy_reference.json")["occupied_keys"])) for r in rows}
    assert keys["ellipsoid_sparse"] < keys["ellipsoid_surface"] < keys["ellipsoid_dense"]
    assert keys["ellipsoid_duplicates"]==keys["ellipsoid_surface"]
    assert len(keys["ellipsoid_outliers"]-keys["ellipsoid_surface"])==6
    assert next(r for r in rows if r["case"]=="ellipsoid_duplicates")["same_geometry_peer"]=="ellipsoid_surface"
    assert sum(r["gap_probe_checks"] for r in rows)==8
    for case in ("torus_surface","hollow_box_surfaces","narrow_corridor_surface","ellipsoid_partial"):
        assert next(r for r in rows if r["case"]==case)["gap_probe_checks"]>0
    report = (out/"report.md").read_text()
    assert "Function time excludes validation and all file I/O" in report
    assert all(r["case"] in report for r in rows)
print("Diverse point clouds: 20 pipelines, all-strata certificates, gap probes, exported functions/gradients and function-only timing passed")
