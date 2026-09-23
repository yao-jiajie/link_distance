"""End-to-end containment and CLI contracts. Uses only the Python standard library."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile


tree, wrap, cloud = map(Path, sys.argv[1:])


def invoke(program, output, arguments, success=True):
    command = [str(program), str(cloud), str(output)]
    if program == tree:
        command += ["--pipeline", "alpha-tet", "--validation-samples", "100"]
    else:
        command += ["--volume-validation-samples", "100"]
    result = subprocess.run(command + arguments, text=True, capture_output=True, timeout=90)
    if (result.returncode == 0) != success:
        raise AssertionError(f"Unexpected exit: {command + arguments}\n{result.stdout}\n{result.stderr}")


with tempfile.TemporaryDirectory(prefix="cgal_voxel_test_") as folder:
    root = Path(folder)
    for index, (mode, parameters) in enumerate([
        ("voxel-factor", []),
        ("absolute", ["--alpha", "3", "--offset", "1"]),
        ("bbox-ratio", ["--alpha-ratio", "0.5", "--offset-ratio", "0.5"]),
    ]):
        output = root / f"tree_{index}"
        invoke(tree, output, ["--voxel-size", "1"] + parameters)
        metrics = json.loads((output / "benchmark.json").read_text())[0]
        assert metrics["parameter_mode"] == mode
        assert metrics["raw_points"] == 9
        assert metrics["occupied_voxels"] == 1
        assert metrics["voxel_points"] == 8
        for name in ("raw_points_outside_voxels", "raw_points_outside_wrap",
                     "raw_points_outside_tree", "false_negative", "false_positive"):
            assert metrics[name] == 0, (name, metrics[name])
        assert not json.loads((output / "validation_report.json").read_text())["failures"]
        rows = list(csv.DictReader((output / "benchmark.csv").open()))
        assert len(rows) == 1 and None not in rows[0]
        assert int(rows[0]["voxel_points"]) == metrics["voxel_points"]
        assert abs(metrics["core_runtime_ms"] + metrics["validation_ms"] - metrics["total_ms"]) < 1
        invoke(wrap, root / f"wrap_{index}.off", ["--voxel-size", "1"] + parameters)

    invalid = [
        ["--voxel-size", "0"],
        ["--voxel-size", "nan"],
        ["--voxel-size", "-0.01"],
        ["--alpha-voxel-factor", "3"],
        ["--preprocess", "voxel"],
        ["--voxel-size", "1", "--preprocess", "jet"],
        ["--voxel-size", "1", "--preprocess", "none"],
        ["--voxel-size", "1", "--alpha", "3"],
        ["--voxel-size", "1", "--alpha-ratio", "2"],
        ["--voxel-size", "1", "--alpha", "3", "--offset", "1", "--alpha-voxel-factor", "3"],
        ["--voxel-size", "1", "--alpha-ratio", "2", "--offset-ratio", "3", "--offset-voxel-factor", "1"],
        ["--alpha", "3", "--offset", "1", "--alpha-ratio", "2", "--offset-ratio", "3"],
    ]
    for index, parameters in enumerate(invalid):
        invoke(tree, root / f"invalid_{index}", parameters, success=False)
        invoke(wrap, root / f"invalid_{index}.off", parameters, success=False)

    # Eight corners alone do not certify the interior of their box. A tight
    # wrap around isolated corners must reject a raw point missed between them,
    # even with random sampling disabled and retry budget exhausted.
    missed = root / "missed"
    tight = ["--voxel-size", "1", "--alpha", "0.1", "--offset", "0.005",
             "--max-offset-retries", "0"]
    invoke(tree, missed, tight + ["--validation-samples", "0"], success=False)
    metrics = json.loads((missed / "benchmark.json").read_text())[0]
    assert metrics["raw_points_outside_wrap"] > 0
    assert metrics["raw_points_outside_tree"] > 0
    invoke(wrap, root / "missed.off", tight + ["--volume-validation-samples", "0"], success=False)

    recovered = root / "recovered"
    recovery = ["--voxel-size", "1", "--alpha", "0.5", "--offset", "0.05",
                "--offset-growth-factor", "10", "--max-offset-retries", "2"]
    invoke(tree, recovered, recovery)
    metrics = json.loads((recovered / "benchmark.json").read_text())[0]
    assert metrics["offset_retries"] == 1 and metrics["offset"] == 0.5
    assert metrics["raw_points_outside_wrap"] == 0 and metrics["raw_points_outside_tree"] == 0
    invoke(wrap, root / "recovered.off", recovery)
    invoke(tree, root / "no_recovery", recovery + ["--auto-offset", "false"], success=False)

    explicit = root / "explicit_factors"
    invoke(tree, explicit, ["--voxel-size", "1", "--alpha-voxel-factor", "2",
                           "--offset-voxel-factor", "0.5", "--pipeline", "alpha_tet"])
    metrics = json.loads((explicit / "benchmark.json").read_text())[0]
    assert metrics["alpha"] == 2 and metrics["requested_offset"] == 0.5
    invoke(wrap, root / "explicit_factors.off", ["--voxel-size", "1",
                                               "--alpha-voxel-factor", "2", "--offset-voxel-factor", "0.5"])

    # Legacy bbox sweep remains available, and cannot silently ignore factors.
    invoke(tree, root / "sweep", ["--voxel-size", "1", "--benchmark-alpha-ratios", "0.5,0.6",
                                 "--offset-ratio", "0.5"])
    invoke(tree, root / "invalid_sweep", ["--voxel-size", "1", "--benchmark-alpha-ratios", "0.5,0.6",
                                         "--offset-ratio", "0.5", "--alpha-voxel-factor", "3"], success=False)

    # Near-coplanar reduction must still contain every raw input sample.
    reduced = root / "reduced"
    invoke(tree, reduced, ["--voxel-size", "1", "--plane-reduction", "conservative",
                          "--plane-angle-deg", "10", "--plane-distance-tol", "0.2"])
    metrics = json.loads((reduced / "benchmark.json").read_text())[0]
    assert metrics["raw_points_outside_tree"] == 0 and metrics["false_negative"] == 0

print("Voxel pipeline CLI and containment tests passed")
