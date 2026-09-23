"""Phase 4 generator, query adapter, matrix export and failure-contract tests."""
import copy
import csv
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo / "scripts"))
from benchmark_scene_clouds import SCENES, inside_scene, interior_points, noisy_points, solid_volume, surface_points
from benchmark_octree_phase4 import interval, mesh_volume, voxel_contains

build = Path(sys.argv[1]).resolve()
for name, scene in SCENES.items():
    points = surface_points(scene, 200, 123)
    assert surface_points(scene, 20, 123) == points[:20]
    assert points == noisy_points(points, 0., 456)
    noisy = noisy_points(points, .003, 456)
    assert all(abs(x-y) <= .003+1e-16 for p, q in zip(points, noisy) for x, y in zip(p, q))
    assert all(inside_scene(p, scene) for p in interior_points(scene, 200, 789))
    assert solid_volume(scene) > 0
    if "gap" in scene:
        lo, hi = scene["gap"]
        assert all(not inside_scene(p, scene) for p in itertools.product(*zip(lo, hi)))
assert math.isclose(solid_volume(SCENES["cube"]), .027)
assert math.isclose(solid_volume(SCENES["thin_wall"]), .00054)
assert math.isclose(solid_volume(SCENES["l"]), .006885)
assert math.isclose(solid_volume(SCENES["u"]), .00972)
assert interval(0, 100)[0] == 0 and interval(0, 100)[1] > 0
assert interval(100, 100)[1] == 1 and interval(100, 100)[0] < 1
assert voxel_contains((.03, .03, .03), {(0, 0, 0)}, .03)
assert voxel_contains((-.03, -.03, -.03), {(-1, -1, -1)}, .03)
assert not voxel_contains((.031, .03, .03), {(0, 0, 0)}, .03)


with tempfile.TemporaryDirectory(prefix="phase4_benchmark_test_") as temporary:
    root = Path(temporary)
    tree = {"length_unit": "m", "halfspace_offset_unit": "m",
            "leaves": [{"id": 0, "normal": [1, 0, 0], "point": [1, 0, 0], "offset": 1},
                       {"id": 1, "normal": [-1, 0, 0], "point": [0, 0, 0], "offset": 0}],
            "nodes": [{"id": 99, "type": "max", "leaf_id": -1, "children": [4, 3]},
                      {"id": 3, "type": "leaf", "leaf_id": 1, "children": []},
                      {"id": 4, "type": "leaf", "leaf_id": 0, "children": []}], "root_id": 99}
    tree_path, query_path, output_path = root / "tree.json", root / "points.xyz", root / "query.json"
    tree_path.write_text(json.dumps(tree))
    query_path.write_text("# meter queries\n-.5 0 0\n0 0 0\n.25 0 0\n1 0 0\n1.5 0 0\n")
    command = [str(build / "halfspace_tree_benchmark"), str(tree_path), str(query_path), str(output_path),
               "--timed-count", "5", "--rounds", "2", "--repeats", "2"]
    subprocess.run(command, capture_output=True, check=True, timeout=30)
    result = json.loads(output_path.read_text())
    assert result["values"] == [.5, 0, -.25, 0, .5] and result["median_us_per_query"] >= 0
    assert len(result["timings_us_per_query"]) == 2 and not result["reference_checked"]
    reference = root / "reference.json"
    changed = copy.deepcopy(tree)
    changed["leaves"][0]["offset"] = .5
    reference.write_text(json.dumps(changed))
    run = subprocess.run(command+["--reference", str(reference)], capture_output=True, timeout=30)
    assert run.returncode == 2
    assert json.loads(output_path.read_text())["reference_decision_mismatches"] == 1
    mutations = []
    bad = copy.deepcopy(tree); bad["length_unit"] = "mm"; mutations.append(bad)
    bad = copy.deepcopy(tree); bad["nodes"][0]["children"] = [99]; mutations.append(bad)
    bad = copy.deepcopy(tree); bad["nodes"][0]["children"] = [101]; mutations.append(bad)
    bad = copy.deepcopy(tree); bad["nodes"][0]["children"] = []; mutations.append(bad)
    bad = copy.deepcopy(tree); bad["nodes"][0]["type"] = "sum"; mutations.append(bad)
    bad = copy.deepcopy(tree); bad["nodes"][1]["leaf_id"] = 2; mutations.append(bad)
    bad = copy.deepcopy(tree); bad["leaves"][1]["id"] = 0; mutations.append(bad)
    bad = copy.deepcopy(tree); bad["nodes"].append(copy.deepcopy(bad["nodes"][0])); mutations.append(bad)
    for bad in mutations:
        tree_path.write_text(json.dumps(bad))
        assert subprocess.run(command, capture_output=True, timeout=30).returncode != 0
    tree_path.write_text(json.dumps(tree))
    for bad in ("", "nan 0 0", "1 2", "1 2 3 4", "inf 0 0"):
        query_path.write_text(bad)
        assert subprocess.run(command, capture_output=True, timeout=30).returncode != 0
    query_path.write_text("0 0 0\n")
    for options in (["--rounds", "0"], ["--unknown", "1"], ["--reference"], ["--rounds", "-1"]):
        assert subprocess.run(command+options, capture_output=True, timeout=30).returncode != 0

    measured = root / "external.xyz"
    measured.write_text("-.01 -.01 .51\n.02 .01 .55\n-.02 .02 .57\n.01 -.02 .52\n")
    suite = [sys.executable, str(repo / "scripts/benchmark_octree_phase4.py"), str(root / "matrix"),
             "--build-dir", str(build), "--scenes", "u", "thin_wall", "--densities", "32",
             "--noise-levels", "0", ".003", "--voxel-sizes", ".04", "--skip-camera", "--include-baselines",
             "--measured-input", str(measured), "--repeats", "2", "--warmups", "1",
             "--validation-samples", "20", "--query-samples", "32", "--truth-samples", "16",
             "--query-rounds", "1", "--query-repeats", "2"]
    run = subprocess.run(suite, capture_output=True, text=True, timeout=110)
    assert run.returncode == 0, run.stdout + run.stderr
    summary = json.loads((root / "matrix/summary.json").read_text())
    assert summary["status"] == "passed" and summary["matrix_complete"] and summary["completed_cases"] == 5
    assert summary["measured_clouds"] == 1 and not summary["failures"]
    for case in summary["cases"]:
        assert len(case["backends"]) == 5
        for backend, record in case["backends"].items():
            assert record["status"] == "passed" and len(record["trials"]) == 2
            assert record["quality"]["raw"]["outside"] == 0
            assert record["query_timing"]["reference_scalar_mismatches"] == 0
            if backend != "alpha-tet":
                assert record["quality"]["occupied_union_certificate"] and record["cross_pruning_decision_mismatch"] == 0
                assert record["quality"]["common_uniform"]["added_samples_vs_voxels"] == 0
            else:
                assert not record["quality"]["occupied_union_certificate"]
                assert mesh_volume(root / "matrix" / case["case"] / backend / "wrap.off") > 0
        if case["scene"] == "measured":
            assert "ideal_interior" not in case["backends"]["octree-rect"]["quality"]
    rows = list(csv.DictReader((root / "matrix/summary.csv").open()))
    assert len(rows) == 25 and all(row["status"] == "passed" for row in rows)
    # Existing results must be preserved, never silently overwritten on rerun.
    old_summary = (root / "matrix/summary.json").read_bytes()
    assert subprocess.run(suite, capture_output=True, timeout=30).returncode != 0
    assert (root / "matrix/summary.json").read_bytes() == old_summary
    for i, bad in enumerate((["--voxel-sizes", "0"], ["--noise-levels", "nan"], ["--densities", "1"],
                             ["--scenes", "u", "u"], ["--query-rounds", "0"], ["--repeats", "0"])):
        command = [sys.executable, str(repo / "scripts/benchmark_octree_phase4.py"), str(root / f"bad_{i}"),
                   "--build-dir", str(build), *bad]
        assert subprocess.run(command, capture_output=True, timeout=30).returncode != 0
        assert not (root / f"bad_{i}").exists()

print("Phase 4 sampling, evaluator, multi-mode exports, strict coverage and CLI tests passed")
