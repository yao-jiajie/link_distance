"""Phase 1/2 geometry and CLI contracts, independently evaluating exported JSON."""
import csv
import itertools
import json
import math
from pathlib import Path
import random
import subprocess
import sys
import tempfile

program = Path(sys.argv[1])


def evaluate(tree, point):
    leaves = {p["id"]: p for p in tree["leaves"]}
    values = {}
    for node in tree["nodes"]:
        if node["type"] == "leaf":
            plane = leaves[node["leaf_id"]]
            value = sum(n * x for n, x in zip(plane["normal"], point)) - plane["offset"]
        else:
            op = min if node["type"] == "min" else max
            value = op(values[i] for i in node["children"])
        values[node["id"]] = value
    return values[tree["root_id"]]


with tempfile.TemporaryDirectory(prefix="octree_pipeline_test_") as folder:
    root = Path(folder)

    def run(name, points, arguments=(), success=True, literal=None, pruning="none"):
        cloud = root / f"{name}.xyz"
        cloud.write_text(literal if literal is not None else "\n".join(" ".join(map(str, p)) for p in points))
        output = root / name
        mode = ["--octree-pruning", pruning] if pruning is not None else []
        result = subprocess.run([str(program), str(cloud), str(output), "--pipeline", "octree",
                                 "--validation-samples", "100", *mode, *arguments],
                                capture_output=True, text=True, timeout=30)
        assert (result.returncode == 0) == success, (name, result.stdout, result.stderr)
        if not success:
            return result
        metrics = json.loads((output / "benchmark.json").read_text())[0]
        validation = json.loads((output / "validation_report.json").read_text())
        assert validation["passed"] and validation["occupied_union_certificate"] and not validation["failures"]
        assert validation["pruning_certificate"] and validation["pruning_tree_comparison_samples"] >= 100
        assert validation["rectangular_partition_certificate"]
        assert not validation["pruning_scalar_equality_required"]
        for key in ["raw_points_outside_voxels", "raw_points_outside_tree", "false_negative", "false_positive",
                    "cluster_mismatch", "raw_simplified_mismatch", "scalar_mismatch", "pruning_tree_decision_mismatch"]:
            assert validation[key] == 0, (name, key)
        cells = json.loads((output / "convex_clusters.json").read_text())["clusters"]
        before = json.loads((output / "logic_tree_raw.json").read_text())
        after = json.loads((output / "logic_tree_simplified.json").read_text())
        count = len(cells)
        assert count == metrics["convex_cells"] == metrics["output_leaf_nodes"]
        fine_count = metrics["leaf_nodes"]
        assert fine_count == metrics["occupied_voxels"] == metrics["convex_cells_before_prune"]
        assert count + 7 * metrics["pruning_merges"] + metrics["rectangular_cells_removed"] == fine_count
        assert metrics["octree_leaf_nodes_after_prune"] == metrics["convex_cells_before_rectangles"]
        assert metrics["tree_nodes_before_prune"] == 7 * fine_count + (fine_count > 1)
        assert metrics["planes_before_prune"] == 6 * fine_count and metrics["planes_after_prune"] == 6 * count
        assert metrics["planes_before_reduce"] == metrics["planes_after_reduce"] == 6 * count
        assert metrics["tree_nodes"] == 7 * count + (count > 1)
        assert metrics["volume_inflation"] == 0 and metrics["pruning_time"] >= 0
        assert math.isclose(metrics["occupied_volume_m3"], fine_count * metrics["voxel_resolution"] ** 3)
        assert math.isclose(metrics["output_volume_m3"], metrics["occupied_volume_m3"])
        if pruning == "none":
            assert metrics["pruning_merges"] == 0 and count == fine_count and metrics["phase"] == 1
        else:
            assert metrics["phase"] == 2 and metrics["pruning_mode"] == (pruning or "rect")
        hierarchy = json.loads((output / "octree.json").read_text())
        nodes = {n["id"]: n for n in hierarchy["nodes"]}
        assert len(nodes) == metrics["octree_nodes_after_prune"] <= metrics["octree_nodes"]
        assert hierarchy["source_octree_nodes"] == metrics["octree_nodes"]
        assert sum(math.prod(c["extent_voxels"]) for c in cells) == fine_count
        assert sorted(v for c in cells for v in c["source_voxel_ids"]) == list(range(fine_count))
        assert sorted(p for n in nodes.values() for p in n["point_indices"]) == list(range(metrics["raw_points"]))
        for node in nodes.values():
            children = [nodes[i] for i in node["children"] if i >= 0]
            if node["convex_cell_id"] >= 0:
                assert not children and len(node["point_indices"]) == node["point_count"]
            else:
                assert children and not node["point_indices"]
                assert sum(c["occupied_count"] for c in children) == node["occupied_count"]
                assert sum(c["point_count"] for c in children) == node["point_count"]
        for cell in cells:
            assert cell["occupied_count"] == math.prod(cell["extent_voxels"])
            if cell["octree_node"] >= 0:
                assert cell["octree_node"] in nodes and cell["width_voxels"] ** 3 == cell["occupied_count"]
            else:
                assert metrics["pruning_mode"] == "rect" and cell["width_voxels"] is None
                assert cell["source_parent_node"] in nodes and 1 < cell["sibling_mask"] < 255
            assert len(cell["source_voxel_ids"]) == cell["occupied_count"]
        assert abs(metrics["core_runtime_ms"] + metrics["validation_ms"] - metrics["total_ms"]) < 0.1
        rows = list(csv.DictReader((output / "benchmark.csv").open()))
        assert len(rows) == 1 and None not in rows[0]
        rng = random.Random(123)
        samples = list(points) + [tuple(rng.uniform(-4, 4) for _ in range(3)) for _ in range(200)]
        for cell in cells:
            samples.extend(itertools.product(*zip(cell["lower"], cell["upper"])))
        for point in samples:
            occupied = any(all(lo <= x <= hi for lo, x, hi in zip(c["lower"], point, c["upper"])) for c in cells)
            size = metrics["voxel_resolution"]
            reference = any(all(k * size <= x <= (k + 1) * size for k, x in zip(key, point))
                            for key in hierarchy["occupied_voxel_indices"])
            assert reference == occupied  # independently verify no voxel added/lost
            assert (evaluate(before, point) <= 0) == occupied
            assert evaluate(before, point) == evaluate(after, point)
        return metrics, after, output

    one, _, _ = run("one", [(-0.005, 0, 0)] * 3, ["--voxel-size", "0.01", "--max-depth", "0"])
    assert one["raw_points"] == 3 and one["convex_cells"] == 1
    block = list(itertools.product((-0.5, 0.5), repeat=3))
    full, tree, _ = run("block", block, ["--voxel-size", "1", "--max-depth", "1"])
    assert full["octree_nodes"] == 9 and full["convex_cells"] == 8  # no Phase 2 pruning
    _, shuffled_tree, _ = run("permuted", block[::-1], ["--voxel-size", "1"])
    assert shuffled_tree == tree
    l_shape = [(0.5, 0.5, 0.5), (1.5, 0.5, 0.5), (0.5, 1.5, 0.5)]
    _, tree, _ = run("l_shape", l_shape, ["--voxel-size", "1"])
    assert evaluate(tree, (1.5, 1.5, 0.5)) > 0
    u_shape = [(x + 0.5, y + 0.5, 0.5) for x in range(3) for y in range(3) if x in (0, 2) or y == 0]
    _, tree, _ = run("u_shape", u_shape, ["--voxel-size", "1"])
    assert evaluate(tree, (1.5, 1.5, 0.5)) > 0
    assert evaluate(tree, (1.5, 2.5, 0.5)) > 0
    for name, points in [("plane", [(x / 10, y / 10, 0) for x in range(-2, 3) for y in range(-2, 3)]),
                         ("line", [(x / 10, 0, 0) for x in range(-3, 4)])]:
        run(name, points, ["--voxel-size", "0.03"])
    run("comments", [(0, 0, 0)], ["--voxel-size", "1"], literal="# meters\n\n0 0 0 # origin\n")

    full2, tree2, _ = run("block_pruned", block, ["--voxel-size", "1"], pruning="exact")
    assert full2["convex_cells"] == 1 and full2["tree_nodes"] == 7 and full2["octree_nodes_after_prune"] == 1
    assert full2["pruning_merges"] == 1
    assert full2["maximum_pruning_scalar_change"] == 1
    assert evaluate(tree2, (0, 0, 0)) == -1  # internal grid planes are removed
    default, default_tree, _ = run("default_pruning", block, ["--voxel-size", "1"], pruning=None)
    assert default_tree == tree2 and default["pruning_mode"] == "rect"
    _, reverse_tree, _ = run("pruned_permutation", block[::-1], ["--voxel-size", "1"], pruning="exact")
    assert reverse_tree == tree2
    cube4 = list(itertools.product((-1.5, -0.5, 0.5, 1.5), repeat=3))
    dense, dense_tree, _ = run("cube4", cube4, ["--voxel-size", "1"], pruning="exact")
    assert dense["convex_cells"] == 1 and dense["pruning_merges"] == 9 and dense["tree_nodes"] == 7
    partial, _, _ = run("cube4_missing", cube4[:-1], ["--voxel-size", "1"], pruning="exact")
    assert partial["convex_cells"] == 14 and partial["pruning_merges"] == 7
    shell = [p for p in cube4 if any(abs(x) == 1.5 for x in p)]
    _, hollow_tree, _ = run("hollow", shell, ["--voxel-size", "1"], pruning="exact")
    assert evaluate(hollow_tree, (0, 0, 0)) > 0
    for name, points in [("l", l_shape), ("u", u_shape)]:
        _, exact_tree, _ = run(f"{name}_exact", points, ["--voxel-size", "1"], pruning="exact")
        assert evaluate(exact_tree, (1.5, 1.5, 0.5)) > 0
    for name, points, expected in [
        ("pair", [(0.5,0.5,0.5),(1.5,0.5,0.5)], 1),
        ("slab", list(itertools.product((0.5,1.5), (0.5,1.5), (0.5,))), 1),
        ("l_rect", l_shape, 2),
        ("plane_rect", list(itertools.product((0.5,1.5,2.5,3.5), (0.5,1.5,2.5,3.5), (0.5,))), 4),
    ]:
        rect, rect_tree, _ = run(name, points, ["--voxel-size", "1"], pruning="rect")
        assert rect["convex_cells"] == expected and rect["rectangular_cells_removed"] > 0
        _, reverse_rect, _ = run(name+"_reverse", points[::-1], ["--voxel-size", "1"], pruning="rect")
        assert rect_tree == reverse_rect
    _, u_rect, _ = run("u_rect", u_shape, ["--voxel-size", "1"], pruning="rect")
    assert evaluate(u_rect, (1.5,1.5,0.5)) > 0 and evaluate(u_rect, (1.5,2.5,0.5)) > 0
    _, hollow_rect, _ = run("hollow_rect", shell, ["--voxel-size", "1"], pruning="rect")
    assert evaluate(hollow_rect, (0,0,0)) > 0
    for i, mode in enumerate(["approximate", "false", "0"]):
        run(f"bad_pruning_{i}", block, ["--voxel-size", "1"], False, pruning=mode)

    for i, args in enumerate([
        [], ["--voxel-size", "0"], ["--voxel-size", "nan"], ["--voxel-size", "inf"],
        ["--voxel-size", "-1"], ["--voxel-size", "1junk"], ["--voxel-size", "1", "--max-depth", "0"],
        ["--voxel-size", "1", "--max-depth", "-1"], ["--voxel-size", "1", "--max-depth", "1.5"],
        ["--voxel-size", "1", "--max-depth", "53"], ["--voxel-size", "1", "--alpha", "1"],
        ["--voxel-size", "1", "--preprocess", "jet"], ["--voxel-size", "1", "--occupancy-error", "0.1"],
        ["--voxel-size", "1", "--voxel-size", "2"], ["--voxel-size"],
    ]):
        run(f"bad_arg_{i}", block, args, False)
    for i, literal in enumerate(["", "# no data", "nan 0 0", "inf 0 0", "1 2", "1 2 3 4", "0 0 0\nbad"]):
        run(f"bad_input_{i}", [], ["--voxel-size", "1"], False, literal)
    # Zero random samples must still validate raw points, all corners, and centers.
    cloud = root / "zero.xyz"
    cloud.write_text("0.5 0.5 0.5\n")
    output = root / "zero"
    subprocess.run([str(program), str(cloud), str(output), "--pipeline", "octree", "--voxel-size", "1",
                    "--validation-samples", "0"], capture_output=True, check=True, timeout=30)
    assert json.loads((output / "validation_report.json").read_text())["tested_points"] == 10

    # Reproduce the read-only estimate on the shipped sparse camera input.
    camera = Path(__file__).resolve().parents[1] / "data" / "realistic_sparse_camera.xyz"
    camera_output = root / "camera_rect"
    subprocess.run([str(program), str(camera), str(camera_output), "--pipeline", "octree",
                    "--voxel-size", "0.03", "--octree-pruning", "rect", "--validation-samples", "100"],
                   capture_output=True, check=True, timeout=30)
    camera_metrics = json.loads((camera_output / "benchmark.json").read_text())[0]
    assert camera_metrics["occupied_voxels"] == 543 and camera_metrics["convex_cells_before_rectangles"] == 529
    assert camera_metrics["convex_cells"] == 270 and camera_metrics["tree_nodes"] == 1891
    assert camera_metrics["planes_after_reduce"] == 1620 and camera_metrics["volume_inflation"] == 0
    report = json.loads((camera_output / "validation_report.json").read_text())
    assert report["passed"] and report["raw_points_outside_tree"] == 0 and report["pruning_tree_decision_mismatch"] == 0

print("Octree Phase 1/2 CLI, pruning, cavity, tree, and containment tests passed")
