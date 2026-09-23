"""Independent exported-tree coverage checks and approximation CLI regressions."""
import itertools
import json
import math
from pathlib import Path
import random
import subprocess
import sys
import tempfile

program = Path(sys.argv[1]).resolve()


def evaluate(tree, point):
    values = {}
    for node in tree["nodes"]:  # Export is postorder.
        if node["type"] == "leaf":
            p = tree["leaves"][node["leaf_id"]]
            values[node["id"]] = sum(a*b for a, b in zip(p["normal"], point))-p["offset"]
        else:
            op = min if node["type"] == "min" else max
            values[node["id"]] = op(values[i] for i in node["children"])
    return values[tree["root_id"]]


with tempfile.TemporaryDirectory(prefix="octree_hull_pipeline_") as temporary:
    root = Path(temporary)
    cloud = root / "l.xyz"
    points = [(.005,.005,.005), (.015,.005,.005), (.005,.015,.005)]
    cloud.write_text("\n".join(" ".join(map(str,p)) for p in points))
    common = [str(program), str(cloud), "", "--pipeline", "octree", "--voxel-size", ".01",
              "--validation-samples", "100"]

    def run(name, options, good=True):
        command = common.copy()
        output = root / name
        command[2] = str(output)
        result = subprocess.run(command+options, capture_output=True, text=True, timeout=30)
        assert (result.returncode == 0) == good, result.stdout+result.stderr
        return output

    flags = ["--octree-approx", "hull", "--max-volume-inflation", ".2", "--max-fill-distance", ".008"]
    hull = run("hull", flags)
    baseline = run("baseline", [])
    none = run("none", ["--octree-approx", "none"])
    for name in ("logic_tree_raw.json", "logic_tree_simplified.json", "convex_clusters.json"):
        assert (baseline/name).read_bytes() == (none/name).read_bytes()
    metrics = json.loads((hull/"benchmark.json").read_text())[0]
    report = json.loads((hull/"validation_report.json").read_text())
    assert report["passed"] and report["source_coverage_and_error_certificate"]
    assert report["raw_points_outside_tree"] == report["source_inside_tree_outside"] == 0
    assert report["scalar_mismatch"] == report["cluster_mismatch"] == 0
    assert not report["known_free_space_checked"]
    assert metrics["convex_cells"] == metrics["non_aabb_cells"] == 1
    assert metrics["planes_after_reduce"] == 7 and metrics["tree_nodes"] == 8
    assert metrics["baseline_tree_nodes"] == 15
    assert math.isclose(metrics["volume_inflation"], 1/6, rel_tol=1e-10)
    assert 0 < metrics["maximum_fill_distance_bound_m"] <= .008
    geometry = json.loads((hull/"convex_clusters.json").read_text())
    assert geometry["geometry_definition"] == "intersection_of_exported_halfspaces"
    c = geometry["clusters"][0]
    assert c["type"] == "convex_polyhedron" and sorted(c["source_voxel_ids"]) == [0,1,2]
    assert any(sum(x != 0 for x in p["normal"]) > 1 for p in c["planes"])
    reference = json.loads((hull/"octree.json").read_text())
    assert reference["representation"] == "fine_occupancy_reference"
    tree = json.loads((hull/"logic_tree_simplified.json").read_text())
    raw_tree = json.loads((hull/"logic_tree_raw.json").read_text())
    old_tree = json.loads((baseline/"logic_tree_simplified.json").read_text())
    for plane, leaf_id in zip(c["planes"], c["leaf_ids"]):
        leaf = tree["leaves"][leaf_id]
        assert plane["normal"] == leaf["normal"] and plane["offset"] == leaf["offset"]
    randomizer = random.Random(41)
    for key in reference["occupied_voxel_indices"]:
        for t in list(itertools.product((0,.25,.5,.75,1), repeat=3)) + [tuple(randomizer.random() for _ in range(3)) for _ in range(30)]:
            # Consistent rounded voxel endpoints, as in the C++ grid.
            point = [key[k]*.01 + t[k]*((key[k]+1)*.01-key[k]*.01) for k in range(3)]
            assert evaluate(tree, point) <= 0
            assert evaluate(tree, point) == evaluate(raw_tree, point)
    added = (.012,.012,.005)
    assert evaluate(old_tree, added) > 0 and evaluate(tree, added) <= 0
    outside = (.019,.019,.005)
    assert evaluate(tree, outside) > 0  # Not merely the enclosing AABB.
    for ratio, distance in (("0",".008"), (".2","0"), (".1",".008"), (".2",".001")):
        path = run("restricted_"+ratio+"_"+distance,
                   ["--octree-approx", "hull", "--max-volume-inflation", ratio, "--max-fill-distance", distance])
        m = json.loads((path/"benchmark.json").read_text())[0]
        assert m["convex_cells"] == 2 and m["volume_inflation"] == 0
    for option in ("--convex-max-voxels", "--convex-test-max-boxes"):
        path = run("capped_"+option, flags+[option,"1"])
        m = json.loads((path/"benchmark.json").read_text())[0]
        assert m["convex_cells"] == 2 and m["rejected_work_limit"] > 0
    run("missing_both", ["--octree-approx", "hull"], False)
    run("missing_distance", flags[:-2], False)
    for bad in ("-1", "nan", "inf", "1x"):
        run("bad_"+bad, flags[:-1]+[bad], False)
    for extras in (["--octree-pruning","none"], ["--octree-boundary","exact"],
                   ["--alpha",".1"], ["--convex-max-voxels","0"], ["--max-fill-distance",".02"]):
        run("bad_flags", flags+extras, False)
    # No legacy/Alpha Wrap metadata gets silently reused in the new backend.
    assert tree["pipeline"] == metrics["pipeline"] == "octree-hull"
    assert not (hull/"boundary.off").exists() and (hull/"convex_cells.off").exists()
print("Octree hull CLI, exact baseline compatibility, and exported geometry tests passed")
