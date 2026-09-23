"""Phase 3 exports: independent face/patch/mesh reconstruction and unchanged trees."""
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

program = Path(sys.argv[1]).resolve()


def cross(a, b):
    return tuple(a[(i+1) % 3]*b[(i+2) % 3] - a[(i+2) % 3]*b[(i+1) % 3] for i in range(3))


def difference(a, b):
    return tuple(x-y for x, y in zip(a, b))


def verify(output):
    metrics = json.loads((output / "benchmark.json").read_text())[0]
    report = json.loads((output / "validation_report.json").read_text())
    assert report["passed"] and report["boundary_certificate"]
    assert report["raw_points_outside_tree"] == report["false_negative"] == report["false_positive"] == 0
    assert metrics["phase"] == 3 and metrics["volume_inflation"] == 0 and metrics["boundary_mode"] == "exact"
    tree = json.loads((output / "octree.json").read_text())
    voxels = [tuple(p) for p in tree["occupied_voxel_indices"]]
    occupied = set(voxels)
    expected = {}
    for voxel_id, voxel in enumerate(voxels):
        for axis, sign in itertools.product(range(3), (-1, 1)):
            neighbor = list(voxel)
            neighbor[axis] += sign
            if tuple(neighbor) not in occupied:
                expected[(axis, sign, voxel[axis]+(sign == 1), voxel[(axis+1) % 3], voxel[(axis+2) % 3])] = voxel_id
    boundary = json.loads((output / "boundary_patches.json").read_text())
    assert boundary["length_unit"] == "m" and not boundary["is_cluster_support_set"]
    faces, patches = boundary["faces"], boundary["patches"]
    assert len(faces) == len(expected) == metrics["boundary_faces_exposed"]
    assert len(patches) == metrics["boundary_patches"]
    face_keys = [(f["axis"], f["sign"], f["plane_index"], f["u"], f["v"]) for f in faces]
    assert len(set(face_keys)) == len(faces) and set(face_keys) == expected.keys()
    clusters = json.loads((output / "convex_clusters.json").read_text())["clusters"]
    owner = {v: c["id"] for c in clusters for v in c["source_voxel_ids"]}
    for f, k in zip(faces, face_keys):
        assert f["source_voxel_id"] == expected[k] and f["convex_cell_id"] == owner[expected[k]]
    covered = []
    for patch in patches:
        u0, u1, v0, v1 = patch["uv_bounds"]
        expected_patch = {(patch["axis"], patch["sign"], patch["plane_index"], u, v)
                          for u in range(u0, u1) for v in range(v0, v1)}
        ids = patch["source_face_ids"]
        assert {face_keys[i] for i in ids} == expected_patch and len(ids) == len(expected_patch)
        assert set(patch["convex_cell_ids"]) == {faces[i]["convex_cell_id"] for i in ids}
        covered.extend(ids)
        corners = patch["corner_indices"]
        n = cross(difference(corners[1], corners[0]), difference(corners[3], corners[0]))
        assert n == tuple(x*len(ids) for x in patch["normal"])
        for p, index in zip(patch["corners_m"], corners):
            assert p == [x*boundary["voxel_size"] for x in index]
            assert sum(x*y for x, y in zip(p, patch["normal"])) == patch["offset"]
    assert sorted(covered) == list(range(len(faces)))
    assert metrics["boundary_faces_before"] == 6*len(voxels)
    assert metrics["internal_faces_removed"] == 6*len(voxels)-len(faces)
    # Independently check OFF face count, winding, area and signed enclosed volume.
    lines = (output / "boundary.off").read_text().splitlines()
    assert lines[0] == "OFF"
    vertex_count, triangle_count, _ = map(int, lines[1].split())
    vertices = [tuple(map(float, row.split())) for row in lines[2:2+vertex_count]]
    triangles = [list(map(int, row.split())) for row in lines[2+vertex_count:]]
    assert len(triangles) == triangle_count == 2*len(faces) == metrics["boundary_triangles"]
    volume6 = 0
    for i, triangle in enumerate(triangles):
        assert triangle[0] == 3
        a, b, c = [vertices[j] for j in triangle[1:]]
        n = cross(difference(b, a), difference(c, a))
        f = faces[i//2]
        assert n[f["axis"]]*f["sign"] > 0
        assert all(n[j] == 0 for j in range(3) if j != f["axis"])
        volume6 += sum(x*y for x, y in zip(a, cross(b, c)))
    assert math.isclose(volume6/6, metrics["occupied_volume_m3"], rel_tol=1e-10, abs_tol=1e-12)


with tempfile.TemporaryDirectory(prefix="octree_boundary_pipeline_") as folder:
    root = Path(folder)
    shapes = {
        "single": [(0.5, 0.5, 0.5)],
        "cube": list(itertools.product((-0.5, 0.5, 1.5), repeat=3)),
        "shell": [p for p in itertools.product((-0.5, 0.5, 1.5), repeat=3) if p != (0.5, 0.5, 0.5)],
        "frame": [(x+.5, y+.5, .5) for x, y in itertools.product(range(3), repeat=2) if (x, y) != (1, 1)],
        "l": [(.5, .5, .5), (1.5, .5, .5), (.5, 1.5, .5)],
        "u": [(x+.5, y+.5, .5) for x, y in itertools.product(range(3), repeat=2) if x != 1 or y == 0],
        "edge_vertex": [(.5, .5, .5), (1.5, 1.5, .5), (2.5, 2.5, 1.5)],
        "negative_decimal": [(-.025, .015, .015), (.005, .015, .015), (.005, .015, .015)],
    }
    for name, points in shapes.items():
        cloud = root / f"{name}.xyz"
        cloud.write_text("\n".join(" ".join(map(str, p)) for p in points))
        for pruning in ("none", "exact", "rect"):
            paths = []
            for mode in ("none", "exact"):
                output = root / f"{name}_{pruning}_{mode}"
                paths.append(output)
                result = subprocess.run([str(program), str(cloud), str(output), "--pipeline", "octree",
                                         "--voxel-size", ".03" if name == "negative_decimal" else "1",
                                         "--octree-pruning", pruning, "--octree-boundary", mode,
                                         "--validation-samples", "0"], capture_output=True, text=True, timeout=30)
                assert result.returncode == 0, result.stdout + result.stderr
                if mode == "exact":
                    verify(output)
                else:
                    report = json.loads((output / "validation_report.json").read_text())
                    assert not report["boundary_enabled"] and report["boundary_certificate"] is None
                    assert not (output / "boundary.off").exists()
            # Stronger than sample consistency: every plane, node, ID, expression,
            # and convex partition is byte-for-byte unchanged with Phase 3 enabled.
            for file in ("logic_tree_raw.json", "logic_tree_simplified.json", "convex_clusters.json",
                         "logic_expression_raw.txt", "logic_expression_simplified.txt"):
                assert (paths[0] / file).read_bytes() == (paths[1] / file).read_bytes(), file
    for bad in (["--octree-boundary", "approx"], ["--octree-boundary"],
                ["--octree-boundary", "none", "--octree-boundary", "exact"]):
        result = subprocess.run([str(program), str(cloud), str(root / "bad"), "--pipeline", "octree",
                                 "--voxel-size", "1", *bad], capture_output=True, timeout=30)
        assert result.returncode != 0
    # Reusing a directory with boundary disabled must not silently present old
    # artifacts as current or destructively delete user-owned output files.
    previous_mesh = (paths[1] / "boundary.off").read_bytes()
    result = subprocess.run([str(program), str(cloud), str(paths[1]), "--pipeline", "octree",
                             "--voxel-size", ".03", "--validation-samples", "0"],
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0 and "previous run" in result.stderr
    assert (paths[1] / "boundary.off").read_bytes() == previous_mesh
    assert json.loads((paths[1] / "benchmark.json").read_text())[0]["boundary_mode"] == "none"

print("Phase 3 CLI, cavity preservation, oriented exports and unchanged tree tests passed")
