#!/usr/bin/env python3
"""Reproducibility, complexity, and non-collision regression for the six-shape cloud."""

import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from generate_panda_complex_scene import generate


def main():
    if len(sys.argv) != 4:
        raise SystemExit("expected link_distance_node, Panda URDF, and complex XYZ")
    node, urdf, cloud = map(Path, sys.argv[1:])
    points, metadata = generate()
    assert len(points) == 9216
    assert len(metadata["components"]) == 6
    assert metadata["recommended_voxel_size_m"] == 0.02
    assert cloud.read_text() == "".join(
        f"{x:.17g} {y:.17g} {z:.17g}\n" for x, y, z in points)
    assert json.loads(cloud.with_suffix(".json").read_text()) == metadata

    with tempfile.TemporaryDirectory(prefix="panda_complex_scene_") as folder:
        output = Path(folder) / "distances.json"
        completed = subprocess.run(
            [
                str(node), str(urdf), str(cloud), str(output),
                "--voxel-size", "0.02",
                "--q", "0,-0.4,0,-2.0,0,1.6,0.8",
                "--env-lse-error", "0.001",
                "--link-lse-error", "0.001",
                "--margin", "0.005",
                "--field-from-base", "0,0,0,0,0,0",
                "--query-workers", "6",
                "--boundary-max-cells", "262144",
                "--boundary-max-expanded-nodes", "32000000",
                "--boundary-max-strata", "2000000",
                "--boundary-max-sign-ops", "40000000000",
            ],
            text=True, capture_output=True, timeout=60,
        )
        assert completed.returncode == 0, completed.stderr
        result = json.loads(output.read_text())

    assert result["point_count"] == len(points)
    assert result["voxel_occupancy"]["occupied_voxels"] == 3989
    assert result["voxel_occupancy"]["merged_boundary_patches"] == 6869
    assert len(result["links"]) == 7
    assert len(result["spheres"]) == 59
    assert all(math.isfinite(link["distance_proxy_m"]) and
               link["distance_proxy_m"] > 0.1 and
               all(math.isfinite(value) for value in link["gradient_dq"])
               for link in result["links"])
    assert all(sphere["clearance_proxy_m"] > 0.0
               for sphere in result["spheres"])
    raw_gap = min(
        math.dist(point, sphere["center_field_m"]) - sphere["radius_m"]
        for sphere in result["spheres"] for point in points
    )
    assert raw_gap > 0.05, raw_gap
    print(f"Complex Panda cloud: {len(points)} points, six shapes, "
          f"{result['voxel_occupancy']['occupied_voxels']} voxels, "
          f"minimum raw sphere-point gap {raw_gap:.3f} m")


if __name__ == "__main__":
    main()
