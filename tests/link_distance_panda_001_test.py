#!/usr/bin/env python3
"""Regression for the certified Panda/camera link-distance build at h=0.01 m."""

import json
import math
import pathlib
import subprocess
import sys
import tempfile


def main():
    if len(sys.argv) != 4:
        raise SystemExit("expected link_distance_node, Panda URDF, and camera XYZ")
    node, urdf, cloud = map(pathlib.Path, sys.argv[1:])
    with tempfile.TemporaryDirectory() as folder:
        output = pathlib.Path(folder) / "panda_001.json"
        completed = subprocess.run(
            [
                str(node), str(urdf), str(cloud), str(output),
                "--voxel-size", "0.01",
                "--q", "0,-0.4,0,-2.0,0,1.6,0.8",
                "--env-lse-error", "0.001",
                "--link-lse-error", "0.001",
                "--margin", "0.005",
                "--field-from-base", "0,0.3,0,0,0,0",
            ],
            text=True,
            capture_output=True,
            timeout=30,
        )
        assert completed.returncode == 0, completed.stderr
        result = json.loads(output.read_text())
        limits = result["environment_direct_limits"]
        assert result["voxel_size_m"] == 0.01
        assert result["voxel_occupancy"]["occupied_voxels"] == 1562
        assert len(result["links"]) == 7
        assert len(result["spheres"]) == 59
        assert limits["max_canonical_cells"] >= 62400
        assert limits["max_expanded_nodes"] >= 7_300_980
        assert limits["max_strata"] >= 480083
        assert limits["max_sign_operations"] > 5_000_000_000
        assert math.isfinite(result["environment_build_ms"])
        assert all(link["distance_proxy_m"] > 0 for link in result["links"])


if __name__ == "__main__":
    main()
