#!/usr/bin/env python3
import json
import math
import pathlib
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET


def main():
    if len(sys.argv) != 7:
        raise SystemExit(
            "expected node, URDF, XYZ, SVG/3D visualizer, and frame validator paths"
        )
    node, urdf, cloud, visualizer, visualizer_3d, frame_validator = map(
        pathlib.Path, sys.argv[1:]
    )
    with tempfile.TemporaryDirectory() as folder:
        output = pathlib.Path(folder) / "link_distances.json"
        completed = subprocess.run(
            [
                str(node),
                str(urdf),
                str(cloud),
                str(output),
                "--voxel-size",
                "0.1",
                "--q",
                "0.3,0.1",
                "--env-lse-error",
                "0.001",
                "--link-lse-error",
                "0.001",
                "--margin",
                "0.005",
            ],
            text=True,
            capture_output=True,
        )
        assert completed.returncode == 0, completed.stderr
        result = json.loads(output.read_text())
        assert "query_phase_ms" not in result
        assert result["schema"] == 1
        assert result["value_kind"] == "smooth_clearance_proxy_not_euclidean_sdf"
        assert result["joint_order"] == ["joint1", "joint2"]
        assert result["configuration"] == [0.3, 0.1]
        assert result["point_count"] > 0
        occupancy = result["voxel_occupancy"]
        assert occupancy["occupied_voxels"] > 0
        assert occupancy["faces_before_common_face_removal"] == 6 * occupancy["occupied_voxels"]
        assert occupancy["faces_before_common_face_removal"] == (
            2 * occupancy["internal_common_faces_removed"]
            + occupancy["exposed_unit_faces"]
        )
        assert occupancy["internal_directed_faces_removed"] == 2 * occupancy["internal_common_faces_removed"]
        assert occupancy["merged_boundary_patches"] == len(occupancy["boundary_patches"])
        assert occupancy["coplanar_exposed_faces_merged"] == (
            occupancy["exposed_unit_faces"] - occupancy["merged_boundary_patches"]
        )
        assert 0 < occupancy["merged_boundary_patches"] <= occupancy["exposed_unit_faces"]
        for patch in occupancy["boundary_patches"]:
            assert patch["axis"] in (0, 1, 2)
            assert patch["sign"] in (-1, 1)
            assert len(patch["corners_field_m"]) == 4
            assert all(
                len(corner) == 3 and all(math.isfinite(value) for value in corner)
                for corner in patch["corners_field_m"]
            )
        assert len(result["links"]) == 2
        assert len(result["spheres"]) == 4
        assert {sphere["link"] for sphere in result["spheres"]} == {
            "base",
            "link1",
            "link2",
        }
        for sphere in result["spheres"]:
            assert len(sphere["center_field_m"]) == 3
            assert len(sphere["center_urdf_link_m"]) == 3
            assert len(sphere["spatial_gradient"]) == 3
            assert sphere["radius_m"] > 0
            assert math.isfinite(sphere["clearance_proxy_m"])
            if sphere["link"] == "base":
                assert sphere["joint_group"] is None
                assert sphere["group_link"] is None
                assert sphere["joint_configuration_index"] is None
                assert sphere["center_group_link_m"] is None
            else:
                expected = int(sphere["link"][-1]) - 1
                assert sphere["joint_group"] == f"joint{expected + 1}"
                assert sphere["group_link"] == sphere["link"]
                assert sphere["joint_configuration_index"] == expected
                assert sphere["center_group_link_m"] == sphere["center_urdf_link_m"]
        for link in result["links"]:
            assert link["joint"] in result["joint_order"]
            assert result["joint_order"][link["joint_configuration_index"]] == link["joint"]
            assert link["source_links"] == [link["name"]]
            assert len(link["gradient_dq"]) == 2
            assert all(math.isfinite(value) for value in link["gradient_dq"])
            assert math.isfinite(link["distance_proxy_m"])
            assert link["distance_proxy_m"] <= link["hard_min_proxy_m"] + 1e-14
            assert (
                link["hard_min_proxy_m"] - link["distance_proxy_m"]
                <= link["link_smoothing_bound_m"] + 1e-14
            )

        profiled_output = pathlib.Path(folder) / "profiled.json"
        profiled = subprocess.run(
            [
                str(node), str(urdf), str(cloud), str(profiled_output),
                "--voxel-size", "0.1", "--q", "0.3,0.1",
                "--env-lse-error", "0.001", "--link-lse-error", "0.001",
                "--margin", "0.005", "--profile-query", "true",
            ], text=True, capture_output=True,
        )
        assert profiled.returncode == 0, profiled.stderr
        profiled_result = json.loads(profiled_output.read_text())
        phases = profiled_result.pop("query_phase_ms")
        assert set(phases) == {"kinematics", "field", "link"}
        assert all(math.isfinite(value) and value >= 0 for value in phases.values())
        for key in ("environment_build_ms", "query_ms"):
            profiled_result.pop(key)
            result.pop(key)
        assert profiled_result == result

        validated = subprocess.run(
            [
                sys.executable,
                str(frame_validator),
                str(output),
                "--urdf",
                str(urdf),
            ],
            text=True,
            capture_output=True,
        )
        assert validated.returncode == 0, validated.stderr
        audit = json.loads(validated.stdout)
        assert audit["passed"] and audit["spheres"] == 4

        svg = pathlib.Path(folder) / "link_distances.svg"
        visualized = subprocess.run(
            [
                sys.executable,
                str(visualizer),
                str(output),
                str(svg),
                "--cloud",
                str(cloud),
            ],
            text=True,
            capture_output=True,
        )
        assert visualized.returncode == 0, visualized.stderr
        root = ET.parse(svg).getroot()
        namespace = {"svg": "http://www.w3.org/2000/svg"}
        circles = root.findall(".//svg:circle", namespace)
        assert len(circles) >= result["point_count"] + len(result["spheres"])

        html = pathlib.Path(folder) / "link_distances_3d.html"
        visualized_3d = subprocess.run(
            [
                sys.executable,
                str(visualizer_3d),
                str(output),
                str(html),
                "--cloud",
                str(cloud),
            ],
            text=True,
            capture_output=True,
        )
        assert visualized_3d.returncode == 0, visualized_3d.stderr
        page = html.read_text()
        assert "id=\"scene-data\"" in page
        assert "Drag: orbit" in page
        assert "Voxel surface" in page
        assert "__SCENE_DATA__" not in page
        assert all(link["name"] in page for link in result["links"])


if __name__ == "__main__":
    main()
