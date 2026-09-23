#!/usr/bin/env python3
"""Paired end-to-end core-time benchmark for identical voxel occupancy across frames."""

import argparse
import json
import math
from pathlib import Path
import statistics
import subprocess
import tempfile


def shifted_inside_same_voxels(source: Path, destination: Path, voxel_size: float):
    rows = []
    for raw in source.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        point = [float(value) for value in line.split()]
        if len(point) != 3 or not all(math.isfinite(value) for value in point):
            raise ValueError(f"Invalid XYZ row: {raw}")
        moved = [
            value + 0.25 * ((math.floor(value / voxel_size) + 0.5) * voxel_size - value)
            for value in point
        ]
        rows.append(" ".join(f"{value:.17g}" for value in moved))
    if not rows:
        raise ValueError("Cloud is empty")
    destination.write_text("\n".join(rows) + "\n")


def remove_occupied_voxels(source: Path, destination: Path, voxel_size: float,
                           count: int):
    points = []
    for raw in source.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if line:
            point = tuple(float(value) for value in line.split())
            if len(point) != 3 or not all(math.isfinite(value) for value in point):
                raise ValueError(f"Invalid XYZ row: {raw}")
            key = tuple(math.floor(value / voxel_size) for value in point)
            points.append((point, key))
    keys = sorted({key for _, key in points})
    if not 0 < count < len(keys):
        raise ValueError("Invalid changed-voxel count")
    begin = (len(keys) - count) // 2
    removed = set(keys[begin:begin + count])
    destination.write_text("".join(
        " ".join(f"{value:.17g}" for value in point) + "\n"
        for point, key in points if key not in removed
    ))


def run(node, urdf, clouds, output, voxel_size, reuse, same_occupancy):
    command = [
        str(node), str(urdf), "-", "-",
        "--stream", "true", "--reuse-environment", str(reuse).lower(),
        "--voxel-size", str(voxel_size), "--q", "0,-0.4,0,-2,0,1.6,0.8",
        "--env-lse-error", "0.001", "--link-lse-error", "0.001",
        "--margin", "0.005", "--boundary-max-cells", "262144",
        "--boundary-max-direct-nodes", "16384",
        "--boundary-max-expanded-nodes", "32000000",
        "--boundary-max-split-checks", "1048576",
        "--boundary-max-strata", "2000000",
        "--boundary-max-sign-ops", "40000000000",
    ]
    outputs = [output.with_name(f"{output.stem}_frame_{index}.json")
               for index in range(len(clouds))]
    manifest = "".join(f"{cloud}\t{frame_output}\n"
                       for cloud, frame_output in zip(clouds, outputs))
    completed = subprocess.run(command, input=manifest, text=True, capture_output=True)
    if completed.returncode:
        raise RuntimeError(completed.stderr)
    timings = [json.loads(line) for line in completed.stdout.splitlines()]
    if len(timings) != len(clouds):
        raise RuntimeError("Node did not report one timing row per frame")
    assert not timings[0]["cache_hit"]
    assert all(frame["cache_hit"] == (reuse and same_occupancy)
               for frame in timings[1:])
    results = [json.loads(path.read_text()) for path in outputs]
    for result in results:
        for field in ("environment_build_ms", "query_ms"):
            result.pop(field)
    return timings[1:], results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("node", type=Path)
    parser.add_argument("urdf", type=Path)
    parser.add_argument("cloud", type=Path)
    parser.add_argument("--voxel-size", type=float, default=0.02)
    parser.add_argument("--pairs", type=int, default=4)
    parser.add_argument("--frames", type=int, default=3)
    parser.add_argument("--changed-voxels", type=int, default=0,
                        help="Remove this many occupied voxels instead of moving points inside their voxels")
    args = parser.parse_args()
    if (args.pairs < 1 or args.frames < 2 or args.changed_voxels < 0 or
            not math.isfinite(args.voxel_size) or args.voxel_size <= 0):
        parser.error("pairs must be >=1, frames >=2, changed voxels >=0, and voxel size positive")

    with tempfile.TemporaryDirectory() as folder:
        folder = Path(folder)
        moved = folder / "changed_points.xyz"
        if args.changed_voxels:
            remove_occupied_voxels(args.cloud, moved, args.voxel_size,
                                   args.changed_voxels)
        else:
            shifted_inside_same_voxels(args.cloud, moved, args.voxel_size)
        clouds = [args.cloud if index % 2 == 0 else moved for index in range(args.frames)]
        pairs = []
        for pair in range(args.pairs):
            values = {}
            for reuse in ((True, False) if pair % 2 else (False, True)):
                output = folder / f"pair_{pair}_{int(reuse)}.json"
                times, result = run(args.node, args.urdf, clouds, output,
                                    args.voxel_size, reuse,
                                    args.changed_voxels == 0)
                values[reuse] = {
                    "core_ms": statistics.median(frame["core_ms"] for frame in times),
                    "core_no_validation_ms": (
                        statistics.median(frame["core_no_validation_ms"] for frame in times)
                        if reuse and not args.changed_voxels else None
                    ),
                    "boundary_ms": statistics.median(frame["boundary_ms"] for frame in times),
                    "direct_ms": statistics.median(frame["direct_ms"] for frame in times),
                    "environment_ms": statistics.median(frame["environment_ms"] for frame in times),
                    "query_ms": statistics.median(frame["query_ms"] for frame in times),
                    "result": result,
                }
            if values[True]["result"] != values[False]["result"]:
                raise RuntimeError("Reused and rebuilt output values differ")
            pairs.append({
                "reuse_core_ms": values[True]["core_ms"],
                "reuse_core_no_validation_ms": values[True]["core_no_validation_ms"],
                "reuse_boundary_ms": values[True]["boundary_ms"],
                "rebuild_boundary_ms": values[False]["boundary_ms"],
                "reuse_direct_ms": values[True]["direct_ms"],
                "rebuild_direct_ms": values[False]["direct_ms"],
                "rebuild_core_ms": values[False]["core_ms"],
                "reuse_environment_ms": values[True]["environment_ms"],
                "rebuild_environment_ms": values[False]["environment_ms"],
                "reuse_query_ms": values[True]["query_ms"],
                "rebuild_query_ms": values[False]["query_ms"],
            })
    print(json.dumps({
        "voxel_size_m": args.voxel_size,
        "changed_voxels": args.changed_voxels,
        "pairs": args.pairs,
        "frames_per_process": args.frames,
        "timed_frames_per_mode": args.pairs * (args.frames - 1),
        "rebuild_core_median_ms": statistics.median(pair["rebuild_core_ms"] for pair in pairs),
        "reuse_core_median_ms": statistics.median(pair["reuse_core_ms"] for pair in pairs),
        "reuse_core_no_validation_median_ms": (
            statistics.median(pair["reuse_core_no_validation_ms"] for pair in pairs)
            if not args.changed_voxels else None
        ),
        "rebuild_boundary_median_ms": statistics.median(pair["rebuild_boundary_ms"] for pair in pairs),
        "reuse_boundary_median_ms": statistics.median(pair["reuse_boundary_ms"] for pair in pairs),
        "rebuild_direct_median_ms": statistics.median(pair["rebuild_direct_ms"] for pair in pairs),
        "reuse_direct_median_ms": statistics.median(pair["reuse_direct_ms"] for pair in pairs),
        "paired_core_delta_median_ms": statistics.median(
            pair["reuse_core_ms"] - pair["rebuild_core_ms"] for pair in pairs),
        "pairs_detail": pairs,
        "values_equal": True,
    }, indent=2))


if __name__ == "__main__":
    main()
