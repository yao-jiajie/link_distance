#!/usr/bin/env python3
"""Generate a deterministic nonconvex C-bracket with curved and planar surfaces."""

import argparse
import json
import math
from pathlib import Path
import random


REPO = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = REPO / "data/c_bracket_surface.xyz"


def proportional_counts(total, weighted_surfaces):
    area_sum = sum(area for _, area in weighted_surfaces)
    exact = [total * area / area_sum for _, area in weighted_surfaces]
    counts = [math.floor(value) for value in exact]
    remainder_order = sorted(
        range(len(counts)), key=lambda index: exact[index] - counts[index], reverse=True)
    for index in remainder_order[:total - sum(counts)]:
        counts[index] += 1
    return counts


def generate(count, seed):
    center = (0.0, 0.0, 0.74)
    outer_radius = 0.22
    inner_radius = 0.10
    height = 0.24
    gap_angle = math.radians(70.0)
    theta_begin = gap_angle * 0.5
    theta_end = 2.0 * math.pi - gap_angle * 0.5
    theta_span = theta_end - theta_begin

    annular_area = 0.5 * (outer_radius**2 - inner_radius**2) * theta_span
    cut_area = (outer_radius - inner_radius) * height
    surfaces = [
        ("outer_curved_wall", outer_radius * theta_span * height),
        ("inner_curved_wall", inner_radius * theta_span * height),
        ("top_planar_face", annular_area),
        ("bottom_planar_face", annular_area),
        ("cut_planar_face_begin", cut_area),
        ("cut_planar_face_end", cut_area),
    ]
    counts = proportional_counts(count, surfaces)
    rng = random.Random(seed)
    points = []

    def point(radius, theta, z):
        return (center[0] + radius * math.cos(theta),
                center[1] + radius * math.sin(theta), center[2] + z)

    for _ in range(counts[0]):
        points.append(point(outer_radius, rng.uniform(theta_begin, theta_end),
                            rng.uniform(-height * 0.5, height * 0.5)))
    for _ in range(counts[1]):
        points.append(point(inner_radius, rng.uniform(theta_begin, theta_end),
                            rng.uniform(-height * 0.5, height * 0.5)))
    for surface_index, z in ((2, height * 0.5), (3, -height * 0.5)):
        for _ in range(counts[surface_index]):
            radius = math.sqrt(inner_radius**2 + rng.random() *
                               (outer_radius**2 - inner_radius**2))
            points.append(point(radius, rng.uniform(theta_begin, theta_end), z))
    for surface_index, theta in ((4, theta_begin), (5, theta_end)):
        for _ in range(counts[surface_index]):
            points.append(point(rng.uniform(inner_radius, outer_radius), theta,
                                rng.uniform(-height * 0.5, height * 0.5)))
    rng.shuffle(points)

    metadata = {
        "schema": 1,
        "name": "c_bracket_surface",
        "description": "Extruded nonconvex C-bracket with cylindrical walls and planar caps/cut faces",
        "length_unit": "m",
        "seed": seed,
        "raw_points": count,
        "sampling": "analytic boundary sampled proportional to surface area",
        "shape": {
            "type": "extruded_annular_sector",
            "center_m": list(center),
            "outer_radius_m": outer_radius,
            "inner_radius_m": inner_radius,
            "height_m": height,
            "gap_angle_deg": math.degrees(gap_angle),
        },
        "surface_point_counts": {name: n for (name, _), n in zip(surfaces, counts)},
        "semantics": (
            "Points sample the complete analytic boundary. The downstream sparse voxel pipeline "
            "occupies only voxels containing samples and does not fill the analytic solid interior."
        ),
        "bounds_m": [
            [min(sample[axis] for sample in points) for axis in range(3)],
            [max(sample[axis] for sample in points) for axis in range(3)],
        ],
    }
    return points, metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--points", type=int, default=4096)
    parser.add_argument("--seed", type=int, default=20260915)
    args = parser.parse_args()
    if args.points < 6:
        parser.error("points must be at least six")
    output = args.output.resolve()
    points, metadata = generate(args.points, args.seed)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("".join(f"{x:.17g} {y:.17g} {z:.17g}\n" for x, y, z in points))
    metadata_path = output.with_suffix(".json")
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"xyz": str(output), "metadata": str(metadata_path),
                      "points": len(points)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
