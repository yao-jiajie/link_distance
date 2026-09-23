#!/usr/bin/env python3
"""Generate two deterministic, disconnected, nonconvex surface point clouds."""

import argparse
import json
import math
from pathlib import Path
import random

from generate_c_bracket_surface import generate as generate_c_bracket


REPO = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = REPO / "data/disconnected_irregular_surface.xyz"


def generate(count, seed):
    if count < 512:
        raise ValueError("point count must be at least 512")

    bracket_count = round(count * 0.5625)
    torus_count = count - bracket_count

    bracket_source, bracket_metadata = generate_c_bracket(bracket_count, seed)
    bracket_shift = (-0.29, 0.0, 0.0)
    bracket = [tuple(point[axis] + bracket_shift[axis] for axis in range(3))
               for point in bracket_source]

    rng = random.Random(seed + 1)
    torus_center = (0.31, 0.0, 0.75)
    major_radius = 0.135
    minor_radius = 0.045
    torus = []
    for _ in range(torus_count):
        theta = rng.uniform(0.0, 2.0 * math.pi)
        phi = rng.uniform(0.0, 2.0 * math.pi)
        # A three-lobed radius and a vertically waving centerline make this
        # component visibly irregular while keeping a toroidal central hole.
        local_major = major_radius * (1.0 + 0.12 * math.cos(3.0 * theta))
        local_minor = minor_radius * (1.0 + 0.10 * math.sin(theta))
        radial = local_major + local_minor * math.cos(phi)
        torus.append((
            torus_center[0] + radial * math.cos(theta),
            torus_center[1] + radial * math.sin(theta),
            torus_center[2] + 0.028 * math.sin(2.0 * theta + 0.4)
            + local_minor * math.sin(phi),
        ))

    points = [(point, "c_bracket") for point in bracket]
    points += [(point, "wavy_torus") for point in torus]
    random.Random(seed + 2).shuffle(points)
    shuffled = [point for point, _ in points]

    metadata = {
        "schema": 1,
        "name": "disconnected_irregular_surface",
        "description": "Disconnected C-bracket and three-lobed wavy torus surface samples",
        "length_unit": "m",
        "seed": seed,
        "raw_points": len(shuffled),
        "sampling": "complete analytic surfaces; deterministic parameter sampling",
        "components": [
            {
                "name": "c_bracket",
                "points": bracket_count,
                "translation_m": list(bracket_shift),
                "shape": bracket_metadata["shape"],
            },
            {
                "name": "wavy_torus",
                "points": torus_count,
                "shape": {
                    "type": "three_lobed_wavy_torus",
                    "center_m": list(torus_center),
                    "major_radius_m": major_radius,
                    "minor_radius_m": minor_radius,
                    "radial_modulation_fraction": 0.12,
                    "vertical_wave_amplitude_m": 0.028,
                },
            },
        ],
        "expected_free_probes_m": [[0.0, 0.0, 0.75]],
        "semantics": (
            "Both components are surface samples. The downstream sparse voxel pipeline "
            "occupies only voxels containing samples and does not fill analytic interiors."
        ),
        "bounds_m": [
            [min(point[axis] for point in shuffled) for axis in range(3)],
            [max(point[axis] for point in shuffled) for axis in range(3)],
        ],
    }
    return shuffled, metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--points", type=int, default=4096)
    parser.add_argument("--seed", type=int, default=20260916)
    args = parser.parse_args()
    try:
        points, metadata = generate(args.points, args.seed)
    except ValueError as error:
        parser.error(str(error))

    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("".join(f"{x:.17g} {y:.17g} {z:.17g}\n" for x, y, z in points))
    metadata_path = output.with_suffix(".json")
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({
        "xyz": str(output),
        "metadata": str(metadata_path),
        "points": len(points),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
