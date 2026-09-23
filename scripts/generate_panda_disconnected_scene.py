#!/usr/bin/env python3
"""Generate three disconnected irregular obstacle surfaces around the fixed Panda pose."""

import argparse
import json
import math
from pathlib import Path
import random

from generate_c_bracket_surface import generate as generate_c_bracket


REPO = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = REPO / "data/panda_disconnected_irregular_scene.xyz"


def add(a, b):
    return tuple(a[index] + b[index] for index in range(3))


def scale(a, value):
    return tuple(component * value for component in a)


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def normalized(a):
    length = math.sqrt(sum(component * component for component in a))
    return tuple(component / length for component in a)


def wavy_torus(count, rng):
    center = (0.80, -0.38, 0.70)
    major_radius = 0.135
    minor_radius = 0.045
    points = []
    for _ in range(count):
        theta = rng.uniform(0.0, 2.0 * math.pi)
        phi = rng.uniform(0.0, 2.0 * math.pi)
        local_major = major_radius * (1.0 + 0.14 * math.cos(3.0 * theta))
        local_minor = minor_radius * (1.0 + 0.12 * math.sin(2.0 * theta))
        radial = local_major + local_minor * math.cos(phi)
        points.append((
            center[0] + radial * math.cos(theta),
            center[1] + radial * math.sin(theta),
            center[2] + 0.026 * math.sin(2.0 * theta + 0.3)
            + local_minor * math.sin(phi),
        ))
    return points, {
        "type": "three_lobed_wavy_torus",
        "center_m": list(center),
        "major_radius_m": major_radius,
        "minor_radius_m": minor_radius,
    }


def tube_frame(parameter):
    x = -0.05 + 0.55 * parameter
    y = 0.43 + 0.055 * math.sin(2.0 * math.pi * parameter)
    z = (0.28 + 0.13 * math.sin(math.pi * parameter)
         + 0.035 * math.sin(5.0 * math.pi * parameter))
    tangent = normalized((
        0.55,
        0.055 * 2.0 * math.pi * math.cos(2.0 * math.pi * parameter),
        0.13 * math.pi * math.cos(math.pi * parameter)
        + 0.035 * 5.0 * math.pi * math.cos(5.0 * math.pi * parameter),
    ))
    normal = normalized(cross(tangent, (0.0, 0.0, 1.0)))
    binormal = normalized(cross(tangent, normal))
    return (x, y, z), tangent, normal, binormal


def bent_tube(count, rng):
    wall_count = count - 2 * max(1, count // 24)
    cap_count = (count - wall_count) // 2
    points = []
    for _ in range(wall_count):
        parameter = rng.random()
        center, _, normal, binormal = tube_frame(parameter)
        phi = rng.uniform(0.0, 2.0 * math.pi)
        radius = 0.047 * (1.0 + 0.16 * math.sin(6.0 * math.pi * parameter))
        radial = add(scale(normal, math.cos(phi)), scale(binormal, math.sin(phi)))
        points.append(add(center, scale(radial, radius)))
    for parameter in (0.0, 1.0):
        center, _, normal, binormal = tube_frame(parameter)
        for _ in range(cap_count):
            radius = 0.047 * math.sqrt(rng.random())
            phi = rng.uniform(0.0, 2.0 * math.pi)
            radial = add(scale(normal, math.cos(phi)), scale(binormal, math.sin(phi)))
            points.append(add(center, scale(radial, radius)))
    while len(points) < count:
        points.append(tube_frame(1.0)[0])
    return points, {
        "type": "capped_variable_radius_bent_tube",
        "centerline_start_m": list(tube_frame(0.0)[0]),
        "centerline_end_m": list(tube_frame(1.0)[0]),
        "nominal_radius_m": 0.047,
    }


def generate(count, seed):
    if count < 768:
        raise ValueError("point count must be at least 768")
    bracket_count = round(count * 0.375)
    torus_count = round(count * 0.3125)
    tube_count = count - bracket_count - torus_count

    bracket_source, bracket_metadata = generate_c_bracket(bracket_count, seed)
    bracket_translation = (-0.48, -0.42, -0.26)
    bracket = [add(point, bracket_translation) for point in bracket_source]
    torus, torus_shape = wavy_torus(torus_count, random.Random(seed + 1))
    tube, tube_shape = bent_tube(tube_count, random.Random(seed + 2))

    labelled = ([(point, "c_bracket") for point in bracket]
                + [(point, "wavy_torus") for point in torus]
                + [(point, "bent_tube") for point in tube])
    random.Random(seed + 3).shuffle(labelled)
    points = [point for point, _ in labelled]
    metadata = {
        "schema": 1,
        "name": "panda_disconnected_irregular_scene",
        "description": (
            "Three disconnected irregular surface clouds placed around the fixed "
            "Panda example pose without moving the robot"
        ),
        "length_unit": "m",
        "seed": seed,
        "raw_points": len(points),
        "recommended_voxel_size_m": 0.04,
        "recommended_configuration": [0, -0.4, 0, -2.0, 0, 1.6, 0.8],
        "recommended_field_from_base": [0, 0, 0, 0, 0, 0],
        "components": [
            {
                "name": "c_bracket",
                "points": len(bracket),
                "translation_m": list(bracket_translation),
                "shape": bracket_metadata["shape"],
            },
            {"name": "wavy_torus", "points": len(torus), "shape": torus_shape},
            {"name": "bent_tube", "points": len(tube), "shape": tube_shape},
        ],
        "sampling": "complete deterministic analytic surface samples",
        "semantics": (
            "All components share one XYZ input. The voxel pipeline occupies only "
            "sample-containing voxels and does not fill analytic interiors."
        ),
        "bounds_m": [
            [min(point[axis] for point in points) for axis in range(3)],
            [max(point[axis] for point in points) for axis in range(3)],
        ],
    }
    return points, metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--points", type=int, default=4608)
    parser.add_argument("--seed", type=int, default=20260921)
    args = parser.parse_args()
    try:
        points, metadata = generate(args.points, args.seed)
    except ValueError as error:
        parser.error(str(error))

    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("".join(
        f"{x:.17g} {y:.17g} {z:.17g}\n" for x, y, z in points))
    metadata_path = output.with_suffix(".json")
    metadata_path.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({
        "xyz": str(output),
        "metadata": str(metadata_path),
        "points": len(points),
        "components": len(metadata["components"]),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
