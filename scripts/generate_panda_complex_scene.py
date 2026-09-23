#!/usr/bin/env python3
"""Generate a denser six-shape XYZ scene around the unchanged Panda pose."""

import argparse
import json
import math
from pathlib import Path
import random

from generate_panda_disconnected_scene import generate as generate_three_shapes


REPO = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = REPO / "data/panda_complex_scene.xyz"
DEFAULT_SEED = 20260922


def disturbed(point, rng, standard_deviation=0.001):
    return tuple(value + rng.gauss(0.0, standard_deviation) for value in point)


def twisted_ribbon(count, rng):
    points = []
    for _ in range(count):
        u = rng.uniform(-1.0, 1.0)
        v = rng.uniform(-1.0, 1.0)
        point = (
            -0.62 + 0.16 * u,
            0.08 + 0.10 * v + 0.025 * math.sin(3.0 * math.pi * u),
            0.72 + 0.12 * u + 0.055 * v * math.sin(2.0 * math.pi * u),
        )
        points.append(disturbed(point, rng))
    return points, {
        "type": "twisted_corrugated_ribbon",
        "center_m": [-0.62, 0.08, 0.72],
        "parameter_domain": "u,v in [-1,1]",
    }


def fluted_column(count, rng):
    wall_count = count - count // 10
    points = []
    for _ in range(wall_count):
        t = rng.random()
        theta = rng.uniform(0.0, 2.0 * math.pi)
        radius = 0.075 * (1.0 + 0.18 * math.sin(7.0 * theta +
                                                 2.0 * math.pi * t)
                          + 0.08 * math.sin(4.0 * math.pi * t))
        point = (0.15 + radius * math.cos(theta),
                 -0.72 + radius * math.sin(theta),
                 0.22 + 0.55 * t)
        points.append(disturbed(point, rng))
    for index in range(count - wall_count):
        theta = rng.uniform(0.0, 2.0 * math.pi)
        radius = 0.075 * math.sqrt(rng.random())
        z = 0.22 if index % 2 == 0 else 0.77
        points.append(disturbed((0.15 + radius * math.cos(theta),
                                 -0.72 + radius * math.sin(theta), z), rng))
    return points, {
        "type": "capped_seven_flute_column",
        "center_xy_m": [0.15, -0.72],
        "height_interval_m": [0.22, 0.77],
        "nominal_radius_m": 0.075,
    }


def folded_canopy(count, rng):
    points = []
    for _ in range(count):
        u = rng.uniform(-1.0, 1.0)
        v = rng.uniform(-1.0, 1.0)
        point = (
            0.86 + 0.16 * u,
            0.50 + 0.12 * v,
            0.72 + 0.10 * u * u + 0.05 * math.cos(3.0 * math.pi * v)
            + 0.035 * math.sin(4.0 * math.pi * u * v),
        )
        points.append(disturbed(point, rng))
    return points, {
        "type": "folded_nonconvex_canopy",
        "center_m": [0.86, 0.50, 0.72],
        "parameter_domain": "u,v in [-1,1]",
    }


def generate(count=9216, seed=DEFAULT_SEED):
    if count < 3840:
        raise ValueError("point count must be at least 3840")
    first_count = count // 2
    extra_count = count - first_count
    counts = (extra_count // 3, extra_count // 3,
              extra_count - 2 * (extra_count // 3))
    first_points, first_metadata = generate_three_shapes(first_count, seed)
    new_shapes = [
        ("twisted_ribbon", twisted_ribbon, counts[0]),
        ("fluted_column", fluted_column, counts[1]),
        ("folded_canopy", folded_canopy, counts[2]),
    ]
    labelled = [(point, "original_three") for point in first_points]
    components = list(first_metadata["components"])
    for index, (name, generator, shape_count) in enumerate(new_shapes):
        points, shape = generator(shape_count, random.Random(seed + 10 + index))
        components.append({"name": name, "points": len(points), "shape": shape})
        labelled.extend((point, name) for point in points)
    random.Random(seed + 20).shuffle(labelled)
    points = [point for point, _ in labelled]
    metadata = {
        "schema": 1,
        "name": "panda_complex_scene",
        "description": "Six irregular obstacle surfaces around the fixed Panda pose",
        "length_unit": "m",
        "seed": seed,
        "raw_points": len(points),
        "recommended_voxel_size_m": 0.02,
        "recommended_configuration": [0, -0.4, 0, -2.0, 0, 1.6, 0.8],
        "recommended_field_from_base": [0, 0, 0, 0, 0, 0],
        "new_shape_noise_standard_deviation_m": 0.001,
        "components": components,
        "sampling": "Deterministic analytic surfaces; new shapes have 1 mm Gaussian coordinate noise",
        "semantics": (
            "One XYZ file contains all six shapes. The voxel pipeline occupies "
            "sample-containing voxels only; it does not fill analytic interiors."
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
    parser.add_argument("--points", type=int, default=9216)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()
    try:
        points, metadata = generate(args.points, args.seed)
    except ValueError as error:
        parser.error(str(error))
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("".join(f"{x:.17g} {y:.17g} {z:.17g}\n"
                              for x, y, z in points))
    metadata_path = output.with_suffix(".json")
    metadata_path.write_text(json.dumps(metadata, indent=2,
                                        ensure_ascii=False) + "\n")
    print(json.dumps({"xyz": str(output), "metadata": str(metadata_path),
                      "points": len(points),
                      "components": len(metadata["components"])},
                     ensure_ascii=False))


if __name__ == "__main__":
    main()
