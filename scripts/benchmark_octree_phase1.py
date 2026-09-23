#!/usr/bin/env python3
"""Timing/size comparison, with optional Phase 2 pruning and Phase 3 boundary modes."""
import argparse
import csv
import json
import math
from pathlib import Path
import statistics
import subprocess
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--binary", type=Path, default=Path("build/build_halfspace_tree"))
    parser.add_argument("--voxel-size", type=float, default=0.03)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--validation-samples", type=int, default=1000)
    parser.add_argument("--compare-pruning", action="store_true", help="compare alpha-tet and octree none/exact/rect")
    parser.add_argument("--compare-boundary", action="store_true", help="compare alpha-tet and octree rect with/without exact boundary export")
    args = parser.parse_args()
    if args.repeats < 1 or args.warmups < 0 or args.validation_samples < 0:
        parser.error("Invalid repeat/warmup/sample count")
    if not math.isfinite(args.voxel_size) or args.voxel_size <= 0:
        parser.error("voxel-size must be finite and positive")
    args.output.mkdir(parents=True, exist_ok=True)
    runs = {"alpha-tet": [], "octree-none": [], "octree-exact": [], "octree-rect": []} if args.compare_pruning else {"alpha-tet": [], "octree": []}
    if args.compare_boundary:
        if not args.compare_pruning:
            runs = {"alpha-tet": [], "octree-rect": []}
        runs["octree-boundary"] = []
    for iteration in range(-args.warmups, args.repeats):
        order = list(runs)
        shift = iteration % len(order)
        order = order[shift:] + order[:shift]
        for backend in order:
            output = args.output / backend
            command = [str(args.binary.resolve()), str(args.input.resolve()), str(output.resolve()),
                       "--pipeline", "alpha-tet" if backend == "alpha-tet" else "octree", "--voxel-size", str(args.voxel_size),
                       "--validation-samples", str(args.validation_samples)]
            if backend == "alpha-tet":
                command += ["--alpha-voxel-factor", "3", "--offset-voxel-factor", "1", "--plane-reduction", "exact"]
            else:
                pruning = "rect" if backend == "octree-boundary" else (backend.removeprefix("octree-") if backend != "octree" else "none")
                command += ["--octree-pruning", pruning]
                command += ["--octree-boundary", "exact" if backend == "octree-boundary" else "none"]
            start = time.perf_counter()
            result = subprocess.run(command, text=True, capture_output=True, timeout=180)
            wall_ms = 1000 * (time.perf_counter() - start)
            if result.returncode:
                raise RuntimeError(f"Failed: {command}\n{result.stdout}\n{result.stderr}")
            metrics = json.loads((output / "benchmark.json").read_text())[0]
            metrics["process_wall_ms"] = wall_ms
            if metrics["raw_points_outside_tree"] != 0 or metrics.get("pruning_tree_decision_mismatch", 0) != 0:
                raise RuntimeError("Benchmark encountered a missed raw point")
            if backend == "octree-boundary":
                report = json.loads((output / "validation_report.json").read_text())
                if not report["passed"] or not report["boundary_certificate"]:
                    raise RuntimeError("Benchmark encountered an invalid boundary")
            if iteration >= 0:
                runs[backend].append(metrics)
    if args.compare_boundary:
        for file in ("convex_clusters.json", "logic_tree_raw.json", "logic_tree_simplified.json"):
            if (args.output / "octree-rect" / file).read_bytes() != (args.output / "octree-boundary" / file).read_bytes():
                raise RuntimeError(f"Phase 3 changed partition/tree: {file}")
    fields = ["core_runtime_ms", "validation_ms", "total_ms", "process_wall_ms"]
    summary = {"input": str(args.input.resolve()), "voxel_size_m": args.voxel_size,
               "repeats": args.repeats, "warmups": args.warmups,
               "validation_samples": args.validation_samples,
               "compare_pruning": args.compare_pruning,
               "compare_boundary": args.compare_boundary,
               "boundary_tree_byte_consistency": True if args.compare_boundary else None,
               "notes": ["Alpha-tet and octree geometries differ; all octree modes preserve the same occupied union.",
                         "Core excludes validation and file I/O. Process wall includes both.",
                         "Octree additionally validates every voxel corner and center.",
                         "Pruned octree additionally compares the unpruned tree on the random samples.",
                         "Volume inflation uses occupied voxel union, not unknown object volume."],
               "backends": {}, "runs": runs}
    for backend, measurements in runs.items():
        last = measurements[-1]
        timings = {}
        for field in fields:
            values = sorted(m[field] for m in measurements)
            timings[field] = {"median": statistics.median(values), "min": values[0],
                              "p95": values[math.ceil(0.95 * len(values)) - 1], "max": values[-1]}
        for field in ("boundary_time", "boundary_extraction_ms", "boundary_merge_ms", "boundary_validation_ms", "boundary_io_ms"):
            if field in last:
                values = sorted(m[field] for m in measurements)
                timings[field] = {"median": statistics.median(values), "min": values[0],
                                  "p95": values[math.ceil(0.95 * len(values)) - 1], "max": values[-1]}
        sizes = {"raw_points": last["raw_points"], "occupied_voxels": last["occupied_voxels"],
                 "convex_cells": last["convex_cells"] if backend != "alpha-tet" else last["convex_clusters"],
                 "planes": last["planes_after_reduce"], "unique_planes": last["unique_planes"],
                 "tree_nodes": last["tree_nodes"], "raw_points_outside_tree": last["raw_points_outside_tree"]}
        for field in ("boundary_faces_before", "internal_faces_removed", "boundary_faces_exposed", "boundary_patches", "boundary_triangles"):
            if field in last:
                sizes[field] = last[field]
        summary["backends"][backend] = {"sizes": sizes, "timings_ms": timings}
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (args.output / "summary.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["backend", "convex_cells", "planes", "unique_planes", "tree_nodes",
                         "core_median_ms", "validation_median_ms", "total_median_ms", "wall_median_ms"])
        for backend, data in summary["backends"].items():
            sizes, timings = data["sizes"], data["timings_ms"]
            row = [backend, sizes["convex_cells"], sizes["planes"], sizes["unique_planes"], sizes["tree_nodes"]]
            row += [timings[field]["median"] for field in fields]
            writer.writerow(row)
            print(", ".join(map(str, row)))
    print(f"Summary: {args.output / 'summary.json'}")


if __name__ == "__main__":
    main()
