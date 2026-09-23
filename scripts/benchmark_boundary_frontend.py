#!/usr/bin/env python3
"""Same-process ordered-map versus contiguous boundary frontend benchmark."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import statistics
import subprocess

from benchmark_boundary_tree import REPO, execute, fixtures, save
from benchmark_octree_phase4 import write_xyz
from benchmark_validation_acceleration import UNCHANGED, compare


def end_to_end(output, binary, baseline, repeats, names):
    """Optional comparison with an executable saved before the source change."""
    rows = []
    modes = ("legacy", "contiguous")
    for name in names:
        folder = output/name/"end_to_end"; folder.mkdir()
        records = {mode: [] for mode in modes}
        for mode in modes:
            (folder/mode).mkdir()
        for trial in range(-1, repeats):
            for mode in (modes if trial % 2 == 0 else modes[::-1]):
                out = folder/mode
                command = [str(baseline if mode == "legacy" else binary), str(output/name/"input.xyz"), str(out),
                    "--pipeline", "octree", "--tree-source", "boundary", "--boundary-expression", "direct",
                    "--boundary-query", "dag", "--voxel-size", ".04", "--validation-samples", "1000",
                    "--distance-field", "lse", "--lse-error", ".001", "--execution-sharing", "structural",
                    "--lse-kernel", "binary", "--validation-cache", "auto", "--sign-propagation", "packed"]
                execute(command, out/f"trial_{trial}.log")
                metrics = json.loads((out/"benchmark.json").read_text())[0]
                assert metrics["passed"] and metrics["strict_sign_certified"]
                if trial >= 0:
                    records[mode].append(metrics)
            compare(folder/"legacy", folder/"contiguous")
            # Also retain all deterministic work, cache, strata and error
            # counters. Only elapsed times may differ in this comparison.
            deterministic = []
            for mode in modes:
                metrics = json.loads((folder/mode/"benchmark.json").read_text())[0]
                deterministic.append({k: v for k, v in metrics.items() if not k.endswith("_ms") and k != "core_build_time"})
            assert deterministic[0] == deterministic[1], (name, "non-timing metrics changed")
        row = {"case": name, "same_outputs": True, "same_non_timing_metrics": True}
        for mode in modes:
            save(folder/mode/"trials.json", records[mode])
            row[mode] = {key: statistics.median(t[key] for t in records[mode]) for key in (
                "boundary_merge_ms", "registry_ms", "geometry_core_ms", "certificate_ms", "validation_ms", "total_ms")}
        row["artifact_sha256"] = {filename: hashlib.sha256((folder/"contiguous"/filename).read_bytes()).hexdigest()
            for filename in UNCHANGED}
        rows.append(row)
    save(output/"end_to_end_summary.json", rows)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--build-dir", type=Path, default=REPO/"build")
    parser.add_argument("--repeats", type=int, default=21)
    parser.add_argument("--baseline-binary", type=Path, help="Optional pre-change build_halfspace_tree for end-to-end output comparison")
    parser.add_argument("--pipeline-repeats", type=int, default=5)
    parser.add_argument("--cases", nargs="+", default=("single", "solid_cuboid", "l_voxels", "u_voxels",
        "staircase", "closed_cavity", "two_disconnected", "thin_wall", "narrow_corridor", "camera",
        "large_slab", "large_solid", "checkerboard", "sparse_disconnected"))
    args = parser.parse_args()
    if not 1 <= args.repeats <= 10000 or args.pipeline_repeats < 1:
        parser.error("Repeats must be in [1,10000]")
    cases = dict(fixtures())
    extra = {
        "large_slab": itertools.product(range(96), range(96), range(1)),
        "large_solid": itertools.product(range(24), repeat=3),
        "checkerboard": (p for p in itertools.product(range(20), repeat=3) if sum(p) % 2 == 0),
    }
    for name, keys in extra.items():
        cases[name] = [tuple(.04*(v+.5) for v in key) for key in keys]
    if set(args.cases)-cases.keys():
        parser.error("Unknown case")
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=False)
    binary = args.build_dir.resolve()/"boundary_frontend_benchmark"
    manifest = {
        "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "sources_sha256": {name: hashlib.sha256((REPO/name).read_bytes()).hexdigest() for name in (
            "src/octree_boundary.cpp", "src/octree_boundary_tree.cpp", "tests/boundary_frontend_reference.hpp",
            "tests/boundary_frontend_benchmark.cpp", "tests/boundary_frontend_checks.hpp")},
        "repeats": args.repeats, "warmups": 2, "voxel_size_m": .04,
        "timing": "same process; alternating legacy/contiguous pairs; median milliseconds; checks outside timing",
        "frontend_scope": "exposed-face extraction, patch merge and support registry; temporary container teardown included; excludes voxelization, certification, tree construction, queries and I/O",
        "merge_scope": "internal merge_ms excludes temporary container teardown in both implementations",
        "correctness": "every ordered face, owner, patch rectangle, source face ID, oriented corner, support coordinate, sign mask and patch ID matched after every pair; both boundaries independently certified",
        "inputs": [],
    }
    if args.baseline_binary:
        args.baseline_binary = args.baseline_binary.resolve()
        manifest["end_to_end"] = {
            "baseline_binary": str(args.baseline_binary),
            "baseline_sha256": hashlib.sha256(args.baseline_binary.read_bytes()).hexdigest(),
            "current_sha256": hashlib.sha256((args.build_dir.resolve()/"build_halfspace_tree").read_bytes()).hexdigest(),
            "repeats": args.pipeline_repeats, "warmups": 1,
            "protocol": "separate processes in alternating pairs; 1 mm LSE budget, 1000 random validation samples; all deterministic artifacts and non-timing metrics checked after every pair",
        }
    rows = []
    for name in args.cases:
        folder = output/name; folder.mkdir(); cloud = folder/"input.xyz"; write_xyz(cloud, cases[name])
        command = [str(binary), str(cloud), ".04", str(args.repeats)]
        save(folder/"command.json", command)
        run = subprocess.run(command, check=True, capture_output=True, text=True, timeout=120)
        row = json.loads(run.stdout); row["case"] = name; row["passed"] = True
        for metric in ("merge_ms", "registry_ms", "frontend_ms"):
            row[metric+"_reduction"] = 1-row["contiguous"][metric]/row["legacy"][metric]
        rows.append(row); save(folder/"trials.json", row)
        manifest["inputs"].append({"case": name, "sha256": hashlib.sha256(cloud.read_bytes()).hexdigest()})
        save(output/"manifest.json", manifest); save(output/"summary.json", rows)
    report = ["# Boundary frontend containers", "",
        f"h=.04 m; 2 warmup pairs + {args.repeats} measured alternating pairs per scene in one process.",
        "Times are medians in milliseconds; all equality checks and certificates are outside timing.", "",
        "| Scene | Voxels / faces / patches | Legacy / contiguous merge ms | Reduction | Legacy / contiguous registry ms | Legacy / contiguous frontend ms | Reduction |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in rows:
        a, b = row["legacy"], row["contiguous"]
        report.append(f'| {row["case"]} | {row["occupied_voxels"]} / {row["faces"]} / {row["patches"]} | '
            f'{a["merge_ms"]:.6f} / {b["merge_ms"]:.6f} | {100*row["merge_ms_reduction"]:.1f}% | '
            f'{a["registry_ms"]:.6f} / {b["registry_ms"]:.6f} | '
            f'{a["frontend_ms"]:.6f} / {b["frontend_ms"]:.6f} | {100*row["frontend_ms_reduction"]:.1f}% |')
    report += ["", "Frontend totals include face extraction and scratch-container destruction. Internal merge timers exclude scratch-container destruction in both modes.",
        "Voxelization, tree construction, certification, queries, final-result destruction and file I/O are excluded from these frontend totals.",
        "Face order, owners, rectangle tiling, source face order and support registry matched exactly in every pair; both boundaries passed full face/patch coverage validation.",
        "The legacy ordered-map implementation is frozen in tests/boundary_frontend_reference.hpp; production has no extra mode or cache configuration.",
        "Raw timing arrays and source/input/binary hashes are retained in JSON. Tiny-scene sub-microsecond timings are sensitive to noise."]
    if args.baseline_binary:
        names = [name for name in args.cases if name in ("l_voxels", "closed_cavity", "camera")]
        pipeline = end_to_end(output, args.build_dir.resolve()/"build_halfspace_tree", args.baseline_binary,
            args.pipeline_repeats, names)
        report += ["", "## End-to-end comparison with the saved pre-change executable", "",
            "Direct DAG, 1 mm LSE budget, full strata certification and distance audit. Median milliseconds.", "",
            "| Scene | Legacy / contiguous geometry core ms | Legacy / contiguous validation ms | Legacy / contiguous total ms | Exact outputs and non-timing counters |",
            "| --- | ---: | ---: | ---: | --- |"]
        for row in pipeline:
            a, b = row["legacy"], row["contiguous"]
            report.append(f'| {row["case"]} | {a["geometry_core_ms"]:.6f} / {b["geometry_core_ms"]:.6f} | '
                f'{a["validation_ms"]:.6f} / {b["validation_ms"]:.6f} | {a["total_ms"]:.6f} / {b["total_ms"]:.6f} | yes |')
        report += ["", "The pipeline total_ms covers construction and validation, excluding input/output I/O and optional query benchmarks. It is dominated by unchanged validation and is sensitive to process/system noise.",
            "All tree, geometry, support, split, occupancy, sign, sampled distance/gradient and ray-validation outputs matched; all non-timing benchmark counters were equal."]
    (output/"report.md").write_text("\n".join(report)+"\n")
    print(output/"report.md")


if __name__ == "__main__":
    main()
