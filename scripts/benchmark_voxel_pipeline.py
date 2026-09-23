"""Sequential, reproducible voxel/raw/Jet comparisons with matched requested scales."""
import argparse
import csv
import json
import math
from pathlib import Path
import statistics
import subprocess
import tempfile
import time


def execute(command, directory, timeout):
    start = time.perf_counter()
    result = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
    wall_ms = (time.perf_counter() - start) * 1000
    (directory / "command.json").write_text(json.dumps(command, indent=2) + "\n")
    (directory / "run.log").write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f"Run failed ({result.returncode}), see {directory / 'run.log'}")
    record = json.loads((directory / "benchmark.json").read_text())
    metrics = record[0] if isinstance(record, list) else record
    metrics["process_wall_ms"] = wall_ms
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path, help="new directory; existing results are not overwritten")
    parser.add_argument("--build-dir", type=Path, default=Path("build"))
    parser.add_argument("--voxel-sizes", type=float, nargs="+", default=[0.01, 0.02, 0.03])
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--validation-samples", type=int, default=1000)
    parser.add_argument("--with-nef", action="store_true", help="also compare exact-nef on each voxel wrap")
    parser.add_argument("--timeout", type=float, default=120.0)
    args = parser.parse_args()
    if (args.repeats < 1 or args.validation_samples < 0 or
            any(not math.isfinite(s) or s <= 0 for s in args.voxel_sizes) or
            not math.isfinite(args.timeout) or args.timeout <= 0):
        parser.error("repeats and voxel sizes must be positive; samples must be nonnegative")
    tree = (args.build_dir / "build_halfspace_tree").resolve()
    voxelize = (args.build_dir / "voxelize_point_cloud").resolve()
    source = args.input.resolve()
    if not source.is_file() or not tree.is_file() or not voxelize.is_file():
        parser.error("input and built executables must exist")
    args.output.mkdir(parents=True, exist_ok=False)
    summary = []
    for size in args.voxel_sizes:
        label = f"voxel_{size:g}m"
        stage = args.output / label / "voxelization"
        subprocess.run([str(voxelize), str(source), str(stage), "--voxel-size", str(size)],
                       check=True, timeout=args.timeout)
        for mode in ("voxel", "none", "jet"):
            destination = args.output / label / mode
            destination.mkdir(parents=True)
            parameters = (["--voxel-size", str(size)] if mode == "voxel" else
                          ["--preprocess", mode, "--alpha", str(3 * size), "--offset", str(size)])
            common = ["--pipeline", "alpha-tet", "--plane-reduction", "exact",
                      "--validation-samples", str(args.validation_samples)] + parameters
            trials = []
            for repeat in range(args.repeats):
                if repeat == 0:
                    trials.append(execute([str(tree), str(source), str(destination)] + common,
                                          destination, args.timeout))
                else:
                    with tempfile.TemporaryDirectory(prefix="cgal_voxel_repeat_") as temporary:
                        path = Path(temporary)
                        trials.append(execute([str(tree), str(source), str(path)] + common,
                                              path, args.timeout))
            # Geometry should not depend on runtime; fail instead of averaging
            # different wraps or silently hiding a validation failure.
            for trial in trials:
                for key in ("wrap_faces", "tetra_count", "convex_clusters", "tree_nodes", "offset"):
                    if trial[key] != trials[0][key]:
                        raise RuntimeError(f"Non-reproducible geometry: {destination}, {key}")
            record = dict(trials[0])
            record.update({key: statistics.median(t[key] for t in trials)
                           for key in record if key.endswith("_ms")})
            record.update(case=label, mode=mode, matched_scale_m=size,
                          repeats=args.repeats, validation_samples=args.validation_samples, status="passed")
            record["core_min_ms"] = min(t["core_runtime_ms"] for t in trials)
            record["core_max_ms"] = max(t["core_runtime_ms"] for t in trials)
            (destination / "trials.json").write_text(json.dumps(trials, indent=2) + "\n")
            summary.append(record)
            print(f"{label} {mode}: points={record['point_count']}->{record['voxel_points'] or record['point_count']} "
                  f"faces={record['wrap_faces']} clusters={record['convex_clusters']} "
                  f"nodes={record['tree_nodes']} core_median={record['core_runtime_ms']:.3f} ms", flush=True)
        if args.with_nef:
            destination = args.output / label / "exact_nef"
            destination.mkdir()
            source_wrap = args.output / label / "voxel" / "wrap.off"
            voxel_record = summary[-3]
            command = [str(tree), str(source_wrap), str(destination), "--pipeline", "exact-nef",
                       "--validation-samples", str(args.validation_samples),
                       "--point-count", str(voxel_record["raw_points"]),
                       "--alpha", str(voxel_record["alpha"]), "--offset", str(voxel_record["offset"])]
            try:
                record = execute(command, destination, args.timeout)
                record.update(case=label, mode="exact-nef", matched_scale_m=size,
                              repeats=1, validation_samples=args.validation_samples, status="passed")
                record["convex_clusters"] = record["convex_parts"]
                # Source preprocessing + wrapping are shared with alpha-tet.
                # These are named stage-sums, not fictitious end-to-end timings.
                record["shared_preprocess_wrap_ms"] = (
                    voxel_record["preprocessing_ms"] + voxel_record["alpha_wrap_ms"])
                record["shared_prefix_plus_nef_ms"] = record["shared_preprocess_wrap_ms"] + record["total_ms"]
                record["timing_scope"] = "Nef total starts at OFF input; prefix sum excludes OFF writing and oracle setup"
                summary.append(record)
                print(f"{label} exact-nef: clusters={record['convex_parts']} "
                      f"nodes={record['tree_nodes']} OFF-to-tree={record['total_ms']:.3f} ms", flush=True)
            except subprocess.TimeoutExpired:
                (destination / "command.json").write_text(json.dumps(command, indent=2) + "\n")
                summary.append(dict(case=label, mode="exact-nef", status="timeout", timeout_s=args.timeout))
                print(f"{label} exact-nef: exceeded {args.timeout:g} s", flush=True)
            except RuntimeError as error:
                # Legacy failures are comparison results, not successful trees.
                # Keep the new pipeline metrics and the complete error log.
                summary.append(dict(case=label, mode="exact-nef", status="failed",
                                    error=str(error), log=str(destination / "run.log")))
                print(f"{label} exact-nef failed; see {destination / 'run.log'}", flush=True)
        # Save progress after each size, including a possible Nef timeout.
        (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        columns = list(dict.fromkeys(key for row in summary for key in row))
        with (args.output / "summary.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=columns)
            writer.writeheader()
            writer.writerows(summary)


if __name__ == "__main__":
    main()
