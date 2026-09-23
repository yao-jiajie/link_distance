#!/usr/bin/env python3
"""Sequential multi-scene accuracy/latency benchmark; never changes a backend."""
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import random
import statistics
import subprocess
import sys
import time

from benchmark_scene_clouds import SCENES, box_volume, bounds, interior_points, noisy_points, solid_volume, surface_points

REPO = Path(__file__).resolve().parents[1]
BACKENDS = ("alpha-tet", "octree-rect")


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_xyz(path):
    points = []
    for number, line in enumerate(path.read_text().splitlines(), 1):
        fields = line.split("#", 1)[0].split()
        if not fields:
            continue
        if len(fields) != 3:
            raise ValueError(f"{path}:{number}: expected exactly three meter XYZ columns")
        point = tuple(map(float, fields))
        if not all(math.isfinite(p) for p in point):
            raise ValueError(f"{path}:{number}: nonfinite coordinate")
        points.append(point)
    if not points:
        raise ValueError(f"Empty cloud: {path}")
    return points


def write_xyz(path, points):
    path.write_text("".join(" ".join(format(x, ".17g") for x in p)+"\n" for p in points))


def stats(values):
    values = sorted(values)
    return {"min": values[0], "median": statistics.median(values),
            "p95": values[math.ceil(.95*len(values))-1], "max": values[-1]}


def execute(command, log, timeout):
    write_json(log.with_suffix(".command.json"), command)
    start = time.perf_counter()
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired as error:
        log.write_text(f"TIMEOUT after {timeout} seconds\n" + str(error.stdout or "") + str(error.stderr or ""))
        raise RuntimeError(f"Timeout; see {log}") from error
    wall_ms = 1000*(time.perf_counter()-start)
    log.write_text(result.stdout+result.stderr)
    if result.returncode:
        raise RuntimeError(f"Exit {result.returncode}; see {log}")
    return wall_ms


def make_clouds(args):
    folder = args.output / "inputs"
    folder.mkdir()
    clouds = []
    for name in args.scenes:
        scene = SCENES[name]
        seed = args.seed + list(SCENES).index(name)*100
        clean = surface_points(scene, max(args.densities), seed)
        for count, noise in itertools.product(args.densities, args.noise_levels):
            label = f"{name}_n{count}_noise{noise:g}"
            points = noisy_points(clean[:count], noise, seed+1)
            path = folder / f"{label}.xyz"
            write_xyz(path, points)
            clouds.append({"name": label, "scene": name, "provenance": "synthetic_surface",
                           "points": points, "path": path, "noise_halfwidth_m": noise,
                           "truth": scene, "surface_seed": seed, "noise_seed": seed+1})
    if not args.skip_camera:
        path = REPO / "data/realistic_sparse_camera.xyz"
        metadata = json.loads((REPO / "data/realistic_sparse_camera_metadata.json").read_text())
        write_json(folder / "camera_metadata.json", metadata)
        destination = folder / "camera_simulated.xyz"
        points = read_xyz(path)
        write_xyz(destination, points)
        clouds.append({"name": "camera_simulated", "scene": "camera", "provenance": "simulated_multiview_RGBD",
                       "path": destination, "points": points, "noise_halfwidth_m": None,
                       "truth": {"center": metadata["object"]["center_m"], "radii": metadata["object"]["radii_m"]}})
    for i, path in enumerate(args.measured_input):
        points = read_xyz(path)
        destination = folder / f"measured_{i}.xyz"
        write_xyz(destination, points)
        clouds.append({"name": f"measured_{i}", "scene": "measured", "provenance": "user_supplied_measured_XYZ",
                       "path": destination, "points": points, "noise_halfwidth_m": None, "truth": None,
                       "original_path": str(path.resolve())})
    manifest = []
    for cloud in clouds:
        item = {k: str(v) if isinstance(v, Path) else v for k, v in cloud.items() if k != "points"}
        item.update(raw_points=len(cloud["points"]), sha256=digest(cloud["path"]), length_unit="m")
        manifest.append(item)
    write_json(folder / "manifest.json", manifest)
    return clouds


def build_case(args, cloud, size, folder, backends):
    records = {b: {"status": "pending", "trials": []} for b in backends}
    for backend in backends:
        (folder / backend).mkdir()
    for iteration in range(-args.warmups, args.repeats):
        shift = iteration % len(backends)
        for backend in backends[shift:]+backends[:shift]:
            record = records[backend]
            if record["status"] == "failed":
                continue
            destination = folder / backend
            command = [str(args.build_dir / "build_halfspace_tree"), str(cloud["path"]), str(destination),
                       "--pipeline", "alpha-tet" if backend == "alpha-tet" else "octree",
                       "--voxel-size", str(size), "--validation-samples", str(args.validation_samples)]
            if backend == "alpha-tet":
                command += ["--alpha-voxel-factor", "3", "--offset-voxel-factor", "1", "--plane-reduction", "exact"]
            else:
                mode = "rect" if backend == "octree-boundary" else backend.removeprefix("octree-")
                command += ["--octree-pruning", mode, "--octree-boundary", "exact" if backend == "octree-boundary" else "none"]
            try:
                wall = execute(command, destination / f"trial_{iteration}.log", args.timeout)
                metrics = json.loads((destination / "benchmark.json").read_text())[0]
                report = json.loads((destination / "validation_report.json").read_text())
                if report["failures"] or metrics["raw_points_outside_tree"] or metrics["false_negative"]:
                    raise RuntimeError("Backend containment/validation failed")
                if backend != "alpha-tet" and (not report["passed"] or metrics["false_positive"]):
                    raise RuntimeError("Octree exact union validation failed")
                metrics["process_wall_ms"] = wall
                fingerprint = {f: digest(destination / f) for f in
                               ("convex_clusters.json", "logic_tree_raw.json", "logic_tree_simplified.json")}
                if "fingerprint" in record and fingerprint != record["fingerprint"]:
                    raise RuntimeError("Repeated run changed exported geometry/tree")
                record.update(status="passed", fingerprint=fingerprint, metrics=metrics)
                if iteration >= 0:
                    record["trials"].append(metrics)
            except (RuntimeError, ValueError, KeyError, OSError) as error:
                record.update(status="failed", error=str(error))
    for backend, record in records.items():
        write_json(folder / backend / "trials.json", record)
        if record["status"] == "passed":
            record["timings_ms"] = {key: stats([r[key] for r in record["trials"]]) for key in
                                    ("core_runtime_ms", "validation_ms", "total_ms", "process_wall_ms")}
    return records


def voxel_contains(point, occupied, size):
    base = tuple(math.floor(x/size) for x in point)
    # Mirror rounded closed grid faces; floor alone is wrong at decimal boundaries.
    for shift in itertools.product((0, -1, 1), repeat=3):
        key = tuple(k+d for k, d in zip(base, shift))
        if key in occupied and all(k*size <= p <= (k+1)*size for k, p in zip(key, point)):
            return True
    return False


def query_set(args, cloud, size, folder):
    octree = json.loads((folder / "octree-rect/octree.json").read_text())
    keys = [tuple(k) for k in octree["occupied_voxel_indices"]]
    low = [min(k[j] for k in keys)*size for j in range(3)]
    high = [(max(k[j] for k in keys)+1)*size for j in range(3)]
    if cloud["truth"]:
        truth_low, truth_high = bounds(cloud["truth"])
        low = [min(x, y) for x, y in zip(low, truth_low)]
        high = [max(x, y) for x, y in zip(high, truth_high)]
    # Extend to include both exported geometries. Do not assume requested offset
    # bounds the actual wrap or that containment recovery never increases it.
    wrap = folder / "alpha-tet/wrap.off"
    if wrap.exists():
        vertices, _ = read_off(wrap)
        low = [min(low[j], min(v[j] for v in vertices)) for j in range(3)]
        high = [max(high[j], max(v[j] for v in vertices)) for j in range(3)]
    low = [x-size for x in low]
    high = [x+size for x in high]
    rng = random.Random(args.seed+500)
    groups = {"uniform": [tuple(rng.uniform(a, b) for a, b in zip(low, high)) for _ in range(args.query_samples)],
              "raw": cloud["points"], "voxel_probes": []}
    for key in keys:
        lo, hi = [k*size for k in key], [(k+1)*size for k in key]
        groups["voxel_probes"].extend(tuple(hi[j] if bits[j] else lo[j] for j in range(3))
                                       for bits in itertools.product((0, 1), repeat=3))
        groups["voxel_probes"].append(tuple(a+(b-a)/2 for a, b in zip(lo, hi)))
    if cloud["truth"]:
        groups["ideal_surface"] = surface_points(cloud["truth"], args.truth_samples, args.seed+501)
        groups["ideal_interior"] = interior_points(cloud["truth"], args.truth_samples, args.seed+502)
        if "gap" in cloud["truth"]:
            a, b = cloud["truth"]["gap"]
            groups["ideal_gap"] = [tuple(rng.uniform(lo, hi) for lo, hi in zip(a, b)) for _ in range(args.truth_samples)]
    queries, ranges = [], {}
    for name, points in groups.items():
        ranges[name] = [len(queries), len(queries)+len(points)]
        queries.extend(points)
    path = folder / "queries.xyz"
    write_xyz(path, queries)
    metadata = {"ranges": ranges, "uniform_domain": [low, high], "domain_volume_m3": box_volume((low, high)),
                "uniform_seed": args.seed+500, "sha256": digest(path), "length_unit": "m"}
    write_json(folder / "queries.json", metadata)
    occupied = set(keys)
    metadata["uniform_voxel_inside"] = [voxel_contains(p, occupied, size) for p in groups["uniform"]]
    return path, metadata


def read_off(path):
    fields = " ".join(line.split("#", 1)[0] for line in path.read_text().splitlines()).split()
    if not fields or fields[0] != "OFF":
        raise ValueError("Expected OFF wrap")
    vertices_count, faces_count = int(fields[1]), int(fields[2])
    position = 4
    vertices = []
    for _ in range(vertices_count):
        vertices.append(tuple(map(float, fields[position:position+3])))
        position += 3
    triangles = []
    for _ in range(faces_count):
        if int(fields[position]) != 3:
            raise ValueError("Expected triangular wrap mesh")
        triangles.append(tuple(map(int, fields[position+1:position+4])))
        position += 4
    if position != len(fields) or not vertices or not triangles:
        raise ValueError("Malformed OFF wrap")
    return vertices, triangles


def mesh_volume(path):
    vertices, triangles = read_off(path)
    origin = vertices[0]
    translated = [tuple(x-o for x, o in zip(p, origin)) for p in vertices]
    terms = []
    for triangle in triangles:
        a, b, c = [translated[i] for i in triangle]
        terms.append(sum(a[j]*(b[(j+1)%3]*c[(j+2)%3]-b[(j+2)%3]*c[(j+1)%3]) for j in range(3))/6)
    volume = math.fsum(terms)
    if not math.isfinite(volume) or volume <= 0:
        raise ValueError("Nonpositive/nonfinite oriented wrap volume")
    return volume


def interval(count, total):
    # Wilson binomial interval: sampling uncertainty, NOT a geometric certificate.
    z = 1.959963984540054
    p = count/total
    center = (p+z*z/(2*total))/(1+z*z/total)
    radius = z*math.sqrt(p*(1-p)/total+z*z/(4*total*total))/(1+z*z/total)
    return [0. if count == 0 else max(0., center-radius), 1. if count == total else min(1., center+radius)]


def analyze(args, cloud, size, folder, records):
    if records["octree-rect"]["status"] != "passed":
        raise RuntimeError("Cannot construct common occupied-voxel reference: octree-rect failed")
    path, queries = query_set(args, cloud, size, folder)
    occupied_volume = records["octree-rect"]["metrics"]["occupied_volume_m3"]
    decisions = {}
    for backend, record in records.items():
        if record["status"] != "passed":
            continue
        destination = folder / backend
        try:
            command = [str(args.build_dir / "halfspace_tree_benchmark"), str(destination / "logic_tree_simplified.json"),
                       str(path), str(destination / "query_evaluation.json"), "--reference", str(destination / "logic_tree_raw.json"),
                       "--timed-count", str(min(256, args.query_samples)), "--rounds", str(args.query_rounds),
                       "--repeats", str(args.query_repeats)]
            execute(command, destination / "query.log", args.timeout)
            evaluation = json.loads((destination / "query_evaluation.json").read_text())
            values = evaluation["values"]
            inside = [x <= 0 for x in values]
            decisions[backend] = inside
            quality = {}
            for group, (start, end) in queries["ranges"].items():
                count = sum(inside[start:end])
                positive = [v for v in values[start:end] if v > 0]
                quality[group] = {"samples": end-start, "inside": count, "outside": end-start-count,
                                  "inside_fraction": count/(end-start),
                                  "max_positive_value_m": max(positive, default=0.),
                                  "positive_within_1e_10_m": sum(v <= 1e-10 for v in positive),
                                  "positive_above_1e_10_m": sum(v > 1e-10 for v in positive)}
            # Strict zero-level coverage is additionally checked here, even when
            # the backend's own numerical validator used a nonzero tolerance.
            if quality["raw"]["outside"]:
                raise RuntimeError(f"Strict exported tree misses {quality['raw']['outside']} raw points")
            if backend != "alpha-tet" and quality["voxel_probes"]["outside"]:
                raise RuntimeError("Exact Octree tree misses an occupied voxel probe")
            uniform = inside[:args.query_samples]
            reference = queries["uniform_voxel_inside"]
            added = sum(c and not v for c, v in zip(uniform, reference))
            missed = sum(v and not c for c, v in zip(uniform, reference))
            if backend != "alpha-tet" and (added or missed):
                raise RuntimeError("Exact Octree tree differs from occupied union on common queries")
            domain = queries["domain_volume_m3"]
            quality["common_uniform"] = {"samples": len(uniform), "added_samples_vs_voxels": added,
                                          "missed_samples_vs_voxels": missed,
                                          "added_volume_estimate_m3": domain*added/len(uniform),
                                          "missed_volume_estimate_m3": domain*missed/len(uniform),
                                          "added_volume_wilson95_m3": [domain*x for x in interval(added, len(uniform))],
                                          "missed_volume_wilson95_m3": [domain*x for x in interval(missed, len(uniform))]}
            volume = mesh_volume(destination / "wrap.off") if backend == "alpha-tet" else record["metrics"]["output_volume_m3"]
            quality.update(output_volume_m3=volume, occupied_volume_m3=occupied_volume,
                           output_volume_over_voxels_minus_one=volume/occupied_volume-1,
                           occupied_union_certificate=backend != "alpha-tet",
                           voxel_probe_coverage_is_sampled=True,
                           volume_method="oriented_wrap_triangles_float" if backend == "alpha-tet" else "disjoint_AABB_sum")
            if cloud["truth"]:
                truth_volume = solid_volume(cloud["truth"])
                quality.update(ideal_solid_volume_m3=truth_volume,
                               output_volume_over_ideal_solid_minus_one=volume/truth_volume-1)
            record["quality"] = quality
            record["query_timing"] = {k: v for k, v in evaluation.items() if k != "values"}
            write_json(destination / "quality.json", quality)
        except (RuntimeError, ValueError, KeyError, OSError) as error:
            record.update(status="failed", error=str(error))
    reference = decisions.get("octree-rect")
    for backend in records:
        if backend.startswith("octree-") and backend in decisions and reference is not None:
            mismatch = sum(a != b for a, b in zip(reference, decisions[backend]))
            records[backend]["cross_pruning_decision_mismatch"] = mismatch
            if mismatch:
                records[backend].update(status="failed", error="Cross-pruning tree decision mismatch")
    if "octree-boundary" in records and records["octree-boundary"]["status"] == "passed":
        if records["octree-boundary"]["fingerprint"] != records["octree-rect"]["fingerprint"]:
            records["octree-boundary"].update(status="failed", error="Boundary stage changed tree or partition")


def save_summary(args, cases):
    failures = [{"case": c["case"], "backend": b, "error": r.get("error", "incomplete")}
                for c in cases for b, r in c["backends"].items() if r["status"] != "passed"]
    summary = {"phase": 4, "length_unit": "m", "status": "failed" if failures else ("passed" if len(cases) == args.expected_cases else "in_progress"),
               "matrix_complete": len(cases) == args.expected_cases, "completed_cases": len(cases),
               "expected_cases": args.expected_cases, "measured_clouds": len(args.measured_input),
               "protocol": {"repeats": args.repeats, "warmups": args.warmups, "seed": args.seed,
                            "validation_samples": args.validation_samples, "query_samples": args.query_samples,
                            "truth_samples": args.truth_samples, "query_rounds": args.query_rounds,
                            "query_repeats": args.query_repeats, "sequential": True,
                            "alpha_voxel_factor": 3, "offset_voxel_factor": 1,
                            "noise_model": "bounded independent coordinate uniform jitter; no outliers",
                            "timing_note": "core excludes validation/I/O; process wall includes both; query adapter is separate",
                            "quality_note": "strict raw checks; corner/center and random coverage are not an Alpha voxel-union certificate",
                            "truth_note": "unobserved solid coverage and gap probes are diagnostics, not point-containment pass criteria"},
               "platform": {"system": platform.platform(), "machine": platform.machine()},
               "binary_sha256": args.binary_hashes, "failures": failures, "cases": cases}
    paired = [c for c in cases if all(c["backends"][b]["status"] == "passed" for b in BACKENDS)]
    if paired:
        speedups = [c["backends"]["alpha-tet"]["timings_ms"]["core_runtime_ms"]["median"] /
                    c["backends"]["octree-rect"]["timings_ms"]["core_runtime_ms"]["median"] for c in paired]
        summary["aggregate"] = {"paired_cases": len(paired), "octree_core_speedup": stats(speedups),
                                "octree_smaller_tree_cases": sum(c["backends"]["octree-rect"]["metrics"]["tree_nodes"] <
                                                                  c["backends"]["alpha-tet"]["metrics"]["tree_nodes"] for c in paired),
                                "octree_faster_query_cases": sum(c["backends"]["octree-rect"]["query_timing"]["median_us_per_query"] <
                                                                 c["backends"]["alpha-tet"]["query_timing"]["median_us_per_query"] for c in paired)}
    write_json(args.output / "summary.json", summary)
    rows = []
    for case in cases:
        for backend, record in case["backends"].items():
            row = {"case": case["case"], "scene": case["scene"], "provenance": case["provenance"],
                   "raw_points": case["raw_points"], "voxel_size_m": case["voxel_size_m"],
                   "noise_halfwidth_m": case["noise_halfwidth_m"], "backend": backend, "status": record["status"]}
            if record["status"] == "passed":
                metrics, quality = record["metrics"], record["quality"]
                row.update(occupied_voxels=metrics["occupied_voxels"],
                           convex_cells=metrics["convex_clusters"] if backend == "alpha-tet" else metrics["convex_cells"],
                           planes=metrics["planes_after_reduce"], tree_nodes=metrics["tree_nodes"],
                           raw_outside=quality["raw"]["outside"], voxel_probes_outside=quality["voxel_probes"]["outside"],
                           voxel_probes_max_positive_m=quality["voxel_probes"]["max_positive_value_m"],
                           voxel_probes_positive_above_1e_10_m=quality["voxel_probes"]["positive_above_1e_10_m"],
                           volume_over_voxels_minus_one=quality["output_volume_over_voxels_minus_one"],
                           query_median_us=record["query_timing"]["median_us_per_query"],
                           actual_alpha=metrics.get("alpha"), actual_offset=metrics.get("offset"),
                           offset_retries=metrics.get("offset_retries"))
                for key, value in record["timings_ms"].items():
                    row[key+"_median"] = value["median"]
                for group in ("ideal_surface", "ideal_interior", "ideal_gap"):
                    row[group+"_inside_fraction"] = quality.get(group, {}).get("inside_fraction")
            else:
                row["error"] = record.get("error")
            rows.append(row)
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with (args.output / "summary.csv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    lines = ["# Phase 4 benchmark", "", f"Completed cases: {len(cases)}/{args.expected_cases}; failures: {len(failures)}.",
             "Core ms excludes validation and I/O. Query us uses the existing C++ recursive evaluator.",
             "Only user-supplied measured inputs are real scans; built-in camera data is simulated.", "",
             "Voxel probes missed uses strict F>0, including tiny numerical residuals; inspect max positive values in CSV/quality JSON.", "",
             "| Case | Backend | Core ms | Query us | Cells | Tree nodes | Voxel probes missed | Solid interior coverage | Gap filled |",
             "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in rows:
        if row["status"] != "passed":
            lines.append(f"| {row['case']} | {row['backend']} FAILED | | | | | | | |")
            continue
        fraction = lambda key: "—" if row.get(key) is None else f"{100*row[key]:.1f}%"
        lines.append(f"| {row['case']} | {row['backend']} | {row['core_runtime_ms_median']:.3f} | {row['query_median_us']:.3f} | "
                     f"{row['convex_cells']} | {row['tree_nodes']} | {row['voxel_probes_outside']} | "
                     f"{fraction('ideal_interior_inside_fraction')} | {fraction('ideal_gap_inside_fraction')} |")
    (args.output / "report.md").write_text("\n".join(lines)+"\n")
    return failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="new directory; existing results will not be overwritten")
    parser.add_argument("--build-dir", type=Path, default=REPO / "build")
    parser.add_argument("--scenes", nargs="+", choices=list(SCENES), default=list(SCENES))
    parser.add_argument("--densities", type=int, nargs="+", default=[256, 1024])
    parser.add_argument("--noise-levels", type=float, nargs="+", default=[0., .003])
    parser.add_argument("--voxel-sizes", type=float, nargs="+", default=[.02, .04])
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--validation-samples", type=int, default=1000)
    parser.add_argument("--query-samples", type=int, default=1024)
    parser.add_argument("--truth-samples", type=int, default=512)
    parser.add_argument("--query-rounds", type=int, default=5)
    parser.add_argument("--query-repeats", type=int, default=7)
    parser.add_argument("--seed", type=int, default=20260908)
    parser.add_argument("--timeout", type=float, default=120.)
    parser.add_argument("--skip-camera", action="store_true")
    parser.add_argument("--measured-input", type=Path, action="append", default=[], help="optional actual XYZ scan, already in meters")
    parser.add_argument("--include-baselines", action="store_true", help="also run octree none/exact/boundary")
    args = parser.parse_args()
    if (min(args.densities) < 4 or min(args.repeats, args.query_samples, args.truth_samples, args.query_rounds, args.query_repeats) < 1
            or min(args.warmups, args.validation_samples) < 0 or not math.isfinite(args.timeout) or args.timeout <= 0
            or any(not math.isfinite(x) or x <= 0 for x in args.voxel_sizes)
            or any(not math.isfinite(x) or x < 0 for x in args.noise_levels)):
        parser.error("Invalid sample, repeat, noise, voxel or timeout value")
    for values in (args.scenes, args.densities, args.noise_levels, args.voxel_sizes):
        if len(set(values)) != len(values):
            parser.error("Repeated matrix values are not allowed")
    for values in (args.noise_levels, args.voxel_sizes):
        if len({f"{v:g}" for v in values}) != len(values):
            parser.error("Matrix values collide in output directory labels; use distinguishable values")
    args.output = args.output.resolve()
    args.build_dir = args.build_dir.resolve()
    if args.output.exists():
        parser.error("Output directory already exists; use a new directory")
    args.binary_hashes = {}
    for binary in ("build_halfspace_tree", "halfspace_tree_benchmark"):
        path = args.build_dir / binary
        if not path.is_file():
            parser.error(f"Build the required target first: {path}")
        args.binary_hashes[binary] = digest(path)
    for path in args.measured_input:
        read_xyz(path)  # reject invalid external input before creating outputs
    args.output.mkdir(parents=True)
    clouds = make_clouds(args)
    args.expected_cases = len(clouds)*len(args.voxel_sizes)
    backends = list(BACKENDS) + (["octree-none", "octree-exact", "octree-boundary"] if args.include_baselines else [])
    cases = []
    for cloud, size in itertools.product(clouds, args.voxel_sizes):
        name = f"{cloud['name']}_h{size:g}"
        print(f"[{len(cases)+1}/{args.expected_cases}] {name}", flush=True)
        folder = args.output / name
        folder.mkdir()
        records = build_case(args, cloud, size, folder, backends)
        try:
            analyze(args, cloud, size, folder, records)
        except (RuntimeError, ValueError, KeyError, OSError) as error:
            for record in records.values():
                if record["status"] == "passed" and "quality" not in record:
                    record.update(status="failed", error="Quality analysis failed: " + str(error))
        cases.append({"case": name, "scene": cloud["scene"], "provenance": cloud["provenance"],
                      "raw_points": len(cloud["points"]), "voxel_size_m": size,
                      "noise_halfwidth_m": cloud["noise_halfwidth_m"], "backends": records})
        failures = save_summary(args, cases)
        for backend, record in records.items():
            print(f"  {backend}: " + (f"core={record['timings_ms']['core_runtime_ms']['median']:.3f} ms, "
                  f"tree={record['metrics']['tree_nodes']}" if record["status"] == "passed" else record["error"]), flush=True)
    print(f"Summary: {args.output / 'summary.json'}; failures={len(failures)}", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
