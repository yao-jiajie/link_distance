#!/usr/bin/env python3
"""Bounded hull vs exact rect; does not require old/new scalar equality."""
import argparse
import hashlib
import json
from pathlib import Path
import random
import statistics
import subprocess
import time

from benchmark_scene_clouds import SCENES, surface_points
from benchmark_octree_phase4 import read_xyz, write_xyz

REPO = Path(__file__).resolve().parents[1]


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False)+"\n")


def execute(command, log):
    write_json(log.with_suffix(".command.json"), command)
    start = time.perf_counter()
    p = subprocess.run(command, capture_output=True, text=True, timeout=120)
    wall = 1000*(time.perf_counter()-start)
    log.write_text(p.stdout+p.stderr)
    if p.returncode:
        raise RuntimeError(f"Benchmark failed: {log}")
    return wall


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--build-dir", type=Path, default=REPO/"build")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmups", type=int, default=1)
    args = parser.parse_args()
    if args.repeats < 1 or args.warmups < 0:
        parser.error("repeats must be positive; warmups must be nonnegative")
    args.output = args.output.resolve()
    args.build_dir = args.build_dir.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    inputs = args.output/"inputs"
    inputs.mkdir()
    cases = [("l_voxels", read_xyz(REPO/"tests/data/convex_l.xyz"), .01, "synthetic_voxel_centers")]
    for name in ("u", "ellipsoid"):
        cases.append((name+"_surface256", surface_points(SCENES[name],256,20260908),.04,"synthetic_surface"))
    cases.append(("camera",read_xyz(REPO/"data/realistic_sparse_camera.xyz"),.04,"simulated_RGBD"))
    # Deliberately include a tighter budget with weak/no merge benefit. These
    # are experiment settings, not robot safety recommendations.
    modes = {"rect": [],
             "hull_tight": ["--octree-approx","hull","--max-volume-inflation",".2"],
             "hull_relaxed": ["--octree-approx","hull","--max-volume-inflation",".2"],
             "alpha-tet": ["--alpha-voxel-factor","3","--offset-voxel-factor","1","--plane-reduction","exact"]}
    manifest = {"length_unit":"m", "time_unit":"ms", "repeats":args.repeats, "warmups":args.warmups,
                "binary_sha256":hashlib.sha256((args.build_dir/"build_halfspace_tree").read_bytes()).hexdigest(),
                "hull_volume_budget":.2, "tight_distance_factor":.5, "relaxed_distance_factor":.8,
                "hull_max_voxels":64,"hull_max_boxes":512,"clouds":[]}
    rows = []
    for name, points, h, provenance in cases:
        print(name, flush=True)
        cloud = inputs/(name+".xyz")
        write_xyz(cloud,points)
        manifest["clouds"].append({"name":name,"raw_points":len(points),"voxel_size":h,"provenance":provenance,
                                   "sha256":hashlib.sha256(cloud.read_bytes()).hexdigest()})
        folder = args.output/name
        folder.mkdir()
        trials = {mode:[] for mode in modes}
        for trial in range(-args.warmups,args.repeats):
            order = list(modes)
            order = order[trial%len(order):]+order[:trial%len(order)]
            for mode in order:
                out = folder/mode
                out.mkdir(exist_ok=True)
                command = [str(args.build_dir/"build_halfspace_tree"),str(cloud),str(out),
                           "--pipeline","alpha-tet" if mode == "alpha-tet" else "octree",
                           "--voxel-size",str(h),"--validation-samples","1000"]+modes[mode]
                if mode.startswith("hull"):
                    command += ["--max-fill-distance",str(h*(.5 if mode == "hull_tight" else .8)),
                                "--convex-max-voxels","64","--convex-test-max-boxes","512"]
                wall = execute(command,out/f"trial_{trial}.log")
                metrics = json.loads((out/"benchmark.json").read_text())[0]
                validation = json.loads((out/"validation_report.json").read_text())
                if validation.get("failures") or (mode != "alpha-tet" and not validation["passed"]):
                    raise RuntimeError(f"Failed validation: {out}")
                if trial >= 0:
                    trials[mode].append({"metrics":metrics,"process_wall_ms":wall})
        # Identical raw, occupied corners, and uniform probes for all trees.
        # The reference geometry is exact rect, NOT an unobserved solid.
        reference = json.loads((folder/"rect"/"octree.json").read_text())
        corners = sorted({tuple((key[k]+((mask>>k)&1))*h for k in range(3))
                          for key in reference["occupied_voxel_indices"] for mask in range(8)})
        rng = random.Random(829)
        bounds = [(min(p[k] for p in corners)-h,max(p[k] for p in corners)+h) for k in range(3)]
        probes = points+corners+[tuple(rng.uniform(lo,hi) for lo,hi in bounds) for _ in range(1000)]
        query_path = folder/"queries.xyz"
        write_xyz(query_path,probes)
        queries = {}
        for mode in modes:
            out = folder/mode
            command = [str(args.build_dir/"halfspace_tree_benchmark"),str(out/"logic_tree_simplified.json"),
                       str(query_path),str(out/"queries.json"),"--reference",str(out/"logic_tree_raw.json")]
            execute(command,out/"queries.log")
            queries[mode] = json.loads((out/"queries.json").read_text())
        for mode, runs in trials.items():
            write_json(folder/mode/"trials.json",runs)
            m = runs[-1]["metrics"]
            values = queries[mode]["values"]
            row = {"case":name,"mode":mode,"raw_points":len(points),"voxel_size":h,
                   "core_ms":statistics.median(r["metrics"]["core_runtime_ms"] for r in runs),
                   "build_validation_ms":statistics.median(r["metrics"]["total_ms"] for r in runs),
                   "wall_ms":statistics.median(r["process_wall_ms"] for r in runs),
                   "convex_cells":m.get("convex_cells",m.get("convex_clusters")),
                   "planes":m.get("planes_after_reduce",m.get("planes_after")),"tree_nodes":m["tree_nodes"],
                   "non_aabb_cells":0 if mode == "rect" else m.get("non_aabb_cells"),"volume_inflation":m.get("volume_inflation"),
                   "fill_distance_bound_m":m.get("maximum_fill_distance_bound_m"),
                   "query_us":queries[mode]["median_us_per_query"],
                   "raw_outside":sum(v>0 for v in values[:len(points)]),
                   "source_corner_outside":sum(v>0 for v in values[len(points):len(points)+len(corners)]),
                   "old_inside_new_outside":sum(a<=0 and b>0 for a,b in zip(queries["rect"]["values"],values)),
                   "added_occupied_samples":sum(a>0 and b<=0 for a,b in zip(queries["rect"]["values"],values)),
                   "passed":True}
            if row["raw_outside"] or (mode.startswith("hull") and (row["source_corner_outside"] or row["old_inside_new_outside"])):
                raise RuntimeError(f"Exported source coverage failed: {row}")
            rows.append(row)
    write_json(args.output/"manifest.json",manifest)
    write_json(args.output/"summary.json",rows)
    lines = ["# Bounded non-AABB hull benchmark", "",
             f"Sequential rotating order, {args.warmups} warmups / {args.repeats} trials. Units: m, ms, μs.", "",
             "Reference = occupied voxel union. Alpha-tet has no comparable volume/distance certificate here.", "",
             "| Case | Mode | Cells | Non-AABB | Plane refs | Tree | Core ms | Build+validation ms | Query μs | Inflation |",
             "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for r in rows:
        inflation = "—" if r["volume_inflation"] is None else f'{100*r["volume_inflation"]:.2f}%'
        non_aabb = "—" if r["non_aabb_cells"] is None else str(r["non_aabb_cells"])
        lines.append(f'| {r["case"]} | {r["mode"]} | {r["convex_cells"]} | {non_aabb} | {r["planes"]} | '
                     f'{r["tree_nodes"]} | {r["core_ms"]:.3f} | {r["build_validation_ms"]:.3f} | {r["query_us"]:.3f} | {inflation} |')
    lines += ["", "Core excludes complete validation and I/O; wall times and per-trial metrics are in summary.json/trials.json.",
              "Raw/simplified equality is checked within each backend. Added occupancy across backends is allowed, not hidden as an equality success.",
              "These are synthetic/simulated cases, not measured scans or hard realtime evidence."]
    (args.output/"report.md").write_text("\n".join(lines)+"\n")
    print(args.output/"report.md",flush=True)


if __name__ == "__main__":
    main()
