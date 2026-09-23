#!/usr/bin/env python3
"""Five-way volume/X/Y/Z/best-axis benchmark, including all candidate costs."""
import argparse
import hashlib
import itertools
import json
import random
from pathlib import Path
import statistics
import subprocess
import time

from benchmark_scene_clouds import SCENES, surface_points
from benchmark_octree_phase4 import read_xyz, write_xyz

REPO = Path(__file__).resolve().parents[1]
MODES = ("volume","sweep-x","sweep-y","sweep-z","best-axis")


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def execute(command, path):
    save(path.with_suffix(".command.json"),command)
    start = time.perf_counter()
    p = subprocess.run(command,capture_output=True,text=True,timeout=120)
    wall = 1000*(time.perf_counter()-start)
    path.write_text(p.stdout+p.stderr)
    if p.returncode:
        raise RuntimeError(f"Failed benchmark: {path}")
    return wall


def fixtures():
    grid = list(itertools.product(range(3),repeat=3))
    rng = random.Random(731)
    indices = {
        "single":[(0,0,0)],
        "solid_cuboid":grid,
        "adjacent_cuboids":list(itertools.product(range(4),range(2),range(2))),
        "l_voxels":[(x,y,0) for x,y in itertools.product(range(3),repeat=2) if x==0 or y==0],
        "u_voxels":[(x,y,0) for x,y in itertools.product(range(3),repeat=2) if x!=1 or y==0],
        "staircase":[(x,y,0) for x,y in itertools.product(range(3),repeat=2) if y<=x],
        "closed_cavity":[p for p in grid if p!=(1,1,1)],
        "narrow_corridor":[p for p in grid if p[0]!=1],
        "two_disconnected":[(0,0,0),(3,2,1)],
        "thin_wall":list(itertools.product(range(1),range(4),range(4))),
        "random_voxel_clusters":[p for p in itertools.product(range(4),repeat=3) if rng.randrange(3)==0],
    }
    cases = [(name,[tuple(.04*(v+.5) for v in p) for p in cloud]) for name,cloud in indices.items()]
    return cases+[("u_surface256",surface_points(SCENES["u"],256,20260908)),
                  ("ellipsoid_surface256",surface_points(SCENES["ellipsoid"],256,20260908)),
                  ("camera",read_xyz(REPO/"data/realistic_sparse_camera.xyz")),
                  ("sparse_disconnected",[(.4*i+.02,.52*i+.02,.68*i+.02) for i in range(30)])]


def command(binary,cloud,out,mode,samples=1000,query=False):
    cmd = [str(binary),str(cloud),str(out),"--pipeline","octree","--voxel-size",".04",
           "--tree-source","volume" if mode=="volume" else "boundary","--validation-samples",str(samples)]
    if mode!="volume":
        cmd += ["--boundary-sweep",mode.removeprefix("sweep-")]
    if query:
        cmd += ["--query-benchmark","true"]
    return cmd


def check_selection(folder):
    best = json.loads((folder/"best-axis"/"axis_candidates.json").read_text())
    candidates = best["candidates"]
    assert [c["axis"] for c in candidates] == ["x","y","z"]
    winner = min(candidates,key=lambda c:(c["fallback"],c["tree_nodes"],c["plane_references"],c["finite_prisms"],c["axis"]))
    assert winner["axis"] == best["selected_axis"]
    for c in candidates:
        single = json.loads((folder/("sweep-"+c["axis"])/"axis_candidates.json").read_text())["candidates"][0]
        for key in ("fallback","fallback_reason","tree_nodes","plane_references","finite_prisms","unique_halfspaces"):
            assert c[key] == single[key], (folder,c["axis"],key)
    filename = "obstacle_tree_simplified.json" if winner["fallback"] else "free_closure_tree_simplified.json"
    assert (folder/"best-axis"/filename).read_bytes() == (folder/("sweep-"+winner["axis"])/filename).read_bytes()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--build-dir",type=Path,default=REPO/"build")
    parser.add_argument("--repeats",type=int,default=5)
    parser.add_argument("--warmups",type=int,default=1)
    parser.add_argument("--cases",nargs="+",help="Optional subset of fixture names; default all 15")
    args = parser.parse_args()
    if args.repeats<1 or args.warmups<0:
        parser.error("Invalid repeat/warmup count")
    cases = fixtures()
    if args.cases:
        if set(args.cases)-{name for name,_ in cases}:
            parser.error("Unknown fixture name")
        cases = [(name,p) for name,p in cases if name in args.cases]
    args.output = args.output.resolve()
    args.output.mkdir(parents=True,exist_ok=False)
    binary = args.build_dir.resolve()/"build_halfspace_tree"
    manifest = {"binary_sha256":hashlib.sha256(binary.read_bytes()).hexdigest(),"repeats":args.repeats,
                "warmups":args.warmups,"voxel_size":.04,"modes":MODES,"inputs":[],
                "ordering":"serial rotated five-way", "best_axis_cost":"all three builds plus selection; shared preprocessing once",
                "baseline_query":"median of four paired baseline query estimates"}
    rows = []
    for name,points in cases:
        print(name,flush=True)
        folder = args.output/name
        folder.mkdir()
        cloud = folder/"input.xyz"
        write_xyz(cloud,points)
        manifest["inputs"].append({"name":name,"raw_points":len(points),"sha256":hashlib.sha256(cloud.read_bytes()).hexdigest()})
        save(args.output/"manifest.json",manifest)
        records = {mode:[] for mode in MODES}
        for mode in records:
            (folder/mode).mkdir()
        for trial in range(-args.warmups,args.repeats):
            rotation = trial%len(MODES)
            for mode in MODES[rotation:]+MODES[:rotation]:
                out = folder/mode
                wall = execute(command(binary,cloud,out,mode),out/f"trial_{trial}.log")
                metrics = json.loads((out/"benchmark.json").read_text())[0]
                report = json.loads((out/"validation_report.json").read_text())
                if not report["passed"]:
                    raise RuntimeError(f"Validation failed: {out}")
                if trial>=0:
                    record = {"metrics":metrics,"process_wall_ms":wall}
                    if mode!="volume":
                        record["axis_candidates"] = json.loads((out/"axis_candidates.json").read_text())
                        assert metrics["core_runtime_ms"] >= metrics["all_candidates_ms"]+metrics["selection_ms"]
                    records[mode].append(record)
        check_selection(folder)
        # Each boundary configuration gets a separate process paired with the
        # baseline and the SAME seeded 1024 queries. Zero guards are timed.
        queries = {}
        for mode in MODES[1:]:
            query_folder = folder/"query_comparison"/mode
            query_folder.mkdir(parents=True)
            execute(command(binary,cloud,query_folder,mode,0,True),query_folder/"query.log")
            queries[mode] = json.loads((query_folder/"benchmark.json").read_text())[0]
        baseline_query = statistics.median(q["baseline_strict_free_query_us"] for q in queries.values())
        for mode,runs in records.items():
            save(folder/mode/"trials.json",runs)
            m = runs[-1]["metrics"]
            query = queries.get(mode)
            median_metric = lambda key: statistics.median(r["metrics"][key] for r in runs) if key in m else None
            row = {"case":name,"mode":mode,"raw_points":len(points),"occupied_voxels":m["occupied_voxels"],
                   "core_ms":statistics.median(r["metrics"]["core_runtime_ms"] for r in runs),
                   "build_validation_ms":statistics.median(r["metrics"]["total_ms"] for r in runs),
                   "wall_ms":statistics.median(r["process_wall_ms"] for r in runs),
                   "tree_nodes":m["tree_nodes"],"regions":m.get("free_prisms_after_merge",m.get("convex_cells")),
                   "plane_refs":m.get("plane_references",m.get("planes_after_reduce")),
                   "unique_planes":m.get("unique_halfspaces",m.get("unique_planes")),
                   "strict_query_us":query["strict_free_query_us"] if query else baseline_query,
                   "paired_baseline_query_us":query["baseline_strict_free_query_us"] if query else baseline_query,
                   "zero_guard_fraction":query["query_zero_guard_fraction"] if query else 0,
                   "boundary_strict_query_us":query["boundary_strict_free_query_us"] if query else None,
                   "selected_axis":m.get("selected_axis","-"),"candidate_count":m.get("candidate_count",0),
                   "all_candidates_ms":median_metric("all_candidates_ms"),"selection_ms":median_metric("selection_ms"),
                   "selected_candidate_build_ms":median_metric("selected_candidate_build_ms"),
                   "fallback":bool(m.get("fallback",False)),"fallback_reason":m.get("fallback_reason",""),
                   "certificate_ms":median_metric("certificate_ms"),"passed":True}
            if row["fallback"]:
                row["regions"] = m["baseline_regions"]
            rows.append(row)
    save(args.output/"manifest.json",manifest)
    save(args.output/"summary.json",rows)
    report = ["# Boundary XYZ and best-axis benchmark","",f"Sequential rotated five-way order; {args.warmups} warmups, {args.repeats} measured trials; h=0.04 m.","",
              "The boundary tree is a CLOSED nonoccupied-domain closure. Strict queries include the zero-value occupancy guard.",
              "Fallback uses the original obstacle rect tree and root>0. Nonoccupied does not mean observed safe free space.","",
              "best-axis minimizes (fallback, nodes, plane refs, finite prisms, axis X/Y/Z), not build time. All three candidate builds and selection are charged.","",
              "| Case | Mode | Selected | Regions | Plane refs | Nodes | Core ms | Build+validation ms | Strict query μs | Fallback |",
              "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |"]
    for r in rows:
        report.append(f'| {r["case"]} | {r["mode"]} | {r["selected_axis"]} | {r["regions"]} | {r["plane_refs"]} | {r["tree_nodes"]} | '
                      f'{r["core_ms"]:.3f} | {r["build_validation_ms"]:.3f} | {r["strict_query_us"]:.3f} | {r["fallback_reason"] or "no"} |')
    report += ["","Core excludes full validation, I/O and query benchmarking. Wall times and individual trials are retained.",
               "Query medians use seven five-round batches over 1024 identical uniform queries; geometry build medians use the configured outer repeats.",
               "The volume query is the median of four paired baseline estimates. Each pair is retained; identical trees can have timing noise.",
               "The query_comparison subfolders also report raw/simplified scalar query time and boundary-point guard stress time.",
               "Phase timings sum all attempted candidates; axis_candidates.json and trials.json retain individual costs and selection certificates.",
               "Finite box count excludes the six unbounded exterior branches. Fallback region count instead refers to occupied rect boxes.",
               "All three best-axis candidates are validated. No hard deadline, CBF gradient or QP is implemented."]
    (args.output/"report.md").write_text("\n".join(report)+"\n")
    print(args.output/"report.md",flush=True)


if __name__=="__main__":
    main()
