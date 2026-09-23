#!/usr/bin/env python3
"""Diverse XYZ-to-hard/LSE-function timings, with separately executed validation."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

from benchmark_boundary_tree import REPO, execute, save

GEOMETRY_FILES = ("free_direct_tree_raw.json", "free_direct_tree_simplified.json", "boundary.off",
    "boundary_patches.json", "boundary_supports.json", "direct_splits.json", "occupancy_reference.json")
STABLE_FILES = GEOMETRY_FILES+("boundary_sign_diagnostics.json", "smooth_distance.json", "distance_samples.csv")


def load(path):
    return json.loads(path.read_text())


def evaluate(tree, point, beta=None):
    """Independent exported-tree formula. Return positive-free value/gradient."""
    values, gradients = {}, {}
    for node in tree["nodes"]:
        if node["type"] == "leaf":
            leaf = tree["leaves"][node["leaf_id"]]
            value = sum(a*b for a,b in zip(leaf["normal"],point))-leaf["offset"]
            gradient = leaf["normal"]
        else:
            children = node["children"]
            args = [values[c] for c in children]
            maximum = node["type"] == "max"
            anchor = (max if maximum else min)(args)
            value, gradient = anchor, (0.,0.,0.)
            if beta is not None:
                sign = 1 if maximum else -1
                weights = [math.exp(sign*beta*(v-anchor)) for v in args]
                total = math.fsum(weights)
                value += sign*math.log(total if maximum else total/len(args))/beta
                gradient = [math.fsum(w/total*gradients[c][a] for w,c in zip(weights,children)) for a in range(3)]
        values[node["id"]] = value
        if beta is not None:
            gradients[node["id"]] = gradient
    return -values[tree["root_id"]], [-v for v in gradients[tree["root_id"]]] if beta is not None else None


def check_exports(out, case, h):
    tree = load(out/"free_direct_tree_raw.json")
    field = load(out/"smooth_distance.json")
    keys = load(out/"occupancy_reference.json")["occupied_keys"]
    patches = load(out/"boundary_patches.json")["patches"]
    boxes = [([min(p[a] for p in patch["corners_m"]) for a in range(3)],
              [max(p[a] for p in patch["corners_m"]) for a in range(3)]) for patch in patches]
    def occupied(point):
        return any(all(k*h <= v <= (k+1)*h for k,v in zip(key,point)) for key in keys)
    with (out/"distance_samples.csv").open() as stream:
        samples = list(csv.DictReader(stream))
    assert samples
    selected = sorted({round(i*(len(samples)-1)/23) for i in range(24)})
    for index in selected:
        row = samples[index]; point = [float(row[a+"_m"]) for a in "xyz"]
        hard,_ = evaluate(tree,point)
        smooth,gradient = evaluate(tree,point,field["beta"])
        distance = min(math.sqrt(sum(max(lo[a]-point[a],0,point[a]-hi[a])**2 for a in range(3))) for lo,hi in boxes)
        reference = -distance if occupied(point) else distance
        assert abs(hard-float(row["hard_proxy_m"])) < 1e-12
        assert abs(smooth-float(row["smooth_proxy_m"])) < 1e-12
        assert abs(reference-float(row["reference_sdf_m"])) < 1e-12
        assert all(abs(g-float(row["gradient_"+a])) < 1e-10 for a,g in zip("xyz",gradient))
        assert -1e-12 <= hard-smooth <= field["scalar_error_bound_m"]+1e-12
        assert math.sqrt(sum(g*g for g in gradient)) <= 1+1e-10
    for point in case["expected_free_probes_m"]:
        assert not occupied(point), (case["name"], "gap/hole occupied", point)
        assert evaluate(tree,point)[0] > 0 and evaluate(tree,point,field["beta"])[0] > 0, (case["name"],point)
    return {"independent_sample_checks": len(selected), "gap_probe_checks": len(case["expected_free_probes_m"])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--input-dir",type=Path,default=REPO/"data/point_cloud_cases")
    parser.add_argument("--build-dir",type=Path,default=REPO/"build")
    parser.add_argument("--cases",nargs="+")
    parser.add_argument("--repeats",type=int,default=5)
    parser.add_argument("--warmups",type=int,default=1)
    parser.add_argument("--validation-samples",type=int,default=1000)
    args = parser.parse_args()
    if args.repeats<1 or args.warmups<0 or args.validation_samples<0:
        parser.error("Invalid trial/sample counts")
    inputs = args.input_dir.resolve(); source = load(inputs/"manifest.json")
    cases = {case["name"]: case for case in source["cases"]}
    names = list(cases) if args.cases is None else args.cases
    if not names or len(set(names))!=len(names) or set(names)-cases.keys():
        parser.error("Case names must be known and unique")
    if source["length_unit"]!="m" or source["recommended_voxel_size_m"]!=.04:
        parser.error("Expected meter fixtures designed for h=.04 m")
    output = args.output.resolve(); output.mkdir(parents=True,exist_ok=False)
    binary = args.build_dir.resolve()/"build_halfspace_tree"
    manifest = {"binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "input_manifest_sha256": hashlib.sha256((inputs/"manifest.json").read_bytes()).hexdigest(),
        "input_dir": str(inputs), "cases": names, "voxel_size_m": .04, "lse_error_budget_m": .001,
        "repeats": args.repeats, "warmups": args.warmups, "validation_samples": args.validation_samples,
        "timed_metric": "core_runtime_ms: in-memory points through voxelization, boundary/tree construction, shared program compilation, LSE and value/gradient workspace preparation",
        "excluded_from_function_timing": "input/output I/O, full certificates, numerical probes, distance audit and query evaluation",
        "validation": "default work caps and all strata retained; independent exported hard/LSE/gradient/distance checks plus explicit gap probes",
        "order": "serial processes; case order rotates each round; each measured process constructs a fresh program"}
    save(output/"manifest.json",manifest)
    records = {name: [] for name in names}; hashes = {}; metrics = {}
    for name in names:
        path = inputs/cases[name]["file"]
        assert hashlib.sha256(path.read_bytes()).hexdigest()==cases[name]["sha256"], path
        (output/name).mkdir()
    for trial in range(-args.warmups,args.repeats):
        offset = trial%len(names)
        for name in names[offset:]+names[:offset]:
            out = output/name; cloud = inputs/cases[name]["file"]
            command = [str(binary),str(cloud),str(out),"--pipeline","octree","--tree-source","boundary",
                "--boundary-expression","direct","--boundary-query","dag","--voxel-size",".04",
                "--distance-field","lse","--lse-error",".001","--execution-sharing","structural",
                "--lse-kernel","binary","--validation-cache","auto","--sign-propagation","packed",
                "--validation-samples",str(args.validation_samples)]
            execute(command,out/f"trial_{trial}.log")
            m = load(out/"benchmark.json")[0]; metrics[name] = m
            assert m["passed"] and m["strict_sign_certified"] and m["sign_certificate_applied"]
            assert m["raw_points"]==cases[name]["raw_points"] and m["smooth_reuses_compiled_program"]==1
            assert m["execution_structural_sharing"]==m["sign_packed_propagation"]==m["smooth_binary_kernel"]==1
            for key in ("raw_points_classified_free","world_sign_errors","true_boundary_sign_errors",
                        "distance_false_free_count","distance_bound_errors","distance_upper_order_errors"):
                assert m[key]==0, (name,key,m[key])
            assert m["distance_gradient_checks"]>0 and m["distance_zero_rays_tested"]>0
            assert m["distance_gradient_max_abs_error"]<1e-4 and m["lse_error_bound_m"]<=.001+1e-12
            current = {f: hashlib.sha256((out/f).read_bytes()).hexdigest() for f in STABLE_FILES}
            if name in hashes:
                assert hashes[name]==current, (name,"nondeterministic output")
            hashes[name] = current
            if trial>=0:
                records[name].append(m); save(out/"trials.json",records[name])
        print(f"Completed round {trial+1}/{args.repeats}",flush=True)
    rows = []
    for name in names:
        case = cases[name]; m = metrics[name]; trials = records[name]
        values = sorted(t["core_runtime_ms"] for t in trials)
        checks = check_exports(output/name,case,.04)
        peer = case.get("same_geometry_as")
        if peer in names:
            assert all(hashes[name][f]==hashes[peer][f] for f in GEOMETRY_FILES), (name,"duplicate/order changed geometry")
        row = {"case": name, "description": case["description"], "raw_points": case["raw_points"],
            "occupied_voxels": int(m["occupied_voxels"]), "execution_nodes": int(m["smooth_execution_nodes"]),
            "strata_count": int(m["strata_count"]), "function_ms": statistics.median(values),
            "function_min_ms": values[0], "function_p95_ms": values[math.ceil(.95*len(values))-1],
            "function_trials_ms": [t["core_runtime_ms"] for t in trials],
            "validation_ms": statistics.median(t["validation_ms"] for t in trials),
            "geometry_ms": statistics.median(t["geometry_core_ms"] for t in trials),
            "gradient_checks": int(m["distance_gradient_checks"]), "zero_rays_tested": int(m["distance_zero_rays_tested"]),
            "deterministic_outputs": True, "same_geometry_peer": peer if peer in names else None,
            "passed": True, **checks}
        rows.append(row)
    save(output/"summary.json",rows)
    report = ["# Diverse point-cloud inputs: time to a queryable function", "",
        f"h=.04 m, 1 mm LSE scalar error budget; {args.warmups} warmup + {args.repeats} measured rounds, rotating case order.",
        "**Function time excludes validation and all file I/O.** Every run still executes full certification and numerical auditing after that timer.",
        "All generated shapes are surface samples; the occupied union is not an interior-filled reconstruction.", "",
        "| Input | Points | Voxels | Execution nodes | Function median ms | Min / p95 ms | Validation ms | Passed |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |"]
    for r in rows:
        report.append(f'| {r["case"]} | {r["raw_points"]} | {r["occupied_voxels"]} | {r["execution_nodes"]} | '
            f'{r["function_ms"]:.6f} | {r["function_min_ms"]:.6f} / {r["function_p95_ms"]:.6f} | {r["validation_ms"]:.3f} | yes |')
    report += ["", "Function includes occupancy, geometry/tree construction, the shared compiled hard program, LSE parameters and query workspaces. It excludes the first query itself.",
        "Each case passed all-strata strict signs, raw-point containment, LSE bounds, gradient and zero-ray auditing, and independent exported-formula checks on up to 24 sampled rows.",
        "Explicit probes check selected concavities, holes and gaps; duplicate/shuffle inputs must export identical geometry when both related cases are selected.",
        "No timing threshold is used as a correctness criterion. Point count alone does not determine construction time: occupied geometry and expression size also matter."]
    (output/"report.md").write_text("\n".join(report)+"\n")
    print(output/"report.md",flush=True)


if __name__ == "__main__":
    main()
