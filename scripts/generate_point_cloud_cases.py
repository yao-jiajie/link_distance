#!/usr/bin/env python3
"""Deterministic surface XYZ inputs for point-cloud-to-function regression tests."""
import argparse
import hashlib
import html
import json
import math
from pathlib import Path
import random
import shutil

from benchmark_octree_phase4 import read_xyz, write_xyz
from benchmark_scene_clouds import SCENES, surface_points

REPO = Path(__file__).resolve().parents[1]
SEED = 20260910


def make_cases():
    cases = {}

    def add(name, description, points, sampling, **metadata):
        cases[name] = {"name": name, "description": description, "points": points,
            "sampling": sampling, "suite_seed": SEED, "expected_free_probes_m": [], **metadata}

    camera = REPO/"data/realistic_sparse_camera.xyz"
    add("camera", "Original noisy multi-view ellipsoid", read_xyz(camera), "simulated_depth_camera",
        source_file="data/realistic_sparse_camera.xyz", seed=20260904)
    for offset, name in enumerate(("cube", "l", "u", "thin_wall")):
        add(name+"_surface", name.replace("_", " ").upper()+" surface", surface_points(SCENES[name],2048,SEED+offset),
            "box_union_surface_area", shape=SCENES[name])

    boxes = {
        "staircase_surface": {"boxes": [((-.18,-.14,.50),(-.06,.14,.62)),
            ((-.06,-.14,.50),(.06,.14,.74)), ((.06,-.14,.50),(.18,.14,.86))]},
        "two_objects_surface": {"boxes": [((-.30,-.10,.52),(-.10,.10,.72)),
            ((.10,-.10,.68),(.30,.10,.88))]},
        "narrow_corridor_surface": {"boxes": [((-.22,-.14,.52),(-.06,.14,.88)),
            ((.06,-.14,.52),(.22,.14,.88))]},
    }
    for offset, (name, shape) in enumerate(boxes.items()):
        add(name, name.removesuffix("_surface").replace("_", " "), surface_points(shape,2048,SEED+10+offset),
            "box_union_surface_area", shape=shape)
    cases["l_surface"]["expected_free_probes_m"] = [[.06,.06,.58]]
    cases["u_surface"]["expected_free_probes_m"] = [[.02,.06,.58]]
    for name in ("two_objects_surface", "narrow_corridor_surface"):
        cases[name]["expected_free_probes_m"] = [[.02,.02,.70]]

    rng = random.Random(SEED+20)
    radius, height, center = .14, .32, (0.,0.,.72)
    cylinder = []
    for _ in range(2048):
        angle = rng.uniform(0,2*math.pi)
        if rng.random() < height/(height+radius):
            radial = radius; z = rng.uniform(-height/2,height/2)
        else:
            radial = radius*math.sqrt(rng.random()); z = rng.choice((-height/2,height/2))
        cylinder.append((radial*math.cos(angle),radial*math.sin(angle),center[2]+z))
    add("cylinder_surface", "Closed cylinder surface", cylinder, "cylinder_surface_area",
        shape={"type": "cylinder", "radius_m": radius, "height_m": height, "center_m": center})

    rng = random.Random(SEED+21)
    major, minor = .16, .06
    torus = []
    for _ in range(2048):
        theta, phi = rng.uniform(0,2*math.pi), rng.uniform(0,2*math.pi)
        radial = major+minor*math.cos(phi)
        torus.append((radial*math.cos(theta),radial*math.sin(theta),.74+minor*math.sin(phi)))
    add("torus_surface", "Torus with an open central hole", torus, "uniform_torus_parameters_not_surface_area",
        shape={"type": "torus", "major_radius_m": major, "minor_radius_m": minor, "center_m": [0,0,.74]},
        expected_free_probes_m=[[.02,.02,.74]])

    outer = {"boxes": [((-.18,-.18,.56),(.18,.18,.92))]}
    inner = {"boxes": [((-.09,-.09,.65),(.09,.09,.83))]}
    hollow = surface_points(outer,1536,SEED+22)+surface_points(inner,512,SEED+23)
    add("hollow_box_surfaces", "Outer and inner surfaces of a hollow box", hollow, "outer_and_inner_box_surfaces",
        shape={"outer": outer, "inner": inner}, expected_free_probes_m=[[.02,.02,.74]])
    rng = random.Random(SEED+24)
    add("plane_patch", "Open planar patch", [(rng.uniform(-.18,.18),rng.uniform(-.14,.14),.74) for _ in range(1024)],
        "uniform_planar_rectangle", expected_free_probes_m=[[.02,.02,.66]])

    ellipsoid = surface_points(SCENES["ellipsoid"],8192,SEED+30)
    clean = ellipsoid[:2048]
    for name, count in (("ellipsoid_sparse",256), ("ellipsoid_surface",2048), ("ellipsoid_dense",8192)):
        add(name, f"Ellipsoid surface, {count} points", ellipsoid[:count], "uniform_directions_scaled_to_ellipsoid",
            shape=SCENES["ellipsoid"], nested_density_family="ellipsoid")
    for millimeters in (2,8):
        rng = random.Random(SEED+31)
        points = [tuple(v+rng.gauss(0,millimeters/1000) for v in p) for p in clean]
        add(f"ellipsoid_noise_{millimeters}mm", f"Ellipsoid with {millimeters} mm coordinate noise", points,
            "ellipsoid_plus_iid_cartesian_gaussian", coordinate_noise_sigma_m=millimeters/1000,
            observed_displacement_rms_m=math.sqrt(sum(sum((a-b)**2 for a,b in zip(p,q)) for p,q in zip(points,clean))/len(clean)),
            base_case="ellipsoid_surface")
    partial = [p for p in ellipsoid if p[0]>=0][:1024]
    add("ellipsoid_partial", "Open positive-X ellipsoid hemisphere", partial, "cropped_ellipsoid_surface",
        crop="x >= 0", expected_free_probes_m=[[-.10,.02,.74]])
    duplicated = clean*3
    random.Random(SEED+32).shuffle(duplicated)
    add("ellipsoid_duplicates", "Ellipsoid shuffled with three copies of every point", duplicated,
        "duplicate_and_shuffle", same_geometry_as="ellipsoid_surface")
    outliers = [(.26,.02,.74),(-.26,.02,.74),(.02,.22,.74),(.02,-.22,.74),(.02,.02,1.06),(.02,.02,.42)]
    add("ellipsoid_outliers", "Ellipsoid with six isolated outliers", clean+outliers,
        "ellipsoid_plus_explicit_outliers", outliers_m=outliers, base_case="ellipsoid_surface")
    return cases


def preview(path, cases):
    """Small dependency-free orthographic scatter gallery of the actual points."""
    width, panel_w, panel_h = 1200, 300, 250
    height = 45+panel_h*math.ceil(len(cases)/4)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        '<text x="16" y="27" font-family="sans-serif" font-size="18">XYZ input gallery — orthographic surface samples (up to 768 displayed per case)</text>']
    for index, case in enumerate(cases):
        x0, y0 = (index%4)*panel_w, 45+(index//4)*panel_h
        parts.append(f'<g transform="translate({x0},{y0})"><rect x="5" y="5" width="290" height="240" rx="7" fill="white" stroke="#cbd5e1"/>')
        parts.append(f'<text x="15" y="26" font-family="sans-serif" font-size="14">{html.escape(case["name"])}</text>')
        parts.append(f'<text x="15" y="44" font-family="sans-serif" font-size="12" fill="#64748b">{len(case["points"])} points · meters · same display scale</text>')
        points = case["points"][::max(1,math.ceil(len(case["points"])/768))]
        for x,y,z in points:
            u = 150+280*(x-y)/math.sqrt(2)
            v = 150-280*(.5*(x+y)/math.sqrt(2)+math.sqrt(.75)*(z-.74))
            parts.append(f'<circle cx="{u:.2f}" cy="{v:.2f}" r=".9" fill="#2563eb" opacity=".65"/>')
        parts.append('</g>')
    parts.append('</svg>')
    path.write_text("\n".join(parts)+"\n")


def write_cases(output, names=None):
    cases = make_cases()
    names = list(cases) if names is None else list(names)
    if not names or len(set(names))!=len(names) or set(names)-cases.keys():
        raise ValueError("Case names must be nonempty, known and unique")
    output.mkdir(parents=True,exist_ok=False)
    entries = []
    for name in names:
        case = cases[name]; path = output/(name+".xyz")
        if "source_file" in case:
            shutil.copyfile(REPO/case["source_file"],path)
        else:
            write_xyz(path,case["points"])
        entries.append({k: v for k,v in case.items() if k!="points"} | {
            "file": path.name, "raw_points": len(case["points"]),
            "bounds_m": [[min(p[a] for p in case["points"]) for a in range(3)],
                         [max(p[a] for p in case["points"]) for a in range(3)]],
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    manifest = {"schema": 1, "length_unit": "m", "seed": SEED, "recommended_voxel_size_m": .04,
        "semantics": "Surface point samples only; voxel occupancy does not fill object interiors. Synthetic cases other than camera have no sensor or visibility model.",
        "generator_sources_sha256": {name: hashlib.sha256((REPO/name).read_bytes()).hexdigest() for name in (
            "scripts/generate_point_cloud_cases.py", "scripts/benchmark_scene_clouds.py", "scripts/benchmark_octree_phase4.py")},
        "cases": entries}
    (output/"manifest.json").write_text(json.dumps(manifest,indent=2,allow_nan=False)+"\n")
    preview(output/"preview.svg",[cases[name] for name in names])
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--cases",nargs="+")
    args = parser.parse_args()
    write_cases(args.output.resolve(),args.cases)
    print(args.output.resolve()/"manifest.json")


if __name__ == "__main__":
    main()
