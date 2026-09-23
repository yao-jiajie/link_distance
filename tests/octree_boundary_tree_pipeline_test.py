"""Read exported strict classifier bundle independently; no C++ oracle calls."""
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile

program = Path(sys.argv[1]).resolve()


def load(path):
    return json.loads(path.read_text())


def evaluate(tree, point):
    values = {}
    for node in tree["nodes"]:
        if node["type"] == "leaf":
            p = tree["leaves"][node["leaf_id"]]
            values[node["id"]] = sum(a*b for a,b in zip(p["normal"],point))-p["offset"]
        else:
            values[node["id"]] = (min if node["type"] == "min" else max)(values[i] for i in node["children"])
    return values[tree["root_id"]]


def verify(output, points):
    manifest = load(output/"strict_free_space.json")
    tree = load(output/manifest["tree"])
    raw_tree = load(output/manifest["raw_tree"])
    reference = load(output/manifest["occupancy_reference"])
    nodes = {n["id"]:n for n in reference["nodes"]}
    def occupied(point, node_id=0):
        n = nodes[node_id]
        if not all(lo<=x<=hi for x,lo,hi in zip(point,n["lower"],n["upper"])):
            return False
        return n["leaf"] or any(occupied(point,i) for i in n["children"] if i>=0)
    # Independent brute force fine-box oracle includes ALL touching leaves.
    fine_boxes = [(n["lower"],n["upper"]) for n in nodes.values() if n["leaf"]]
    def expected(point):
        return not any(all(lo<=p<=hi for p,lo,hi in zip(point,*b)) for b in fine_boxes)
    def strict(point):
        value = evaluate(tree,point)
        if manifest["classifier"] == "root_gt_zero":
            return value>0
        assert manifest["classifier"] == "root_lt_zero_or_root_eq_zero_and_not_occupied"
        return value<0 or (value==0 and not occupied(point))
    coords = []
    for axis in range(3):
        planes = sorted({b[s][axis] for b in fine_boxes for s in (0,1)})
        coords.append([planes[0]-1]+sorted(planes+[(a+b)/2 for a,b in zip(planes,planes[1:])])+[planes[-1]+1])
    seen_free_zero = seen_occupied_zero = False
    for point in itertools.product(*coords):
        value = evaluate(tree,point)
        assert value == evaluate(raw_tree,point)
        assert strict(point) == expected(point), (point,manifest)
        seen_free_zero |= value==0 and expected(point)
        seen_occupied_zero |= value==0 and not expected(point)
    for p in points:
        assert not strict(p)
    m = load(output/"benchmark.json")[0]
    v = load(output/"validation_report.json")
    assert v["passed"] and v["integer_coverage_certificate"]
    assert m["false_free"] == m["missed_free"] == m["raw_points_classified_free"] == m["scalar_mismatches"] == 0
    assert m["tree_nodes"] == len(tree["nodes"])
    candidates = load(output/"axis_candidates.json")
    choices = candidates["candidates"]
    assert len(choices) == (3 if m["sweep_axis"]=="best-axis" else 1)
    assert len(choices) == m["candidate_count"] == v["validated_candidates"]
    assert v["candidate_point_checks"] == len(choices)*v["tested_points"]
    winner = min(choices,key=lambda c:(c["fallback"],c["tree_nodes"],c["plane_references"],c["finite_prisms"],c["axis"]))
    assert choices[candidates["selected_index"]] == winner
    assert winner["axis"] == candidates["selected_axis"] == m["selected_axis"] == manifest["selected_axis"]
    assert winner["tree_nodes"] == m["tree_nodes"]
    assert all(c["certificate_passed"] and c["probes_passed"] for c in choices)
    assert m["core_runtime_ms"] >= m["all_candidates_ms"]+m["selection_ms"]
    assert m["all_candidates_ms"] >= sum(c["build_ms"] for c in choices)
    assert abs(m["tree_build_simplify_ms"]-sum(c["tree_ms"] for c in choices)) < 1e-8
    supports = load(output/"boundary_supports.json")
    patches = load(output/"boundary_patches.json")["patches"]
    for plane in supports["planes"]:
        assert plane["source_patch_ids"]
        for p in plane["source_patch_ids"]:
            patch = patches[p]
            assert patch["axis"] == plane["axis"] and patch["plane_index"] == plane["grid_coordinate"]
    if not manifest["fallback"]:
        assert manifest["requires_zero_guard"] and seen_occupied_zero
        assert len(supports["leaf_sources"]) == len(tree["leaves"])
        for link in supports["leaf_sources"]:
            plane = supports["planes"][link["support_id"]]
            leaf = tree["leaves"][link["leaf_id"]]
            assert leaf["normal"] == [link["coefficient_sign"] if a==plane["axis"] else 0 for a in range(3)]
            assert leaf["offset"] == link["coefficient_sign"]*plane["coordinate_m"]
    else:
        assert not manifest["requires_zero_guard"] and not supports["leaf_sources"]
    return seen_free_zero


with tempfile.TemporaryDirectory(prefix="boundary_tree_pipeline_") as temporary:
    root = Path(temporary)
    shapes = {
        "single":[(.5,.5,.5)],
        "u":[(x+.5,y+.5,.5) for x,y in itertools.product(range(3),repeat=2) if x!=1 or y==0],
        "shell":[(x+.5,y+.5,z+.5) for x,y,z in itertools.product(range(3),repeat=3) if (x,y,z)!=(1,1,1)],
        "edge":[(.5,.5,.5),(1.5,1.5,.5)],
        "negative_decimal":[(-.015,-.015,-.015),(.015,.015,.015)],
    }
    def run(name, points, flags, good=True):
        cloud = root/(name+".xyz")
        cloud.write_text("\n".join(" ".join(map(str,p)) for p in points))
        out = root/name
        p = subprocess.run([str(program),str(cloud),str(out),"--pipeline","octree","--voxel-size",
                            ".03" if name.startswith("negative_decimal") else "1","--validation-samples","20"]+flags,
                           capture_output=True,text=True,timeout=30)
        assert (p.returncode==0)==good,p.stdout+p.stderr
        return out
    for name,points in shapes.items():
        out = run(name,points,["--tree-source","boundary"])
        free_zero = verify(out,points)
        if name=="u":
            assert free_zero  # Strict <0 alone would incorrectly remove internal free seams.
        outputs = {}
        for axis in ("x","y","z","best-axis"):
            outputs[axis] = run(name+"_"+axis,points,["--tree-source","boundary","--boundary-sweep",axis])
            verify(outputs[axis],points)
        filename = "free_closure_tree_simplified.json"
        assert (out/filename).read_bytes() == (outputs["x"]/filename).read_bytes()
        winner = load(outputs["best-axis"]/"strict_free_space.json")["selected_axis"]
        assert (outputs["best-axis"]/filename).read_bytes() == (outputs[winner]/filename).read_bytes()
    for option in ("--boundary-max-cells","--boundary-max-events","--boundary-max-prisms"):
        for axis in ("x","y","z","best-axis"):
            out = run("cap_"+option+axis,shapes["edge"],["--tree-source","boundary",option,"1","--boundary-sweep",axis])
            assert load(out/"strict_free_space.json")["fallback"]
            verify(out,shapes["edge"])
    old = run("old",shapes["u"],[])
    explicit = run("explicit",shapes["u"],["--tree-source","volume"])
    for name in ("logic_tree_raw.json","logic_tree_simplified.json","convex_clusters.json"):
        assert (old/name).read_bytes()==(explicit/name).read_bytes()
    for flags in (["--boundary-sweep","w"],["--octree-approx","hull"],["--boundary-max-cells","0"],
                  ["--tree-source","volume"],["--alpha","1"],["--boundary-merge-passes","0"]):
        run("bad",shapes["single"],["--tree-source","boundary"]+flags,False)
print("Boundary strict-query exports, provenance, strata and fallback passed")
