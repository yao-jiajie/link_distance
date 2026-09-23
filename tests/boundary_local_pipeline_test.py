"""Independent world-coordinate verification of the local export contract."""
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile

binary = Path(sys.argv[1]).resolve()
def load(path):
    return json.loads(path.read_text())
def evaluate(tree,p):
    values = {}
    for n in tree["nodes"]:
        if n["type"]=="leaf":
            leaf = tree["leaves"][n["leaf_id"]]
            values[n["id"]] = sum(a*b for a,b in zip(leaf["normal"],p))-leaf["offset"]
        else:
            values[n["id"]] = (min if n["type"]=="min" else max)(values[c] for c in n["children"])
    return values[tree["root_id"]]

with tempfile.TemporaryDirectory(prefix="boundary_local_pipeline_") as directory:
    root = Path(directory)
    shapes = {
        "single":[(.5,.5,.5)],
        "gap":[(.5,.5,.5),(2.5,.5,.5)],
        "l_cavity":[(x+.5,y+.5,z+.5) for x,y,z in itertools.product(range(-2,2),range(-1,3),range(-1,2))
                    if not(z==0 and ((x==-1 and y==0) or (x==0 and y in (0,1))))],
    }
    def run(name,points,mode,extra=(),good=True):
        cloud = root/(name+".xyz"); cloud.write_text("\n".join(" ".join(map(str,p)) for p in points))
        out = root/name
        result = subprocess.run([str(binary),str(cloud),str(out),"--pipeline","octree","--tree-source","boundary",
                                 "--voxel-size","1","--boundary-sweep","y","--boundary-expression",mode,
                                 "--boundary-diagnostics","true","--validation-samples","20",*extra],
                                capture_output=True,text=True,timeout=30)
        assert (result.returncode==0)==good,result.stdout+result.stderr
        return out
    for name,points in shapes.items():
        flat = run(name+"_flat",points,"flat")
        local = run(name+"_local",points,"local")
        for suffix in ("raw","simplified"):
            f = "free_closure_tree_"+suffix+".json"
            assert (flat/f).read_bytes()==(local/f).read_bytes()
        manifest = load(local/"strict_free_space.json")
        report = load(local/"boundary_sign_diagnostics.json")
        metrics = load(local/"benchmark.json")[0]
        assert metrics["sign_packed_propagation"]==1
        assert metrics["sign_propagation_batches"]<=metrics["sign_strata_count"]
        d = report["local"]
        assert d["applied"] and d["closure_set_certified"] and d["unsafe_free_count"]==d["boundary_sign_errors"]==0
        assert d["strict_sign_certified"] == (name!="gap")
        assert manifest["requires_zero_guard"] != d["strict_sign_certified"]
        assert manifest["strict_sign_certified"] == d["strict_sign_certified"]
        if name=="gap":
            assert len(report["adjacency_faces"])==4 and all(f["exterior"] for f in report["adjacency_faces"])
            assert d["witnesses"] and sum(d["free_zeros_by_dimension"])>0
        if name=="l_cavity":
            assert metrics["logical_merge_count"]>0 and sum(report["flat"]["free_zeros_by_dimension"])>0
        tree = load(local/manifest["tree"]); raw = load(local/manifest["raw_tree"])
        reference = load(local/manifest["occupancy_reference"])
        boxes = [(n["lower"],n["upper"]) for n in reference["nodes"] if n["leaf"]]
        patches = load(local/"boundary_patches.json")["patches"]
        coords=[]
        for axis in range(3):
            c=sorted({b[s][axis] for b in boxes for s in (0,1)})
            coords.append([c[0]-1]+sorted(c+[(a+b)/2 for a,b in zip(c,c[1:])])+[c[-1]+1])
        for p in itertools.product(*coords):
            occupied=any(all(lo<=v<=hi for v,lo,hi in zip(p,*b)) for b in boxes)
            boundary=any(p[f["axis"]]==f["corners_m"][0][f["axis"]] and
                         all(min(c[a] for c in f["corners_m"])<=p[a]<=max(c[a] for c in f["corners_m"])
                             for a in range(3)) for f in patches)
            q=evaluate(tree,p); assert q==evaluate(raw,p)
            if d["strict_sign_certified"]:
                assert (-1 if q<0 else 1 if q>0 else 0)==(0 if boundary else 1 if occupied else -1),(name,p,q)
            free=q<0 or (q==0 and manifest["requires_zero_guard"] and not occupied)
            assert free==(not occupied)
        assert metrics["core_runtime_ms"]>=metrics["all_candidates_ms"]+metrics["selection_ms"]
        assert metrics["constructed_tree_nodes"]==len(raw["nodes"])
        assert metrics["expanded_node_references"]>=metrics["simplified_tree_nodes"]
    fallback=run("fallback",shapes["gap"],"local",["--boundary-max-cells","1"])
    m=load(fallback/"strict_free_space.json")
    assert m["fallback"] and m["classifier"]=="root_gt_zero" and not m["strict_sign_certified"]
    capped=run("node_cap",shapes["gap"],"local",["--boundary-max-local-nodes","1"])
    assert load(capped/"boundary_sign_diagnostics.json")["stop_reason"]=="initial_expanded_node_limit"
    run("sign_cap",shapes["gap"],"local",["--boundary-max-sign-ops","1"],False)
    run("bad_sign_propagation",shapes["single"],"local",["--sign-propagation","unknown"],False)
    run("bad_mode",shapes["single"],"unknown",good=False)
print("Local CLI exports, strict signs, safe incomplete mode and work caps passed")
