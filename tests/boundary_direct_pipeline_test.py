"""Independent exported-tree verification; no occupancy guard or flat reference required."""
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile

binary=Path(sys.argv[1]).resolve()
def load(path):
    return json.loads(path.read_text())
def evaluate(tree,p):
    values={}
    for n in tree["nodes"]:
        assert n["type"] in ("leaf","min","max")
        if n["type"]=="leaf":
            leaf=tree["leaves"][n["leaf_id"]]
            values[n["id"]]=sum(a*b for a,b in zip(leaf["normal"],p))-leaf["offset"]
        else:
            values[n["id"]]=(min if n["type"]=="min" else max)(values[c] for c in n["children"])
    return values[tree["root_id"]]

with tempfile.TemporaryDirectory(prefix="boundary_direct_cli_") as directory:
    root=Path(directory)
    grid=list(itertools.product(range(3),repeat=3))
    shapes={"cuboid":grid, "l":[(0,0,0),(1,0,0),(0,1,0)],
        "u":[(x,y,0) for x,y in itertools.product(range(3),repeat=2) if x!=1 or y==0],
        "stair":[(x,y,0) for x,y in itertools.product(range(3),repeat=2) if y<=x],
        "cavity":[p for p in grid if p!=(1,1,1)], "disconnected":[(0,0,0),(3,2,1)],
        "gap":[(0,0,0),(2,0,0)], "edge":[(0,0,0),(1,1,0)], "vertex":[(0,0,0),(1,1,1)],
        "negative":[(-2,-1,-1),(-1,-1,-1),(-2,0,-1)]}
    def run(name,keys,mode="direct",extra=(),good=True):
        cloud=root/(name+".xyz"); cloud.write_text("\n".join(" ".join(str(v+.5) for v in p) for p in keys))
        out=root/name
        p=subprocess.run([str(binary),str(cloud),str(out),"--pipeline","octree","--tree-source","boundary",
                 "--boundary-expression",mode,"--voxel-size","1","--validation-samples","20",
                 "--boundary-occupancy","voxels",*extra],capture_output=True,text=True,timeout=60)
        assert (p.returncode==0)==good,p.stdout+p.stderr
        if not good: assert not (out/"strict_free_space.json").exists()
        return out
    for name,keys in shapes.items():
        out=run(name,keys)
        manifest=load(out/"strict_free_space.json"); d=load(out/"boundary_sign_diagnostics.json"); m=load(out/"benchmark.json")[0]
        assert manifest["classifier"]=="root_lt_zero" and not manifest["requires_zero_guard"]
        assert not manifest["occupancy_reference_required_for_query"] and manifest["strict_sign_certified"]
        assert manifest["query_backend"]==m["query_backend"]=="dag"
        assert m["execution_structural_sharing"]==1
        assert m["sign_packed_propagation"]==1 and m["sign_propagation_batches"]==(m["strata_count"]+63)//64
        assert m["query_node_visits"]+m["compiled_shared_nodes"]==m["tree_nodes"] and m["evaluator_compile_ms"]>0
        assert m["compiled_program_bytes"]>0 and m["query_workspace_bytes"]==8*m["query_node_visits"]
        recursive=run(name+"_recursive",keys,extra=["--boundary-query","recursive","--boundary-query-cross-check","true"])
        rm=load(recursive/"benchmark.json")[0]
        assert rm["scalar_cross_check_errors"]==0 and rm["scalar_cross_checks"]==rm["tested_points"]
        assert rm["query_node_visits"]==rm["expanded_node_references"]
        for f in ("free_direct_tree_raw.json","free_direct_tree_simplified.json","direct_splits.json","boundary_sign_diagnostics.json"):
            assert (out/f).read_bytes()==(recursive/f).read_bytes()
        assert d["applied"] and d["strict_sign_certified"] and not d["witnesses"]
        assert sum(d["free_zeros_by_dimension"]+d["occupied_zeros_by_dimension"])==0
        assert m["octree_nodes"]==m["free_rectangles"]==m["free_prisms"]==m["exterior_free_branches"]==m["simplifier_applied"]==0
        for f in ("free_prisms.json","axis_candidates.json","free_closure_tree_raw.json"): assert not (out/f).exists()
        tree=load(out/manifest["tree"]); raw=load(out/manifest["raw_tree"])
        assert raw["nodes"]==tree["nodes"] and raw["leaves"]==tree["leaves"]
        assert len(tree["nodes"])==m["tree_nodes"]==m["constructed_tree_nodes"]==m["simplified_tree_nodes"]
        reg=load(out/"boundary_supports.json"); splits=load(out/"direct_splits.json")
        patches=load(out/"boundary_patches.json")["patches"]
        assert splits["includes_unbounded_cells"]
        for s in reg["supports"]:
            assert s["patches"] and all(patches[i]["axis"]==s["axis"] and
                patches[i]["corners_m"][0][s["axis"]]==s["coordinate_m"] for i in s["patches"])
        for source in reg["leaf_sources"]:
            leaf=tree["leaves"][source["leaf_id"]]; support=reg["supports"][source["support"]]
            assert leaf["normal"]==[source["sign"] if a==support["axis"] else 0 for a in range(3)]
            assert leaf["offset"]==source["sign"]*support["coordinate_m"]
        coords=[]
        for a in range(3):
            c=sorted({p[a]+offset for p in keys for offset in (0,1)})
            coords.append([c[0]-1]+sorted(c+[(x+y)/2 for x,y in zip(c,c[1:])])+[c[-1]+1])
        for p in itertools.product(*coords):
            occupied=any(all(k<=v<=k+1 for k,v in zip(key,p)) for key in keys)
            boundary=any(p[f["axis"]]==f["corners_m"][0][f["axis"]] and
                         all(min(c[a] for c in f["corners_m"])<=p[a]<=max(c[a] for c in f["corners_m"]) for a in range(3)) for f in patches)
            q=evaluate(tree,p); expected=0 if boundary else 1 if occupied else -1
            assert (int(q>0)-int(q<0))==expected,(name,p,q,expected)
        if name in ("l","gap"):
            flat=run(name+"_flat",keys,"flat",["--boundary-sweep","y"])
            ft=load(flat/"free_closure_tree_raw.json")
            zeros=[p for p in itertools.product(*coords) if evaluate(ft,p)==0 and
                   not any(all(k<=v<=k+1 for k,v in zip(key,p)) for key in keys)]
            assert zeros and all(evaluate(tree,p)<0 for p in zeros)
    for flag in ("--boundary-max-cells","--boundary-max-direct-nodes","--boundary-max-expanded-nodes",
                 "--boundary-max-split-checks","--boundary-max-direct-depth","--boundary-max-strata",
                 "--boundary-max-sign-ops","--boundary-max-probe-ops"):
        run("cap_"+flag[2:],shapes["l"],extra=[flag,"1"],good=False)
    run("query_cap",shapes["l"],extra=["--query-benchmark","true","--boundary-max-query-ops","1"],good=False)
    run("no_diagnostics",shapes["l"],extra=["--boundary-diagnostics","false"],good=False)
    run("no_sweep",shapes["l"],extra=["--boundary-sweep","x"],good=False)
    run("duplicate_mode",shapes["l"],extra=["--boundary-expression","flat"],good=False)
    run("bad_query_backend",shapes["l"],extra=["--boundary-query","unknown"],good=False)
    run("bad_sharing",shapes["l"],extra=["--execution-sharing","unknown"],good=False)
    run("bad_sign_propagation",shapes["l"],extra=["--sign-propagation","unknown"],good=False)
    run("bad_cross_check",shapes["l"],extra=["--boundary-query-cross-check","yes"],good=False)
    audit=run("dag_audit",shapes["l"],extra=["--boundary-query-cross-check","true"])
    am=load(audit/"benchmark.json")[0]
    assert am["scalar_cross_checks"]==am["tested_points"] and am["scalar_cross_check_errors"]==0
print("Direct CLI: strict signs on 10 shapes, real support provenance, flat-zero repair, no guard and failure caps passed")
