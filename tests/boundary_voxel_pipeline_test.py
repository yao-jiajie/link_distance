"""Compare unchanged flat trees and independently query both occupancy bundles."""
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile

binary=Path(sys.argv[1]).resolve()
def load(path): return json.loads(path.read_text())
def evaluate(tree,p):
    values={}
    for n in tree["nodes"]:
        if n["type"]=="leaf":
            leaf=tree["leaves"][n["leaf_id"]]
            values[n["id"]]=sum(x*y for x,y in zip(p,leaf["normal"]))-leaf["offset"]
        else: values[n["id"]]=(min if n["type"]=="min" else max)(values[c] for c in n["children"])
    return values[tree["root_id"]]

with tempfile.TemporaryDirectory(prefix="boundary_voxel_pipeline_") as directory:
    root=Path(directory)
    shapes={"single":[(.5,.5,.5)],"gap":[(.5,.5,.5),(2.5,.5,.5)],
            "edge":[(.5,.5,.5),(1.5,1.5,.5)],"vertex":[(.5,.5,.5),(1.5,1.5,1.5)],
            "u":[(x+.5,y+.5,.5) for x,y in itertools.product(range(3),repeat=2) if x!=1 or y==0],
            "hollow":[tuple(x+.5 for x in p) for p in itertools.product(range(3),repeat=3) if p!=(1,1,1)],
            "decimal":[(-.015,-.015,-.015),(.015,.015,.015)]}
    def run(name,points,mode,axis="y",extra=(),good=True):
        cloud=root/(name+".xyz"); cloud.write_text("\n".join(" ".join(map(str,p)) for p in points))
        out=root/name
        result=subprocess.run([str(binary),str(cloud),str(out),"--pipeline","octree","--tree-source","boundary",
                               "--voxel-size",".03" if name.startswith("decimal") else "1",
                               "--boundary-occupancy",mode,"--boundary-sweep",axis,"--boundary-diagnostics","true",
                               "--boundary-validation-reference","voxels","--validation-samples","20"]+list(extra),
                              capture_output=True,text=True,timeout=60)
        assert (result.returncode==0)==good,result.stdout+result.stderr
        return out
    def compare(a,b):
        ma,mb=(load(p/"benchmark.json")[0] for p in (a,b))
        manifest=load(b/"strict_free_space.json")
        assert mb["occupancy_backend"]=="voxels" and mb["validation_reference"]=="fine_voxels"
        assert manifest["occupancy_backend"]=="voxels" and manifest["occupancy_schema"]==2
        for k in ("fallback","tree_nodes","plane_references","raw_points","occupied_voxels","tested_points",
                  "false_free","missed_free","scalar_mismatches","raw_points_classified_free","candidate_count",
                  "artificial_zero_faces","artificial_zero_edges","artificial_zero_vertices","strict_sign_certified"):
            assert ma[k]==mb[k],(k,ma[k],mb[k])
        assert mb["octree_built"]==mb["fallback"]
        if not mb["fallback"]:
            assert mb["octree_nodes"]==mb["lazy_fallback_ms"]==mb["baseline_partition_ms"]==0
        else: assert mb["octree_nodes"]>0 and mb["lazy_fallback_ms"]>0
        for name in (manifest["tree"],manifest["raw_tree"],"free_prisms.json","boundary_supports.json","boundary.off"):
            assert (a/name).read_bytes()==(b/name).read_bytes(),name
        ref=load(b/manifest["occupancy_reference"])
        assert ref["root_id"] is None and ref["nodes"]==[] and ref["closed_boxes"]
        old_ref=load(a/"occupancy_reference.json")
        assert ref["occupied_voxel_indices"]==old_ref["occupied_voxel_indices"]
        h=ref["voxel_size"]
        boxes=[([k*h for k in key],[(k+1)*h for k in key]) for key in ref["occupied_voxel_indices"]]
        assert sorted(boxes)==sorted((n["lower"],n["upper"]) for n in old_ref["nodes"] if n["leaf"])
        tree=load(b/manifest["tree"]); raw=load(b/manifest["raw_tree"])
        coords=[]
        for axis in range(3):
            cs=sorted({bounds[axis] for box in boxes for bounds in box})
            coords.append([cs[0]-h]+sorted(cs+[(x+y)/2 for x,y in zip(cs,cs[1:])])+[cs[-1]+h])
        for p in itertools.product(*coords):
            occupied=any(all(lo<=x<=hi for x,lo,hi in zip(p,*box)) for box in boxes)
            q=evaluate(tree,p)
            assert q==evaluate(raw,p)
            actual=q>0 if manifest["fallback"] else q<0 or (q==0 and not occupied)
            assert actual== (not occupied),(p,q)
        if (a/"query_samples.json").exists():
            assert (a/"query_samples.json").read_bytes()==(b/"query_samples.json").read_bytes()
            assert ma["query_zero_guard_fraction"]==mb["query_zero_guard_fraction"]
            assert ma["artificial_seam_query_count"]==mb["artificial_seam_query_count"]
    for name,points in shapes.items():
        for axis in ("x","y","z","best-axis"):
            extra=("--query-benchmark","true") if name=="gap" and axis=="y" else ()
            compare(run(name+axis+"old",points,"octree",axis,extra),run(name+axis+"new",points,"voxels",axis,extra))
    for cap in ("--boundary-max-cells","--boundary-max-events","--boundary-max-prisms"):
        for axis in ("x","best-axis"):
            compare(run("old"+cap+axis,shapes["edge"],"octree",axis,(cap,"1")),
                    run("new"+cap+axis,shapes["edge"],"voxels",axis,(cap,"1")))
    # No hidden hierarchy in normal construction/validation: depth zero is valid.
    out=run("nohierarchy",shapes["gap"],"voxels",extra=("--max-depth","0"))
    assert load(out/"benchmark.json")[0]["octree_nodes"]==0
    run("fallbackdepth",shapes["gap"],"voxels",extra=("--max-depth","0","--boundary-max-cells","1"),good=False)
    run("badlocal",shapes["single"],"voxels",extra=("--boundary-expression","local"),good=False)
    run("badmode",shapes["single"],"invalid",good=False)
print("Voxel flat: byte-identical geometry, closed occupancy exports, paired queries and lazy fallback passed")
