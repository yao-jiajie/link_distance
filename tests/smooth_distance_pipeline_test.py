"""Independent reference distance and exported LSE function checks."""
import csv
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

binary=Path(sys.argv[1]).resolve()
def load(p): return json.loads(p.read_text())
def evaluate(tree,p,beta=None,with_gradient=False):
    assert not with_gradient or beta is not None
    values={}; gradients={}
    for n in tree["nodes"]:
        if n["type"]=="leaf":
            f=tree["leaves"][n["leaf_id"]]; v=sum(a*b for a,b in zip(f["normal"],p))-f["offset"]
            if with_gradient: gradients[n["id"]]=f["normal"]
        else:
            a=[values[c] for c in n["children"]]; anchor=(max if n["type"]=="max" else min)(a)
            if beta is None: v=anchor
            else:
                sign=1 if n["type"]=="max" else -1
                weights=[math.exp(sign*beta*(x-anchor)) for x in a]
                total=math.fsum(weights)
                v=anchor+sign*math.log(total if sign==1 else total/len(a))/beta
                if with_gradient:
                    # Independent formula: include ALL occurrences, normalize
                    # first, then combine (not the C++ anchor-seeded algorithm).
                    gradients[n["id"]]=[math.fsum(w/total*gradients[c][axis]
                        for w,c in zip(weights,n["children"])) for axis in range(3)]
        values[n["id"]]=v
    value=-values[tree["root_id"]]
    return (value,[-v for v in gradients[tree["root_id"]]]) if with_gradient else value

def path_weight(tree):
    weights={}
    for n in tree["nodes"]:
        weights[n["id"]]=0 if n["type"]=="leaf" else max(weights[c] for c in n["children"])+math.log(len(n["children"]))
    return weights[tree["root_id"]]

with tempfile.TemporaryDirectory(prefix="smooth_distance_") as directory:
    root=Path(directory)
    grid=list(itertools.product(range(3),repeat=3))
    cases={"cube":[(0,0,0)],"l":[(0,0,0),(1,0,0),(0,1,0)],
           "u":[(x,y,0) for x,y in itertools.product(range(3),repeat=2) if x!=1 or y==0],
           "cavity":[p for p in grid if p!=(1,1,1)],"negative":[(-2,0,-1),(-1,0,-1),(-2,1,-1)],
           "edge":[(0,0,0),(1,1,0)],"vertex":[(0,0,0),(1,1,1)],"adjacent":[(0,0,0),(1,0,0)],
           "thin_wall":[(x,y,0) for x,y in itertools.product(range(3),repeat=2)]}
    def run(name,keys,extra=(),good=True):
        cloud=root/(name+".xyz"); cloud.write_text("\n".join(" ".join(str(.04*(v+.5)) for v in p) for p in keys))
        out=root/name
        p=subprocess.run([str(binary),str(cloud),str(out),"--pipeline","octree","--tree-source","boundary",
            "--boundary-expression","direct","--voxel-size",".04","--validation-samples","20",*extra],capture_output=True,text=True,timeout=60)
        assert (p.returncode==0)==good,p.stdout+p.stderr
        if not good: assert not (out/"smooth_distance.json").exists() and not (out/"strict_free_space.json").exists()
        return out
    for name,keys in cases.items():
        exact=run(name+"_exact",keys); out=run(name,keys,["--distance-field","lse","--lse-error",".001",
            "--validation-cache","exact","--query-benchmark","true"])
        for f in ("free_direct_tree_raw.json","free_direct_tree_simplified.json","strict_free_space.json","boundary.off","boundary_sign_diagnostics.json"):
            assert (exact/f).read_bytes()==(out/f).read_bytes(),f
        assert not (exact/"smooth_distance.json").exists()
        field=load(out/"smooth_distance.json"); tree=load(out/field["tree"]); d=load(out/"distance_validation.json"); m=load(out/"benchmark.json")[0]
        assert d["passed"] and not field["is_exact_euclidean_sdf"] and not field["zero_level_set_exact"] and not field["sdf_accuracy_certified"]
        weight=path_weight(tree)
        assert abs(weight-field["path_log_arity_weight"])<1e-13
        assert abs(weight/field["beta"]-field["scalar_error_bound_m"])<1e-15
        assert field["scalar_error_bound_m"]<=.001+1e-15
        metrics=d["metrics"]; assert metrics["false_free_count"]==metrics["upper_order_errors"]==metrics["bound_errors"]==0
        assert metrics["gradient_checks"]>0 and metrics["gradient_max_abs_error"]<1e-4
        assert m["smooth_value_query_us"]>0 and m["smooth_gradient_query_us"]>0 and m["reference_sdf_query_us"]>0
        assert m["smooth_binary_kernel"]==1 and m["smooth_binary_nodes"]>0
        assert m["smooth_reuses_compiled_program"]==m["execution_program_instances"]==1
        assert m["validation_cache_selected"]==1 and m["validation_cache_requested_auto"]==0
        assert m["sign_packed_propagation"]==1 and m["sign_propagation_batches"]==math.ceil(m["strata_count"]/64)
        assert metrics["probe_cache_enabled"]==1 and metrics["probe_cache_hits"]>0
        assert metrics["probe_evaluations"]+metrics["probe_cache_hits"]==metrics["distance_probes"]
        assert metrics["hard_evaluations"]+metrics["hard_cache_hits"]==metrics["probe_evaluations"]
        assert m["core_runtime_ms"]>=m["smooth_compile_ms"] and m["validation_ms"]>=m["sdf_validation_ms"]
        patches=load(out/"boundary_patches.json")["patches"]
        bounds=[([min(c[a] for c in f["corners_m"]) for a in range(3)],
                 [max(c[a] for c in f["corners_m"]) for a in range(3)]) for f in patches]
        rows=list(csv.DictReader((out/"distance_samples.csv").open()))
        for row in rows:
            p=[float(row[a+"_m"]) for a in "xyz"]
            filled=any(all(k*.04<=v<=(k+1)*.04 for k,v in zip(key,p)) for key in keys)
            unsigned=min(math.sqrt(sum(max(lo[a]-p[a],0,p[a]-hi[a])**2 for a in range(3))) for lo,hi in bounds)
            reference=-unsigned if filled else unsigned
            hard=evaluate(tree,p); smooth,gradient=evaluate(tree,p,field["beta"],True)
            assert abs(reference-float(row["reference_sdf_m"]))<1e-13
            assert abs(hard-float(row["hard_proxy_m"]))<1e-13 and abs(smooth-float(row["smooth_proxy_m"]))<1e-13
            assert smooth<=hard+1e-14 and hard-smooth<=field["scalar_error_bound_m"]+1e-13
            assert all(abs(v-float(row["gradient_"+a]))<1e-10 for a,v in zip("xyz",gradient)),(name,p,gradient,row)
            assert math.sqrt(math.fsum(v*v for v in gradient))<=1+1e-12
        # Independent exported evaluator finite differences, including roots
        # near ties: no access to in-process C++ values or derivative workspace.
        for row in rows[::max(1,len(rows)//12)]:
            p=[float(row[a+"_m"]) for a in "xyz"]; step=1e-3/field["beta"]
            for axis,a in enumerate("xyz"):
                lo=p.copy(); hi=p.copy(); lo[axis]-=step; hi[axis]+=step
                fd=(evaluate(tree,hi,field["beta"])-evaluate(tree,lo,field["beta"]))/(hi[axis]-lo[axis])
                assert abs(fd-float(row["gradient_"+a]))<2e-6,(name,p,axis,fd,row)
        if name=="cube":
            assert metrics["hard_sdf_error_max_m"]>0 and metrics["zero_ray_max_offset_m"]>0, metrics
    beta_out=run("beta",cases["cube"],["--distance-field","lse","--lse-beta","1000"])
    beta_field=load(beta_out/"smooth_distance.json"); beta_tree=load(beta_out/beta_field["tree"])
    assert beta_field["beta"]==1000 and abs(beta_field["scalar_error_bound_m"]-path_weight(beta_tree)/1000)<1e-15
    for row in csv.DictReader((beta_out/"distance_samples.csv").open()):
        value,gradient=evaluate(beta_tree,[float(row[a+"_m"]) for a in "xyz"],1000,True)
        assert abs(value-float(row["smooth_proxy_m"]))<1e-13
        assert all(abs(v-float(row["gradient_"+a]))<1e-10 for a,v in zip("xyz",gradient))
    no_extra=run("no_extra",cases["cube"],["--distance-field","lse","--lse-error",".001",
        "--sdf-gradient-samples","0","--sdf-zero-rays","0"])
    no_extra_metrics=load(no_extra/"distance_validation.json")["metrics"]
    no_extra_benchmark=load(no_extra/"benchmark.json")[0]
    assert no_extra_metrics["gradient_checks"]==no_extra_metrics["zero_rays_tested"]==0
    assert no_extra_metrics["distance_probes"]>0 and load(no_extra/"validation_report.json")["strict_sign_certified"]
    assert no_extra_benchmark["validation_cache_requested_auto"]==1 and no_extra_benchmark["validation_cache_selected"]==0
    for i,extra in enumerate((["--distance-field","lse"], ["--lse-beta","1000"],
        ["--lse-kernel","binary"], ["--distance-field","lse","--lse-error",".001","--lse-kernel","unknown"],
        ["--distance-field","lse","--lse-beta","1000","--lse-error",".01"],
        ["--distance-field","lse","--lse-beta","-1"], ["--distance-field","lse","--lse-beta","nan"],
        ["--distance-field","lse","--lse-error",".001","--validation-cache","unknown"],
        ["--sign-propagation","unknown"],
        ["--distance-field","lse","--lse-error",".001","--sdf-max-work","1"],
        ["--distance-field","lse","--lse-error",".001","--query-benchmark","true","--boundary-max-query-ops","1"])):
        run("bad_"+str(i),cases["cube"],extra,False)
print("Smooth CLI: unchanged exact exports, independent finite-patch distances/LSE, gradients, parameters and budgets passed")
