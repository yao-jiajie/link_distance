"""End-to-end equivalence checks for validation caching and packed signs."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

repo=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix="validation_acceleration_") as temporary:
    out=Path(temporary)/"results"
    subprocess.run([sys.executable,str(repo/"scripts/benchmark_validation_acceleration.py"),str(out),
        "--build-dir",sys.argv[1],"--cases","l_voxels","camera","--warmups","0","--repeats","1",
        "--validation-samples","20"],check=True,timeout=100)
    rows=json.loads((out/"summary.json").read_text())
    assert len(rows)==6 and all(row["passed"] and row["same_outputs"] for row in rows)
    for case in {row["case"] for row in rows}:
        modes={row["mode"]:row for row in rows if row["case"]==case}
        baseline=modes["scalar_uncached"]; cached=modes["scalar_cached"]; packed=modes["packed_cached"]
        assert baseline["sign_batches"]==baseline["strata"] and baseline["probe_cache_hits"]==baseline["hard_cache_hits"]==0
        assert packed["sign_batches"]==(packed["strata"]+63)//64
        assert cached["probe_cache_hits"]>0 and cached["probe_evaluations"]<cached["distance_probes"]
        assert cached["hard_cache_hits"]>0 and packed["hard_cache_hits"]>0
        assert baseline["strata"]==cached["strata"]==packed["strata"]
        assert baseline["distance_probes"]==cached["distance_probes"]==packed["distance_probes"]
    camera=out/"camera"; automatic=camera/"automatic"
    subprocess.run([str(Path(sys.argv[1]).resolve()/"build_halfspace_tree"),str(camera/"input.xyz"),str(automatic),
        "--pipeline","octree","--tree-source","boundary","--boundary-expression","direct","--boundary-query","dag",
        "--voxel-size",".04","--validation-samples","20","--distance-field","lse","--lse-error",".001",
        "--execution-sharing","structural","--lse-kernel","binary"],check=True,capture_output=True,text=True,timeout=60)
    automatic_metrics=json.loads((automatic/"benchmark.json").read_text())[0]
    assert automatic_metrics["validation_cache_requested_auto"]==automatic_metrics["validation_cache_selected"]==1
    assert automatic_metrics["distance_probe_cache_hits"]>0 and automatic_metrics["sign_packed_propagation"]==1
print("Validation acceleration: exact probe reuse and 64-way all-strata propagation preserve outputs")
