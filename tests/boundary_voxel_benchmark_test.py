"""Paired benchmark must retain byte equality and explicitly label lazy fallback."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

repo=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix="boundary_voxel_benchmark_") as directory:
    output=Path(directory)/"results"
    subprocess.run([sys.executable,str(repo/"scripts/benchmark_boundary_voxels.py"),str(output),
                    "--build-dir",sys.argv[1],"--repeats","1","--warmups","0","--query-repeats","1",
                    "--cases","single","aligned_gap","sparse_disconnected"],check=True,timeout=100)
    rows=json.loads((output/"summary.json").read_text())
    assert len(rows)==6
    for name in ("single","aligned_gap","sparse_disconnected"):
        group={r["occupancy"]:r for r in rows if r["case"]==name}
        a,b=group["octree"],group["voxels"]
        assert a["passed"] and b["passed"] and a["tree_sha256"]==b["tree_sha256"]
        assert b["byte_identical_geometry"] and b["byte_identical_query_sets"]
        assert (b["octree_nodes"]>0)==b["fallback"]==(name=="sparse_disconnected")
        assert b["core_ms"]<=b["build_validation_ms"] and b["rect_partition_ms"]==0
        assert a["artificial_zeros_fev"]==b["artificial_zeros_fev"]
        if name=="aligned_gap":
            assert b["seam_query_count"]>0 and b["seam_query_us"]>0
    assert (output/"report.md").is_file()
print("Octree removal paired benchmark, query corpora and fallback reporting passed")
