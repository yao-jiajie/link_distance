"""Check the synthetic geometry, sampling relationships and reproducible XYZ files."""
from collections import Counter
import hashlib
import math
from pathlib import Path
import sys
import tempfile

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO/"scripts"))
from generate_point_cloud_cases import make_cases, write_cases
from benchmark_octree_phase4 import read_xyz


def on_box(point, box):
    lo,hi = box
    return all(a<=v<=b for a,v,b in zip(lo,point,hi)) and any(v in (a,b) for a,v,b in zip(lo,point,hi))


cases = make_cases()
assert len(cases)==20 and cases==make_cases()
for name,case in cases.items():
    points = case["points"]
    assert points and all(len(p)==3 and all(math.isfinite(v) for v in p) for p in points)
    assert all(-.4<p[0]<.4 and -.4<p[1]<.4 and .3<p[2]<1.2 for p in points),name
    if case["sampling"]=="box_union_surface_area":
        boxes = case["shape"]["boxes"]
        assert all(any(on_box(p,b) for b in boxes) for p in points),name
for p in cases["cylinder_surface"]["points"]:
    radial = math.hypot(p[0],p[1]); z = abs(p[2]-.72)
    assert radial<=.14+1e-14 and z<=.16+1e-14
    assert abs(radial-.14)<1e-14 or abs(z-.16)<1e-14
for x,y,z in cases["torus_surface"]["points"]:
    assert abs((math.hypot(x,y)-.16)**2+(z-.74)**2-.06**2)<1e-14
for p in cases["hollow_box_surfaces"]["points"]:
    shape = cases["hollow_box_surfaces"]["shape"]
    assert any(on_box(p,shape[side]["boxes"][0]) for side in ("outer","inner"))
assert all(p[2]==.74 for p in cases["plane_patch"]["points"])
dense = cases["ellipsoid_dense"]["points"]
clean = cases["ellipsoid_surface"]["points"]
assert clean==dense[:2048] and cases["ellipsoid_sparse"]["points"]==dense[:256]
for x,y,z in dense:
    assert abs((x/.18)**2+(y/.13)**2+((z-.75)/.24)**2-1)<1e-14
assert len(cases["ellipsoid_partial"]["points"])==1024
assert all(p[0]>=0 for p in cases["ellipsoid_partial"]["points"])
assert Counter(cases["ellipsoid_duplicates"]["points"])==Counter({p: 3*n for p,n in Counter(clean).items()})
assert cases["ellipsoid_duplicates"]["points"]!=clean*3
for a,b,c in zip(clean,cases["ellipsoid_noise_2mm"]["points"],cases["ellipsoid_noise_8mm"]["points"]):
    assert all(abs((z-x)-4*(y-x))<1e-14 for x,y,z in zip(a,b,c))
outliers = cases["ellipsoid_outliers"]
assert outliers["points"][:2048]==clean and len(outliers["points"])==2054
for x,y,z in outliers["outliers_m"]:
    assert (x/.18)**2+(y/.13)**2+((z-.75)/.24)**2>1
assert cases["camera"]["points"]==read_xyz(REPO/"data/realistic_sparse_camera.xyz")

with tempfile.TemporaryDirectory(prefix="point_cloud_fixtures_") as temporary:
    root = Path(temporary)
    a = write_cases(root/"a"); b = write_cases(root/"b")
    assert a==b and a["length_unit"]=="m" and len(a["cases"])==20
    assert (root/"a/preview.svg").read_bytes()==(root/"b/preview.svg").read_bytes()
    for case in a["cases"]:
        path = root/"a"/case["file"]
        assert path.read_bytes()==(root/"b"/case["file"]).read_bytes()
        assert hashlib.sha256(path.read_bytes()).hexdigest()==case["sha256"]
        assert read_xyz(path)==cases[case["name"]]["points"]
    try:
        write_cases(root/"a")
    except FileExistsError:
        pass
    else:
        raise AssertionError("Generator overwrote an existing input directory")
    for names in (["unknown"],["camera","camera"],[]):
        try:
            write_cases(root/"bad",names)
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid fixture selection accepted")
print("Point-cloud inputs: 20 deterministic meter-scale shapes, density/noise relationships and XYZ round trips passed")
