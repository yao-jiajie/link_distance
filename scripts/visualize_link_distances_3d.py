#!/usr/bin/env python3
"""Create a self-contained interactive 3D HTML for link-distance results."""

import argparse
import json
import math
from pathlib import Path


def load_xyz(path):
    points = []
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        fields = line.split()
        if len(fields) != 3:
            raise ValueError(f"{path}:{line_number}: expected exactly three columns")
        point = [float(value) for value in fields]
        if not all(math.isfinite(value) for value in point):
            raise ValueError(f"{path}:{line_number}: nonfinite coordinate")
        points.append(point)
    if not points:
        raise ValueError(f"{path}: empty point cloud")
    return points


def resolve_cloud(result_path, result, override):
    if override is not None:
        return override
    candidate = Path(result["cloud"])
    if candidate.exists():
        return candidate
    candidate = result_path.parent / candidate
    if candidate.exists():
        return candidate
    raise FileNotFoundError("point cloud not found; pass --cloud explicitly")


def field_to_base_transform(result):
    """Return a point transform from the field frame to the robot base frame."""
    matrix = result.get("field_from_base")
    if matrix is None:
        matrix = [1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    if len(matrix) != 16 or not all(math.isfinite(float(value)) for value in matrix):
        raise ValueError("invalid field_from_base transform")
    matrix = [float(value) for value in matrix]
    translation = (matrix[3], matrix[7], matrix[11])

    def transform(point):
        delta = [float(point[index]) - translation[index] for index in range(3)]
        # field_T_base stores R; its rigid inverse uses R transpose.
        return [
            sum(matrix[row * 4 + column] * delta[row] for row in range(3))
            for column in range(3)
        ]

    return transform


HTML = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Robot link clearance 3D</title>
<style>
  *{box-sizing:border-box} html,body{width:100%;height:100%;margin:0;overflow:hidden;background:#0b1020;color:#e6edf7;font-family:Inter,system-ui,sans-serif}
  #app{display:grid;grid-template-columns:minmax(0,1fr) 330px;width:100%;height:100%}
  #viewport{position:relative;min-width:0;touch-action:none} canvas{display:block;width:100%;height:100%;cursor:grab} canvas.dragging{cursor:grabbing}
  #panel{overflow:auto;background:#121a2b;border-left:1px solid #26344f;padding:18px 16px}
  h1{font-size:18px;margin:0 0 5px}.sub{color:#93a4bf;font-size:12px;line-height:1.45;margin-bottom:14px}
  .controls{display:grid;grid-template-columns:1fr 1fr;gap:8px;margin-bottom:14px}.controls label{font-size:12px;color:#b8c5d9;background:#182238;border:1px solid #2b3b5a;border-radius:7px;padding:7px;cursor:pointer}
  .legend{display:grid;grid-template-columns:12px 1fr;gap:7px 8px;align-items:center;font-size:12px;color:#b8c5d9;margin:10px 0 16px}.swatch{width:12px;height:12px;border-radius:50%}
  .link{border:1px solid #293955;border-radius:8px;padding:9px 10px;margin:7px 0;background:#172137;cursor:pointer}.link:hover,.link.selected{border-color:#7ca9ff;background:#1c2b47}.link-head{display:flex;justify-content:space-between;gap:8px;font-size:13px}.distance{font-variant-numeric:tabular-nums;font-weight:650}.link-meta{color:#8496b3;font-size:11px;margin-top:4px}
  #help{position:absolute;left:14px;bottom:12px;background:#0e1729d9;border:1px solid #31415f;border-radius:8px;padding:8px 10px;color:#bac8dd;font-size:12px;pointer-events:none}
  #tooltip{display:none;position:absolute;z-index:5;max-width:260px;background:#111a2ddd;border:1px solid #54719f;border-radius:8px;padding:9px 11px;font-size:12px;line-height:1.5;pointer-events:none;box-shadow:0 8px 30px #0008}
  #badge{position:absolute;left:14px;top:12px;background:#0e1729d9;border:1px solid #31415f;border-radius:7px;padding:7px 9px;color:#b7c7df;font-size:12px;pointer-events:none}
  button{width:100%;border:1px solid #365079;background:#1a2944;color:#dbe8fa;border-radius:7px;padding:8px;cursor:pointer;margin:3px 0 10px}button:hover{background:#23385d}
  @media(max-width:800px){#app{grid-template-columns:1fr}#panel{display:none}}
</style>
</head>
<body><div id="app"><div id="viewport"><canvas id="canvas"></canvas><div id="badge"></div><div id="tooltip"></div><div id="help">Drag: orbit · Wheel: zoom · Right/Shift-drag: pan · Click sphere: inspect · Double-click: reset</div></div><aside id="panel"><h1 id="title"></h1><div class="sub" id="meta"></div><button id="reset">Reset camera</button><div class="controls"><label><input id="showPoints" type="checkbox" checked> Point cloud</label><label><input id="showVoxels" type="checkbox" checked> Voxel surface</label><label><input id="showSpheres" type="checkbox" checked> Spheres</label><label><input id="showLabels" type="checkbox" checked> Labels</label><label><input id="showAxes" type="checkbox" checked> Base axes</label></div><div class="legend"><span class="swatch" style="background:#4e90dc"></span><span>merged occupied-voxel boundary</span><span class="swatch" style="background:#d62828"></span><span>proxy collision ≤ 0</span><span class="swatch" style="background:#f4a261"></span><span>near zero</span><span class="swatch" style="background:#2a9d8f"></span><span>positive clearance</span></div><div id="links"></div></aside></div>
<script id="scene-data" type="application/json">__SCENE_DATA__</script>
<script>
"use strict";
const scene=JSON.parse(document.getElementById("scene-data").textContent);
const canvas=document.getElementById("canvas"),ctx=canvas.getContext("2d"),viewport=document.getElementById("viewport"),tooltip=document.getElementById("tooltip");
const checks={points:document.getElementById("showPoints"),voxels:document.getElementById("showVoxels"),spheres:document.getElementById("showSpheres"),labels:document.getElementById("showLabels"),axes:document.getElementById("showAxes")};
const vadd=(a,b)=>[a[0]+b[0],a[1]+b[1],a[2]+b[2]],vsub=(a,b)=>[a[0]-b[0],a[1]-b[1],a[2]-b[2]],vmul=(a,s)=>[a[0]*s,a[1]*s,a[2]*s];
const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2],cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const norm=a=>Math.hypot(a[0],a[1],a[2]),normalize=a=>{const n=norm(a)||1;return vmul(a,1/n)};
const voxelCorners=scene.voxelPatches.flatMap(p=>p.corners);const all=scene.points.concat(scene.spheres.map(s=>s.center),voxelCorners);let lo=[Infinity,Infinity,Infinity],hi=[-Infinity,-Infinity,-Infinity];
for(const p of all)for(let k=0;k<3;k++){lo[k]=Math.min(lo[k],p[k]);hi[k]=Math.max(hi[k],p[k])}for(const s of scene.spheres)for(let k=0;k<3;k++){lo[k]=Math.min(lo[k],s.center[k]-s.radius);hi[k]=Math.max(hi[k],s.center[k]+s.radius)}
const initialTarget=vmul(vadd(lo,hi),.5),extent=norm(vsub(hi,lo))||1;
let target=initialTarget.slice(),yaw=.8,pitch=.42,distance=extent*1.35,selected=null,drag=null,dragMoved=false,screenSpheres=[];
const fov=Math.PI/4;
function basis(){const cp=Math.cos(pitch),dir=[cp*Math.cos(yaw),cp*Math.sin(yaw),Math.sin(pitch)],position=vadd(target,vmul(dir,distance)),forward=normalize(vsub(target,position)),right=normalize(cross(forward,[0,0,1])),up=normalize(cross(right,forward));return{position,forward,right,up}}
function color(clearance){const scale=scene.clearanceScale,clamp=x=>Math.max(0,Math.min(1,x)),mix=(a,b,t)=>a.map((x,i)=>Math.round(x+t*(b[i]-x)));let rgb;if(clearance<=0)rgb=mix([244,162,97],[214,40,40],clamp(-clearance/scale));else rgb=mix([244,162,97],[42,157,143],clamp(clearance/scale));return `rgb(${rgb.join(",")})`}
function cameraProject(p,b,w,h,focal){const r=vsub(p,b.position),z=dot(r,b.forward);if(z<=1e-5)return null;return{x:w/2+focal*dot(r,b.right)/z,y:h/2-focal*dot(r,b.up)/z,z}}
function resize(){const dpr=Math.min(devicePixelRatio||1,2),w=viewport.clientWidth,h=viewport.clientHeight;if(canvas.width!==Math.round(w*dpr)||canvas.height!==Math.round(h*dpr)){canvas.width=Math.round(w*dpr);canvas.height=Math.round(h*dpr);ctx.setTransform(dpr,0,0,dpr,0,0)}return{w,h}}
function line3(a,b,color,width,cam,w,h,focal){const pa=cameraProject(a,cam,w,h,focal),pb=cameraProject(b,cam,w,h,focal);if(!pa||!pb)return;ctx.strokeStyle=color;ctx.lineWidth=width;ctx.beginPath();ctx.moveTo(pa.x,pa.y);ctx.lineTo(pb.x,pb.y);ctx.stroke()}
function render(){const{w,h}=resize(),cam=basis(),focal=h/(2*Math.tan(fov/2));ctx.clearRect(0,0,w,h);const bg=ctx.createLinearGradient(0,0,0,h);bg.addColorStop(0,"#0a1020");bg.addColorStop(1,"#111b2d");ctx.fillStyle=bg;ctx.fillRect(0,0,w,h);
 if(checks.axes.checked){const L=extent*.16;line3([0,0,0],[L,0,0],"#ff5b64",2,cam,w,h,focal);line3([0,0,0],[0,L,0],"#56d68b",2,cam,w,h,focal);line3([0,0,0],[0,0,L],"#5b9dff",2,cam,w,h,focal)}
 if(checks.voxels.checked){const faces=[];for(const patch of scene.voxelPatches){const q=patch.corners.map(p=>cameraProject(p,cam,w,h,focal));if(q.every(Boolean))faces.push({q,depth:q.reduce((sum,p)=>sum+p.z,0)/q.length,axis:patch.axis})}faces.sort((a,b)=>b.depth-a.depth);const fills=["rgba(78,144,220,.20)","rgba(87,164,224,.23)","rgba(102,181,213,.25)"];for(const face of faces){ctx.beginPath();ctx.moveTo(face.q[0].x,face.q[0].y);for(let i=1;i<face.q.length;i++)ctx.lineTo(face.q[i].x,face.q[i].y);ctx.closePath();ctx.fillStyle=fills[face.axis];ctx.fill();ctx.strokeStyle="rgba(126,190,244,.58)";ctx.lineWidth=.75;ctx.stroke()}}
 if(checks.points.checked){const projected=[];for(const p of scene.points){const q=cameraProject(p,cam,w,h,focal);if(q)projected.push(q)}projected.sort((a,b)=>b.z-a.z);ctx.fillStyle="rgba(151,169,195,.48)";for(const p of projected){const r=Math.max(.7,Math.min(2.2,2.4*distance/p.z));ctx.fillRect(p.x-r/2,p.y-r/2,r,r)}}
 screenSpheres=[];if(checks.spheres.checked){for(const s of scene.spheres){const p=cameraProject(s.center,cam,w,h,focal);if(p){p.r=Math.max(2,focal*s.radius/p.z);p.sphere=s;screenSpheres.push(p)}}screenSpheres.sort((a,b)=>b.z-a.z);for(const p of screenSpheres){const c=color(p.sphere.clearance),g=ctx.createRadialGradient(p.x-p.r*.28,p.y-p.r*.32,Math.max(1,p.r*.04),p.x,p.y,p.r);g.addColorStop(0,"rgba(255,255,255,.8)");g.addColorStop(.18,c);g.addColorStop(1,"rgba(8,16,29,.7)");ctx.globalAlpha=.7;ctx.fillStyle=g;ctx.beginPath();ctx.arc(p.x,p.y,p.r,0,Math.PI*2);ctx.fill();ctx.globalAlpha=1;ctx.strokeStyle=selected===p.sphere.id?"#ffe66d":c;ctx.lineWidth=selected===p.sphere.id?4:(p.sphere.nearest?2.4:1.1);ctx.stroke();if(checks.labels.checked&&p.sphere.nearest){ctx.font="12px system-ui";ctx.lineWidth=3;ctx.strokeStyle="#07101d";ctx.fillStyle="#e8f0fc";const t=`${p.sphere.link}: ${(p.sphere.linkDistance*1000).toFixed(1)} mm`;ctx.strokeText(t,p.x+5,p.y-p.r-5);ctx.fillText(t,p.x+5,p.y-p.r-5)}}}
 requestAnimationFrame(()=>{});
}
function reset(){target=initialTarget.slice();yaw=.8;pitch=.42;distance=extent*1.35;selected=null;tooltip.style.display="none";render()}
function focusSphere(id){const s=scene.spheres.find(x=>x.id===id);if(!s)return;target=s.center.slice();distance=Math.max(s.radius*8,extent*.28);selected=id;render()}
function inspect(event){let hit=null;for(const p of screenSpheres){const d=Math.hypot(event.offsetX-p.x,event.offsetY-p.y);if(d<=p.r&&(hit===null||p.z<hit.z))hit=p}if(!hit){selected=null;tooltip.style.display="none";render();return}const s=hit.sphere;selected=s.id;tooltip.style.display="block";tooltip.style.left=`${Math.min(viewport.clientWidth-270,event.offsetX+12)}px`;tooltip.style.top=`${Math.max(8,event.offsetY-18)}px`;tooltip.innerHTML=`<b>${s.link}</b><br>sphere ${s.id}, radius ${(s.radius*1000).toFixed(1)} mm<br>sphere clearance ${(s.clearance*1000).toFixed(3)} mm<br>link smooth ${(s.linkDistance*1000).toFixed(3)} mm`;render()}
canvas.addEventListener("mousedown",e=>{drag={x:e.clientX,y:e.clientY,pan:e.button===2||e.shiftKey};dragMoved=false;canvas.classList.add("dragging")});window.addEventListener("mouseup",()=>{drag=null;canvas.classList.remove("dragging")});window.addEventListener("mousemove",e=>{if(!drag)return;const dx=e.clientX-drag.x,dy=e.clientY-drag.y;if(Math.abs(dx)+Math.abs(dy)>1)dragMoved=true;drag.x=e.clientX;drag.y=e.clientY;if(drag.pan){const b=basis(),scale=distance/Math.max(300,viewport.clientHeight);target=vadd(target,vadd(vmul(b.right,-dx*scale),vmul(b.up,dy*scale)))}else{yaw-=dx*.008;pitch=Math.max(-1.45,Math.min(1.45,pitch+dy*.008))}tooltip.style.display="none";render()});
canvas.addEventListener("wheel",e=>{e.preventDefault();distance=Math.max(extent*.01,distance*Math.exp(e.deltaY*.001));render()},{passive:false});canvas.addEventListener("click",e=>{if(!dragMoved)inspect(e)});canvas.addEventListener("dblclick",reset);canvas.addEventListener("contextmenu",e=>e.preventDefault());window.addEventListener("resize",render);for(const c of Object.values(checks))c.addEventListener("change",render);document.getElementById("reset").addEventListener("click",reset);
document.getElementById("title").textContent=`${scene.robot} · 3D clearance`;document.getElementById("meta").textContent=`h=${(scene.voxelSize*1000).toFixed(2)} mm · ${scene.links.length} joint-link groups · ${scene.points.length} points · ${scene.occupancy.occupiedVoxels} occupied voxels · ${scene.voxelPatches.length} merged patches · ${scene.spheres.length} spheres · q=[${scene.configuration.join(", ")}]`;document.getElementById("badge").innerHTML=`common faces removed: ${scene.occupancy.commonFacesRemoved} · exposed faces: ${scene.occupancy.exposedFaces} → ${scene.voxelPatches.length} patches<br>env ε ≤ ${(scene.envError*1000).toPrecision(3)} mm · link ε ≤ ${(scene.linkError*1000).toPrecision(3)} mm`;
const links=document.getElementById("links");for(const link of scene.links){const row=document.createElement("div");row.className="link";const merged=link.sourceLinks.length>1?` · merged ${link.sourceLinks.join(" + ")}`:"";row.innerHTML=`<div class="link-head"><span>${link.name}</span><span class="distance" style="color:${color(link.distance)}">${(link.distance*1000).toFixed(2)} mm</span></div><div class="link-meta">${link.joint} · ${link.sphereCount} spheres · hard ${(link.hard*1000).toFixed(2)} mm${merged}</div>`;row.onclick=()=>{document.querySelectorAll(".link").forEach(x=>x.classList.remove("selected"));row.classList.add("selected");focusSphere(link.nearest)};links.appendChild(row)}
render();
</script></body></html>
"""


def main():
    parser = argparse.ArgumentParser(
        description="Create an offline interactive 3D link-clearance viewer"
    )
    parser.add_argument("result", type=Path, help="link_distance_node JSON")
    parser.add_argument("output", nargs="?", type=Path, help="output HTML")
    parser.add_argument("--cloud", type=Path, help="override XYZ path")
    parser.add_argument("--max-points", type=int, default=20000)
    parser.add_argument("--clearance-scale", type=float, default=0.1)
    args = parser.parse_args()
    if not args.result.is_file():
        parser.error(
            f"result JSON does not exist: {args.result}; run link_distance_node first"
        )
    if args.max_points <= 0:
        parser.error("max-points must be positive")
    if not math.isfinite(args.clearance_scale) or args.clearance_scale <= 0:
        parser.error("clearance-scale must be finite and positive")

    result = json.loads(args.result.read_text())
    if not result.get("spheres"):
        raise ValueError("result has no sphere diagnostics; rerun link_distance_node")
    cloud_path = resolve_cloud(args.result, result, args.cloud)
    field_to_base = field_to_base_transform(result)
    points = [field_to_base(point) for point in load_xyz(cloud_path)]
    if len(points) > args.max_points:
        step = len(points) / args.max_points
        points = [points[int(index * step)] for index in range(args.max_points)]

    links = []
    link_by_name = {}
    for link in result["links"]:
        item = {
            "name": link["name"],
            "joint": link["joint"],
            "sourceLinks": link["source_links"],
            "distance": link["distance_proxy_m"],
            "hard": link["hard_min_proxy_m"],
            "nearest": link["nearest_sphere"],
            "sphereCount": link["sphere_count"],
        }
        links.append(item)
        for source_link in item["sourceLinks"]:
            link_by_name[source_link] = item
    spheres = []
    for sphere in result["spheres"]:
        link = link_by_name.get(sphere["link"])
        spheres.append(
            {
                "id": sphere["id"],
                "link": link["name"] if link else sphere["link"],
                "sourceLink": sphere["link"],
                "center": field_to_base(sphere["center_field_m"]),
                "radius": sphere["radius_m"],
                "clearance": sphere["clearance_proxy_m"],
                "nearest": bool(link) and sphere["id"] == link["nearest"],
                "linkDistance": link["distance"] if link else sphere["clearance_proxy_m"],
            }
        )
    occupancy = result.get("voxel_occupancy")
    if not occupancy:
        raise ValueError("result has no voxel occupancy; rerun link_distance_node")
    voxel_patches = []
    for patch in occupancy.get("boundary_patches", []):
        corners = patch.get("corners_field_m", [])
        if len(corners) != 4 or any(len(corner) != 3 for corner in corners):
            raise ValueError("invalid merged voxel boundary patch")
        voxel_patches.append(
            {
                "axis": patch["axis"],
                "sign": patch["sign"],
                "corners": [field_to_base(corner) for corner in corners],
            }
        )
    scene = {
        "robot": result.get("robot", "robot"),
        "voxelSize": result["voxel_size_m"],
        "configuration": result.get("configuration", []),
        "points": points,
        "spheres": spheres,
        "links": links,
        "voxelPatches": voxel_patches,
        "occupancy": {
            "occupiedVoxels": occupancy["occupied_voxels"],
            "commonFacesRemoved": occupancy["internal_common_faces_removed"],
            "exposedFaces": occupancy["exposed_unit_faces"],
        },
        "envError": result["environment_scalar_smoothing_bound_m"],
        "linkError": result["maximum_link_smoothing_bound_m"],
        "clearanceScale": args.clearance_scale,
    }
    encoded = json.dumps(scene, separators=(",", ":"), ensure_ascii=False).replace(
        "</", "<\\/"
    )
    output = args.output or args.result.with_name(args.result.stem + "_3d.html")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(HTML.replace("__SCENE_DATA__", encoded))
    print(f"Wrote {output} ({len(points)} points, {len(spheres)} spheres)")


if __name__ == "__main__":
    main()
