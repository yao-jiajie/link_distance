#!/usr/bin/env python3
"""Render point cloud -> voxels -> hard zero set -> LSE zero set as one SVG."""

import argparse
import html
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np


REPO = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = REPO / "data/c_bracket_surface.xyz"
DEFAULT_ARTIFACTS = REPO / "output/c_bracket_visualization_50mm"
DEFAULT_OUTPUT = REPO / "docs/assets/point_cloud_pipeline_demo.svg"


def load_json(path):
    return json.loads(path.read_text())


def load_xyz(path):
    points = []
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        fields = line.split()
        if not fields or fields[0].startswith("#"):
            continue
        if len(fields) < 3:
            raise ValueError(f"{path}:{line_number}: expected at least three columns")
        point = tuple(float(value) for value in fields[:3])
        if not all(math.isfinite(value) for value in point):
            raise ValueError(f"{path}:{line_number}: nonfinite coordinate")
        points.append(point)
    if not points:
        raise ValueError(f"{path}: empty point cloud")
    return np.asarray(points, dtype=np.float64)


def evaluate_lse(tree, points, beta):
    """Evaluate the exported ordered DAG and return positive-free d_beta."""
    points = np.asarray(points, dtype=np.float64)
    uses = Counter(child for node in tree["nodes"] for child in node["children"])
    values = {}
    with np.errstate(under="ignore"):
        for node in tree["nodes"]:
            node_id = node["id"]
            if node["type"] == "leaf":
                leaf = tree["leaves"][node["leaf_id"]]
                normal = np.asarray(leaf["normal"], dtype=np.float64)
                values[node_id] = points @ normal - leaf["offset"]
                continue

            children = node["children"]
            arguments = [values[child] for child in children]
            maximum = node["type"] == "max"
            anchor = arguments[0].copy()
            for argument in arguments[1:]:
                (np.maximum if maximum else np.minimum)(anchor, argument, out=anchor)
            sign = 1.0 if maximum else -1.0
            total = np.zeros_like(anchor)
            for argument in arguments:
                total += np.exp(sign * beta * (argument - anchor))
            if not maximum:
                total /= len(arguments)
            values[node_id] = anchor + sign * np.log(total) / beta
            for child in children:
                uses[child] -= 1
                if uses[child] == 0 and child != tree["root_id"]:
                    del values[child]
    return -values[tree["root_id"]]


def evaluate_lse_with_gradient(tree, points, beta):
    """Evaluate positive-free d_beta and its world-coordinate analytic gradient."""
    points = np.asarray(points, dtype=np.float64)
    uses = Counter(child for node in tree["nodes"] for child in node["children"])
    values = {}
    gradients = {}
    with np.errstate(under="ignore"):
        for node in tree["nodes"]:
            node_id = node["id"]
            if node["type"] == "leaf":
                leaf = tree["leaves"][node["leaf_id"]]
                normal = np.asarray(leaf["normal"], dtype=np.float64)
                values[node_id] = points @ normal - leaf["offset"]
                gradients[node_id] = np.broadcast_to(normal, (len(points), 3))
                continue

            children = node["children"]
            arguments = [values[child] for child in children]
            maximum = node["type"] == "max"
            anchor = arguments[0].copy()
            for argument in arguments[1:]:
                (np.maximum if maximum else np.minimum)(anchor, argument, out=anchor)
            sign = 1.0 if maximum else -1.0
            weights = [np.exp(sign * beta * (argument - anchor)) for argument in arguments]
            total = np.zeros_like(anchor)
            gradient = np.zeros((len(points), 3), dtype=np.float64)
            for weight, child in zip(weights, children):
                total += weight
                gradient += weight[:, None] * gradients[child]
            normalizer = total if maximum else total / len(arguments)
            values[node_id] = anchor + sign * np.log(normalizer) / beta
            gradients[node_id] = gradient / total[:, None]
            for child in children:
                uses[child] -= 1
                if uses[child] == 0 and child != tree["root_id"]:
                    del values[child]
                    del gradients[child]
    root = tree["root_id"]
    return -values[root], -gradients[root]


def grid_zero_points(tree, beta, lower, upper, resolution, bisections):
    axes = [np.linspace(lower[a], upper[a], resolution) for a in range(3)]
    x, y, z = np.meshgrid(axes[0], axes[1], axes[2], indexing="ij")
    points = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
    values = evaluate_lse(tree, points, beta).reshape((resolution,) * 3)
    roots = []
    for axis in range(3):
        lo_slice = [slice(None)] * 3
        hi_slice = [slice(None)] * 3
        lo_slice[axis] = slice(0, -1)
        hi_slice[axis] = slice(1, None)
        lo_values = values[tuple(lo_slice)]
        hi_values = values[tuple(hi_slice)]
        crossing = np.signbit(lo_values) != np.signbit(hi_values)
        indices = np.argwhere(crossing)
        if not len(indices):
            continue
        lo = np.column_stack(tuple(axes[a][indices[:, a]] for a in range(3)))
        hi_indices = indices.copy()
        hi_indices[:, axis] += 1
        hi = np.column_stack(tuple(axes[a][hi_indices[:, a]] for a in range(3)))
        value_lo = lo_values[crossing].copy()
        for _ in range(bisections):
            middle = (lo + hi) * 0.5
            value_middle = evaluate_lse(tree, middle, beta)
            same_side = np.signbit(value_middle) == np.signbit(value_lo)
            lo[same_side] = middle[same_side]
            value_lo[same_side] = value_middle[same_side]
            hi[~same_side] = middle[~same_side]
        roots.append((lo + hi) * 0.5)
    if not roots:
        raise RuntimeError("No LSE zero crossing found in the sampling domain")
    result = np.concatenate(roots)
    residual = float(np.max(np.abs(evaluate_lse(tree, result, beta))))
    return result, residual


def exposed_voxel_faces(keys, h):
    occupied = {tuple(key) for key in keys}
    faces = []
    for key in occupied:
        for axis in range(3):
            for sign in (-1, 1):
                neighbour = list(key)
                neighbour[axis] += sign
                if tuple(neighbour) in occupied:
                    continue
                low = np.asarray(key, dtype=np.float64) * h
                high = low + h
                others = [a for a in range(3) if a != axis]
                vertices = []
                for first, second in ((0, 0), (1, 0), (1, 1), (0, 1)):
                    point = low.copy()
                    point[axis] = high[axis] if sign > 0 else low[axis]
                    point[others[0]] = high[others[0]] if first else low[others[0]]
                    point[others[1]] = high[others[1]] if second else low[others[1]]
                    vertices.append(point)
                normal = np.zeros(3)
                normal[axis] = sign
                faces.append((np.asarray(vertices), normal))
    return faces


def trace_cell_union(cells):
    """Return simplified, directed boundary loops for a union of grid cells."""
    edges = set()
    for u, v in cells:
        if (u, v - 1) not in cells:
            edges.add(((u, v), (u + 1, v)))
        if (u + 1, v) not in cells:
            edges.add(((u + 1, v), (u + 1, v + 1)))
        if (u, v + 1) not in cells:
            edges.add(((u + 1, v + 1), (u, v + 1)))
        if (u - 1, v) not in cells:
            edges.add(((u, v + 1), (u, v)))

    outgoing = {}
    for begin, end in edges:
        outgoing.setdefault(begin, []).append(end)
    direction_id = {(1, 0): 0, (0, 1): 1, (-1, 0): 2, (0, -1): 3}
    turn_rank = {1: 0, 0: 1, 3: 2, 2: 3}  # left, straight, right, reverse
    unused = set(edges)
    loops = []
    while unused:
        begin, current = min(unused)
        unused.remove((begin, current))
        loop = [begin]
        previous = begin
        while current != begin:
            loop.append(current)
            candidates = [end for end in outgoing.get(current, ()) if (current, end) in unused]
            if not candidates:
                raise ValueError("Open boundary while tracing a coplanar patch union")
            incoming = direction_id[(current[0] - previous[0], current[1] - previous[1])]

            def candidate_rank(end):
                outgoing_id = direction_id[(end[0] - current[0], end[1] - current[1])]
                return turn_rank[(outgoing_id - incoming) % 4]

            following = min(candidates, key=candidate_rank)
            unused.remove((current, following))
            previous, current = current, following

        simplified = []
        for index, point in enumerate(loop):
            previous = loop[index - 1]
            following = loop[(index + 1) % len(loop)]
            before = (point[0] - previous[0], point[1] - previous[1])
            after = (following[0] - point[0], following[1] - point[1])
            if before[0] * after[1] != before[1] * after[0]:
                simplified.append(point)
        if len(simplified) < 4:
            raise ValueError("Degenerate boundary loop in a coplanar patch union")
        loops.append(simplified)

    twice_area = sum(sum(loop[i][0] * loop[(i + 1) % len(loop)][1]
                         - loop[(i + 1) % len(loop)][0] * loop[i][1]
                         for i in range(len(loop))) for loop in loops)
    if twice_area != 2 * len(cells):
        raise ValueError("Coplanar patch union contour does not preserve its cell area")
    return loops


def hard_patch_unions(patches, h):
    """Group coplanar patches and extract only the boundary of each planar union."""
    groups = {}
    for patch in patches:
        key = (int(patch["axis"]), int(patch["sign"]), int(patch["plane_index"]))
        group = groups.setdefault(key, {"faces": [], "cells": set()})
        group["faces"].append(np.asarray(patch["corners_m"], dtype=np.float64))
        u0, u1, v0, v1 = (int(value) for value in patch["uv_bounds"])
        group["cells"].update((u, v) for v in range(v0, v1) for u in range(u0, u1))

    result = []
    for (axis, sign, plane), group in groups.items():
        cells = group["cells"]
        def point(u, v):
            index = np.zeros(3, dtype=np.float64)
            index[axis] = plane
            index[(axis + 1) % 3] = u
            index[(axis + 2) % 3] = v
            return index * h

        loops = [[point(u, v) for u, v in loop] for loop in trace_cell_union(cells)]
        normal = np.zeros(3, dtype=np.float64)
        normal[axis] = sign
        result.append({"faces": group["faces"], "normal": normal, "loops": loops})
    return result


def rgb(color):
    return "#" + "".join(f"{max(0, min(255, round(value))):02x}" for value in color)


def shade(base, amount):
    amount = max(-0.35, min(0.35, amount))
    if amount >= 0:
        return rgb(tuple(value + (255 - value) * amount for value in base))
    return rgb(tuple(value * (1 + amount) for value in base))


def svg_polygon(points, fill, stroke, stroke_width, opacity=1.0):
    coordinates = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
    return (f'<polygon points="{coordinates}" fill="{fill}" stroke="{stroke}" '
            f'stroke-width="{stroke_width}" stroke-linejoin="round" opacity="{opacity}"/>')


def svg_line(begin, end, stroke, stroke_width, opacity=1.0, dash=None):
    dash_attribute = f' stroke-dasharray="{dash}"' if dash else ''
    return (f'<line x1="{begin[0]:.2f}" y1="{begin[1]:.2f}" '
            f'x2="{end[0]:.2f}" y2="{end[1]:.2f}" stroke="{stroke}" '
            f'stroke-width="{stroke_width}" stroke-linecap="round" opacity="{opacity}"'
            f'{dash_attribute}/>')


def render(args):
    input_points = load_xyz(args.input)
    occupancy = load_json(args.artifacts / "occupancy_reference.json")
    boundary = load_json(args.artifacts / "boundary_patches.json")
    tree = load_json(args.artifacts / "free_direct_tree_raw.json")
    field = load_json(args.artifacts / "smooth_distance.json")
    h = float(occupancy["voxel_size"])
    keys = occupancy["occupied_keys"]
    beta = float(field["beta"])
    error = float(field["scalar_error_bound_m"])
    case_name = args.input.stem
    case_labels = {
        "realistic_sparse_camera": "相机样例",
        "u_surface": "U 型非凸样例",
        "two_objects_surface": "两个不连通物体样例",
        "disconnected_irregular_surface": "两个不连通不规则物体样例",
        "c_bracket_surface": "曲面＋直面 C 型非凸样例",
    }
    case_label = case_labels.get(case_name, case_name)
    sampling_label = "相机可见表面采样" if case_name == "realistic_sparse_camera" else "输入表面点采样"

    key_array = np.asarray(keys, dtype=np.float64)
    # The LSE zero set may move into free space. Include the full scalar error
    # budget as sampling headroom instead of assuming a sub-voxel displacement.
    padding = max(0.4 * h, 1.1 * error)
    lower = np.min(key_array, axis=0) * h - padding
    upper = (np.max(key_array, axis=0) + 1) * h + padding
    lse_points, residual = grid_zero_points(
        tree, beta, lower, upper, args.resolution, args.bisections)
    voxel_faces = exposed_voxel_faces(keys, h)
    patch_unions = hard_patch_unions(boundary["patches"], h)

    view_by_case = {
        # Look almost along -y so the empty x slab between the two components
        # remains visible instead of disappearing under perspective overlap.
        "two_objects_surface": (0.0, -1.70, 0.72),
        "disconnected_irregular_surface": (0.45, -1.70, 1.05),
    }
    view = np.asarray(
        view_by_case.get(case_name, (1.35, -1.70, 1.05)), dtype=np.float64)
    view /= np.linalg.norm(view)
    right = np.asarray((-view[1], view[0], 0.0))
    right /= np.linalg.norm(right)
    up = np.cross(view, right)
    center = (lower + upper) * 0.5

    def raw_projection(points):
        relative = np.asarray(points) - center
        return relative @ right, -(relative @ up), relative @ view

    corners = np.asarray([[x, y, z] for x in (lower[0], upper[0])
                          for y in (lower[1], upper[1])
                          for z in (lower[2], upper[2])])
    corner_u, corner_v, _ = raw_projection(corners)
    panel_width, panel_height = 350, 395
    draw_width, draw_height = 304, 286
    scale = min(draw_width / np.ptp(corner_u), draw_height / np.ptp(corner_v))
    panel_x = [20, 395, 770, 1145]
    draw_center_y = 292

    def project(points, index):
        u, v, depth = raw_projection(points)
        return (panel_x[index] + panel_width * 0.5 + scale * u,
                draw_center_y + scale * v, depth)

    width, height = 1515, 555
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
           f'viewBox="0 0 {width} {height}">',
           '<rect width="100%" height="100%" fill="#f4f7fb"/>',
           '<text x="28" y="36" font-family="sans-serif" font-size="24" font-weight="700" fill="#152033">'
           '点云到隐式函数的几何演变</text>',
           '<text x="28" y="62" font-family="sans-serif" font-size="13" fill="#526277">'
           f'{html.escape(case_label)} · 同一视角与尺度 · 几何单位：m</text>']

    titles = ("1  原始点云", "2  占据体素", "3  硬零面", "4  LSE 零面对比")
    subtitles = (sampling_label, "仅包含点所在的体素", "q(x) = 0",
                 "绿色：LSE · 橙虚线：硬零面")
    for index, x0 in enumerate(panel_x):
        svg.append(f'<rect x="{x0}" y="78" width="{panel_width}" height="{panel_height}" rx="12" '
                   'fill="white" stroke="#d9e1eb"/>')
        svg.append(f'<text x="{x0 + 18}" y="110" font-family="sans-serif" font-size="18" '
                   f'font-weight="700" fill="#172033">{titles[index]}</text>')
        svg.append(f'<text x="{x0 + 18}" y="132" font-family="sans-serif" font-size="12" '
                   f'fill="#65758b">{subtitles[index]}</text>')
        svg.append(f'<clipPath id="panel-{index}"><rect x="{x0 + 8}" y="140" '
                   f'width="{panel_width - 16}" height="292" rx="6"/></clipPath>')
        if index < 3:
            svg.append(f'<text x="{x0 + panel_width + 11}" y="282" font-family="sans-serif" '
                       'font-size="25" fill="#9aa8ba">→</text>')

    # Raw point cloud.
    sx, sy, depth = project(input_points, 0)
    order = np.argsort(depth)
    depth_range = max(float(np.ptp(depth)), 1e-12)
    svg.append('<g clip-path="url(#panel-0)">')
    for index in order:
        t = (depth[index] - np.min(depth)) / depth_range
        color = shade((43, 108, 176), 0.28 - 0.30 * t)
        svg.append(f'<circle cx="{sx[index]:.2f}" cy="{sy[index]:.2f}" r="1.35" '
                   f'fill="{color}" opacity="0.78"/>')
    svg.append('</g>')

    def draw_faces(faces, panel, base, stroke, line_width):
        visible = []
        for vertices, normal in faces:
            facing = float(normal @ view)
            if facing <= 0.025:
                continue
            px, py, face_depth = project(vertices, panel)
            visible.append((float(np.mean(face_depth)), facing, np.column_stack((px, py))))
        svg.append(f'<g clip-path="url(#panel-{panel})">')
        for _, facing, points in sorted(visible, key=lambda item: item[0]):
            fill = shade(base, 0.55 * facing - 0.28) if stroke is None \
                else shade(base, 0.20 * facing - 0.07)
            # A same-color stroke overlaps adjacent rasterized faces by less
            # than one screen pixel and hides antialiasing cracks without
            # exposing the unit-face or patch partition.
            outline = fill if stroke is None else stroke
            svg.append(svg_polygon(points, fill, outline, line_width, 0.96))
        svg.append('</g>')
        return len(visible)

    def visible_union_outline_pieces(groups, patches):
        """Find contour portions that are nearest along the viewing ray."""
        occluders = []
        for patch in patches:
            normal = np.asarray(patch["normal"], dtype=np.float64)
            if float(normal @ view) <= 0.025:
                continue
            vertices = np.asarray(patch["corners_m"], dtype=np.float64)
            axis = int(patch["axis"])
            occluders.append((axis, vertices[0, axis], np.min(vertices, axis=0),
                              np.max(vertices, axis=0)))

        def point_visible(point):
            relative = point - center
            u = float(relative @ right)
            v = -float(relative @ up)
            depth = float(relative @ view)
            ray_base = center + u * right - v * up
            nearest = -math.inf
            for axis, plane, low, high in occluders:
                ray_depth = (plane - ray_base[axis]) / view[axis]
                intersection = ray_base + ray_depth * view
                others = [coordinate for coordinate in range(3) if coordinate != axis]
                if all(low[coordinate] - 1e-10 <= intersection[coordinate]
                       <= high[coordinate] + 1e-10 for coordinate in others):
                    nearest = max(nearest, ray_depth)
            return depth >= nearest - 1e-8

        pieces = []
        source_segments = 0
        visible_groups = 0
        for group in groups:
            if float(group["normal"] @ view) <= 0.025:
                continue
            visible_groups += 1
            for loop in group["loops"]:
                for begin, end in zip(loop, np.roll(loop, -1, axis=0)):
                    source_segments += 1
                    screen_u, screen_v, _ = raw_projection(np.asarray((begin, end)))
                    screen_begin = scale * np.asarray((screen_u[0], screen_v[0]))
                    screen_end = scale * np.asarray((screen_u[1], screen_v[1]))
                    bins = max(1, int(math.ceil(np.linalg.norm(screen_end - screen_begin) / 0.75)))
                    run_begin = None
                    for index in range(bins):
                        alpha = (index + 0.5) / bins
                        visible = point_visible((1.0 - alpha) * begin + alpha * end)
                        if visible and run_begin is None:
                            run_begin = index
                        if run_begin is not None and (not visible or index + 1 == bins):
                            run_end = index if not visible else index + 1
                            first_alpha = run_begin / bins
                            last_alpha = run_end / bins
                            first = (1.0 - first_alpha) * begin + first_alpha * end
                            last = (1.0 - last_alpha) * begin + last_alpha * end
                            pieces.append((first, last))
                            run_begin = None
        return visible_groups, source_segments, pieces

    def draw_outline_pieces(pieces, panel, stroke, line_width, opacity=0.96, dash=None):
        svg.append(f'<g clip-path="url(#panel-{panel})">')
        for begin, end in pieces:
            px, py, _ = project(np.asarray((begin, end)), panel)
            svg.append(svg_line((px[0], py[0]), (px[1], py[1]), stroke,
                                line_width, opacity, dash))
        svg.append('</g>')

    visible_voxels = draw_faces(voxel_faces, 1, (72, 145, 211), "#e8f2fb", 0.42)
    visible_hard_faces = draw_faces(voxel_faces, 2, (239, 136, 61), None, 0.72)
    visible_patch_groups, hard_outline_segments, hard_outline_pieces = \
        visible_union_outline_pieces(patch_unions, boundary["patches"])
    draw_outline_pieces(hard_outline_pieces, 2, "#a65324", 0.48)

    # Keep only the nearest LSE root in each small screen bin. Analytic-gradient
    # lighting turns the samples into a readable surface instead of a uniform point cloud.
    sx, sy, depth = project(lse_points, 3)
    nearest = {}
    for index, (x, y, d) in enumerate(zip(sx, sy, depth)):
        bucket = (round(x / 1.15), round(y / 1.15))
        if bucket not in nearest or d > depth[nearest[bucket]]:
            nearest[bucket] = index
    visible_roots = np.asarray(list(nearest.values()), dtype=np.int64)
    visible_values, visible_gradients = evaluate_lse_with_gradient(
        tree, lse_points[visible_roots], beta)
    gradient_norms = np.linalg.norm(visible_gradients, axis=1)
    normals = visible_gradients / np.maximum(gradient_norms[:, None], 1e-12)
    light = np.asarray((0.25, -0.45, 1.0), dtype=np.float64)
    light /= np.linalg.norm(light)
    illumination = 0.30 + 0.70 * np.maximum(normals @ light, 0.0)
    order = visible_roots[np.argsort(depth[visible_roots])]
    root_to_visible = {root: index for index, root in enumerate(visible_roots)}
    svg.append('<g clip-path="url(#panel-3)">')
    for index in order:
        lit = illumination[root_to_visible[index]]
        base = np.asarray((35, 180, 125), dtype=np.float64)
        color = rgb(base * (0.58 + 0.55 * lit) + (255.0 - base) * (0.10 * lit))
        svg.append(f'<circle cx="{sx[index]:.2f}" cy="{sy[index]:.2f}" r="1.65" '
                   f'fill="{color}"/>')
    svg.append('</g>')
    draw_outline_pieces(hard_outline_pieces, 3, "#d66f2f", 0.72, 0.78, "3 2")

    footers = (
        f"{len(input_points)} 个 XYZ 点",
        f"{len(keys)} 个体素 · h = {h:g} m",
        f"{len(boundary['patches'])} 个矩形 patch · 精确边界",
        f"β = {beta:.1f} m⁻¹ · 误差界 {1000 * error:g} mm",
    )
    for index, footer in enumerate(footers):
        svg.append(f'<text x="{panel_x[index] + panel_width / 2}" y="453" text-anchor="middle" '
                   f'font-family="sans-serif" font-size="12" fill="#526277">{html.escape(footer)}</text>')

    svg.extend([
        '<line x1="28" y1="500" x2="1487" y2="500" stroke="#d9e1eb"/>',
        f'<text x="28" y="526" font-family="sans-serif" font-size="12" fill="#526277">'
        f'LSE 面：{args.resolution}³ 网格检测符号变化，边上二分 {args.bisections} 次；'
        f'{len(lse_points)} 个求根点，可见采样 {len(visible_roots)} 个，最大求值残差 '
        f'{residual:.2e} m；绿色按解析梯度着色。</text>',
        '<text x="1487" y="526" text-anchor="end" font-family="sans-serif" font-size="12" '
        'fill="#526277">橙色虚线为硬零面位置参照。</text>',
        f'<!-- visible voxel faces: {visible_voxels}; visible hard unit faces: {visible_hard_faces}; '
        f'visible coplanar groups: {visible_patch_groups}; hard outline segments: '
        f'{hard_outline_segments}; visible outline pieces: {len(hard_outline_pieces)}; '
        f'LSE gradient norm range: {float(np.min(gradient_norms)):.6g} to '
        f'{float(np.max(gradient_norms)):.6g} -->',
        '</svg>',
    ])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(svg) + "\n")
    return {
        "output": str(args.output), "raw_points": len(input_points), "occupied_voxels": len(keys),
        "hard_patches": len(boundary["patches"]), "lse_root_points": len(lse_points),
        "lse_visible_points": len(visible_roots), "lse_max_residual_m": residual,
        "hard_visible_coplanar_groups": visible_patch_groups,
        "hard_outline_segments": hard_outline_segments,
        "hard_visible_outline_pieces": len(hard_outline_pieces),
        "lse_visible_max_residual_m": float(np.max(np.abs(visible_values))),
        "lse_visible_gradient_norm_min": float(np.min(gradient_norms)),
        "lse_visible_gradient_norm_max": float(np.max(gradient_norms)),
        "resolution": args.resolution, "bisections": args.bisections,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--artifacts", type=Path, default=DEFAULT_ARTIFACTS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--resolution", type=int, default=64)
    parser.add_argument("--bisections", type=int, default=8)
    args = parser.parse_args()
    if args.resolution < 8 or not 0 <= args.bisections <= 20:
        parser.error("resolution must be >= 8 and bisections must be in [0, 20]")
    args.input = args.input.resolve()
    args.artifacts = args.artifacts.resolve()
    args.output = args.output.resolve()
    print(json.dumps(render(args), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
