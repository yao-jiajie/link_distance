#!/usr/bin/env python3
"""Render a link_distance_node JSON result and its XYZ cloud as an SVG."""

import argparse
import html
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
        point = tuple(float(value) for value in fields)
        if not all(math.isfinite(value) for value in point):
            raise ValueError(f"{path}:{line_number}: nonfinite coordinate")
        points.append(point)
    if not points:
        raise ValueError(f"{path}: empty point cloud")
    return points


def dot(lhs, rhs):
    return sum(a * b for a, b in zip(lhs, rhs))


def view_basis(azimuth_degrees, elevation_degrees):
    azimuth = math.radians(azimuth_degrees)
    elevation = math.radians(elevation_degrees)
    horizontal = (-math.sin(azimuth), math.cos(azimuth), 0.0)
    vertical = (
        -math.sin(elevation) * math.cos(azimuth),
        -math.sin(elevation) * math.sin(azimuth),
        math.cos(elevation),
    )
    depth = (
        math.cos(elevation) * math.cos(azimuth),
        math.cos(elevation) * math.sin(azimuth),
        math.sin(elevation),
    )
    return horizontal, vertical, depth


def interpolate_color(left, right, amount):
    amount = max(0.0, min(1.0, amount))
    rgb = tuple(round(a + amount * (b - a)) for a, b in zip(left, right))
    return "#{:02x}{:02x}{:02x}".format(*rgb)


def clearance_color(clearance, scale):
    if clearance <= 0.0:
        amount = min(1.0, -clearance / scale)
        return interpolate_color((244, 162, 97), (214, 40, 40), amount)
    amount = min(1.0, clearance / scale)
    return interpolate_color((244, 162, 97), (42, 157, 143), amount)


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
        matrix = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    if len(matrix) != 16 or not all(math.isfinite(float(value)) for value in matrix):
        raise ValueError("invalid field_from_base transform")
    matrix = [float(value) for value in matrix]
    translation = (matrix[3], matrix[7], matrix[11])

    def transform(point):
        delta = tuple(float(point[index]) - translation[index] for index in range(3))
        # field_T_base stores R; its rigid inverse uses R transpose.
        return tuple(
            sum(matrix[row * 4 + column] * delta[row] for row in range(3))
            for column in range(3)
        )

    return transform


def main():
    parser = argparse.ArgumentParser(
        description="Visualize point cloud and spherized robot link clearances"
    )
    parser.add_argument("result", type=Path, help="link_distance_node JSON")
    parser.add_argument("output", nargs="?", type=Path, help="output SVG")
    parser.add_argument("--cloud", type=Path, help="override XYZ path")
    parser.add_argument("--width", type=int, default=1200)
    parser.add_argument("--height", type=int, default=850)
    parser.add_argument("--azimuth", type=float, default=45.0)
    parser.add_argument("--elevation", type=float, default=25.0)
    parser.add_argument("--max-points", type=int, default=12000)
    parser.add_argument(
        "--clearance-scale",
        type=float,
        default=0.1,
        help="meters mapped from orange to green/red saturation",
    )
    args = parser.parse_args()
    if args.width < 200 or args.height < 200 or args.max_points <= 0:
        parser.error("width/height must be >=200 and max-points must be positive")
    if not math.isfinite(args.clearance_scale) or args.clearance_scale <= 0:
        parser.error("clearance-scale must be finite and positive")
    if not args.result.is_file():
        parser.error(
            f"result JSON does not exist: {args.result}; run link_distance_node first"
        )

    result = json.loads(args.result.read_text())
    spheres = result.get("spheres", [])
    if not spheres:
        raise ValueError("result has no sphere diagnostics; rerun link_distance_node")
    cloud_path = resolve_cloud(args.result, result, args.cloud)
    field_to_base = field_to_base_transform(result)
    points = [field_to_base(point) for point in load_xyz(cloud_path)]
    if len(points) > args.max_points:
        step = len(points) / args.max_points
        points = [points[int(index * step)] for index in range(args.max_points)]

    horizontal, vertical, depth_axis = view_basis(args.azimuth, args.elevation)

    def project(point):
        return dot(point, horizontal), dot(point, vertical), dot(point, depth_axis)

    projected_points = [project(point) for point in points]
    occupancy = result.get("voxel_occupancy")
    if not occupancy:
        raise ValueError("result has no voxel occupancy; rerun link_distance_node")
    projected_voxel_patches = []
    for patch in occupancy.get("boundary_patches", []):
        corners = patch.get("corners_field_m", [])
        if len(corners) != 4 or any(len(corner) != 3 for corner in corners):
            raise ValueError("invalid merged voxel boundary patch")
        projected_voxel_patches.append(
            ([project(field_to_base(corner)) for corner in corners], patch)
        )
    projected_spheres = []
    for sphere in spheres:
        center = field_to_base(sphere["center_field_m"])
        if len(center) != 3 or not all(math.isfinite(value) for value in center):
            raise ValueError("invalid sphere center")
        radius = float(sphere["radius_m"])
        if not math.isfinite(radius) or radius <= 0:
            raise ValueError("invalid sphere radius")
        projected_spheres.append((project(center), radius, sphere))

    u_values = [point[0] for point in projected_points]
    v_values = [point[1] for point in projected_points]
    for corners, _ in projected_voxel_patches:
        u_values.extend(point[0] for point in corners)
        v_values.extend(point[1] for point in corners)
    for (u, v, _), radius, _ in projected_spheres:
        u_values.extend((u - radius, u + radius))
        v_values.extend((v - radius, v + radius))
    u_min, u_max = min(u_values), max(u_values)
    v_min, v_max = min(v_values), max(v_values)
    if u_max == u_min:
        u_min -= 0.5
        u_max += 0.5
    if v_max == v_min:
        v_min -= 0.5
        v_max += 0.5

    margin_left, margin_right, margin_top, margin_bottom = 55, 245, 70, 55
    plot_width = args.width - margin_left - margin_right
    plot_height = args.height - margin_top - margin_bottom
    scale = min(plot_width / (u_max - u_min), plot_height / (v_max - v_min))
    used_width = scale * (u_max - u_min)
    used_height = scale * (v_max - v_min)
    offset_x = margin_left + (plot_width - used_width) / 2
    offset_y = margin_top + (plot_height - used_height) / 2

    def screen(u, v):
        return offset_x + (u - u_min) * scale, offset_y + (v_max - v) * scale

    output_path = args.output or args.result.with_suffix(".svg")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    nearest = {link["nearest_sphere"]: link for link in result["links"]}
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{args.width}" height="{args.height}" viewBox="0 0 {args.width} {args.height}">',
        '<rect width="100%" height="100%" fill="#f7f8fa"/>',
        '<style>text{font-family:system-ui,sans-serif}.label{font-size:12px;paint-order:stroke;stroke:white;stroke-width:3px;stroke-linejoin:round}.meta{font-size:13px;fill:#344054}</style>',
        f'<text x="{margin_left}" y="32" font-size="20" font-weight="600" fill="#172b4d">{html.escape(result.get("robot", "robot"))} link clearance proxy</text>',
        f'<text x="{margin_left}" y="53" class="meta">h={result["voxel_size_m"] * 1000:g} mm, points={result.get("point_count", len(points))}, occupied voxels={occupancy["occupied_voxels"]}, merged patches={occupancy["merged_boundary_patches"]}, q={html.escape(str(result.get("configuration", [])))}</text>',
        '<g id="merged-voxel-boundary">',
    ]
    for corners, patch in sorted(
        projected_voxel_patches,
        key=lambda item: sum(point[2] for point in item[0]) / len(item[0]),
    ):
        coordinates = " ".join(
            f"{screen(point[0], point[1])[0]:.3f},{screen(point[0], point[1])[1]:.3f}"
            for point in corners
        )
        colors = ("#4e90dc", "#57a4e0", "#66b5d5")
        lines.append(
            f'<polygon points="{coordinates}" fill="{colors[patch["axis"]]}" '
            'fill-opacity="0.20" stroke="#4a86bd" stroke-opacity="0.62" '
            'stroke-width="0.7" stroke-linejoin="round"/>'
        )
    lines.extend([
        '</g>',
        '<g id="point-cloud" fill="#5f6b7a" fill-opacity="0.32">',
    ])
    for u, v, _ in sorted(projected_points, key=lambda value: value[2]):
        x, y = screen(u, v)
        lines.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="1.15"/>')
    lines.append("</g>")

    lines.append('<g id="robot-spheres">')
    for (u, v, depth), radius, sphere in sorted(
        projected_spheres, key=lambda value: value[0][2]
    ):
        x, y = screen(u, v)
        clearance = float(sphere["clearance_proxy_m"])
        color = clearance_color(clearance, args.clearance_scale)
        is_nearest = sphere["id"] in nearest
        stroke_width = 3.0 if is_nearest else 1.2
        lines.append(
            f'<circle class="robot-sphere" data-sphere="{sphere["id"]}" '
            f'data-link="{html.escape(sphere["link"])}" cx="{x:.3f}" cy="{y:.3f}" '
            f'r="{radius * scale:.3f}" fill="{color}" fill-opacity="0.42" '
            f'stroke="{color}" stroke-width="{stroke_width}"/>'
        )
        if is_nearest:
            link = nearest[sphere["id"]]
            label = f'{link["name"]}: {link["distance_proxy_m"] * 1000:.1f} mm'
            lines.append(
                f'<text class="label" x="{x + 6:.3f}" y="{y - radius * scale - 5:.3f}" fill="#172b4d">{html.escape(label)}</text>'
            )
    lines.append("</g>")

    legend_x = args.width - margin_right + 28
    lines.extend(
        [
            f'<g id="legend" transform="translate({legend_x},90)">',
            '<text x="0" y="0" font-size="15" font-weight="600" fill="#172b4d">Sphere clearance</text>',
            '<rect x="2" y="18" width="16" height="16" fill="#4e90dc" fill-opacity="0.3" stroke="#4a86bd"/><text x="28" y="31" class="meta">merged voxel boundary</text>',
            '<circle cx="10" cy="52" r="8" fill="#d62828" fill-opacity="0.6"/><text x="28" y="57" class="meta">collision (≤ 0)</text>',
            '<circle cx="10" cy="80" r="8" fill="#f4a261" fill-opacity="0.6"/><text x="28" y="85" class="meta">near zero</text>',
            '<circle cx="10" cy="108" r="8" fill="#2a9d8f" fill-opacity="0.6"/><text x="28" y="113" class="meta">positive clearance</text>',
            '<circle cx="10" cy="142" r="8" fill="none" stroke="#172b4d" stroke-width="3"/><text x="28" y="147" class="meta">link nearest sphere</text>',
            f'<text x="0" y="181" class="meta">common faces removed: {occupancy["internal_common_faces_removed"]}</text>',
            f'<text x="0" y="203" class="meta">exposed faces: {occupancy["exposed_unit_faces"]} → {occupancy["merged_boundary_patches"]} patches</text>',
            f'<text x="0" y="235" class="meta">env ε ≤ {result["environment_scalar_smoothing_bound_m"] * 1000:.3g} mm</text>',
            f'<text x="0" y="257" class="meta">link ε ≤ {result["maximum_link_smoothing_bound_m"] * 1000:.3g} mm</text>',
            '<text x="0" y="295" class="meta">Orthographic projection</text>',
            f'<text x="0" y="317" class="meta">az={args.azimuth:g}°, el={args.elevation:g}°</text>',
            "</g>",
            "</svg>",
        ]
    )
    output_path.write_text("\n".join(lines) + "\n")
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
