#!/usr/bin/env python3
"""Independently audit every collision sphere's movable-joint reference frame."""

import argparse
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path


def vector(text, default):
    values = [float(value) for value in (text or default).split()]
    if len(values) != 3 or not all(math.isfinite(value) for value in values):
        raise ValueError(f"invalid three-vector: {text}")
    return values


def identity():
    return [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]


def multiply(left, right):
    return [[sum(left[row][k] * right[k][column] for k in range(4))
             for column in range(4)] for row in range(4)]


def apply(transform, point):
    return [sum(transform[row][column] * point[column] for column in range(3))
            + transform[row][3] for row in range(3)]


def rotation_matrix(axis, angle):
    norm = math.sqrt(sum(value * value for value in axis))
    if not norm > 0.0:
        raise ValueError("zero joint axis")
    x, y, z = (value / norm for value in axis)
    c, s, one = math.cos(angle), math.sin(angle), 1.0 - math.cos(angle)
    return [
        [c + x * x * one, x * y * one - z * s, x * z * one + y * s],
        [y * x * one + z * s, c + y * y * one, y * z * one - x * s],
        [z * x * one - y * s, z * y * one + x * s, c + z * z * one],
    ]


def rigid(xyz, rotation):
    result = identity()
    for row in range(3):
        for column in range(3):
            result[row][column] = rotation[row][column]
        result[row][3] = xyz[row]
    return result


def origin(element):
    node = element.find("origin")
    if node is None:
        return identity()
    xyz = vector(node.get("xyz"), "0 0 0")
    roll, pitch, yaw = vector(node.get("rpy"), "0 0 0")
    rx = rotation_matrix([1, 0, 0], roll)
    ry = rotation_matrix([0, 1, 0], pitch)
    rz = rotation_matrix([0, 0, 1], yaw)
    rzyx = [[sum(rz[row][k] * ry[k][column] for k in range(3))
             for column in range(3)] for row in range(3)]
    rotation = [[sum(rzyx[row][k] * rx[k][column] for k in range(3))
                 for column in range(3)] for row in range(3)]
    return rigid(xyz, rotation)


def joint_motion(joint_type, axis, value):
    if joint_type in ("revolute", "continuous"):
        return rigid([0, 0, 0], rotation_matrix(axis, value))
    if joint_type == "prismatic":
        norm = math.sqrt(sum(component * component for component in axis))
        return rigid([value * component / norm for component in axis],
                     identity_rotation())
    return identity()


def identity_rotation():
    return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]


def error(left, right):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right)))


def audit(result_path, urdf_path, tolerance):
    result = json.loads(result_path.read_text())
    robot = ET.parse(urdf_path).getroot()
    link_elements = {link.get("name"): link for link in robot.findall("link")}
    if len(link_elements) != len(robot.findall("link")):
        raise ValueError("duplicate or unnamed URDF link")

    joints = []
    parent_by_child = {}
    movable_names = []
    for element in robot.findall("joint"):
        joint_type = element.get("type")
        parent = element.find("parent").get("link")
        child = element.find("child").get("link")
        item = {
            "name": element.get("name"), "type": joint_type,
            "parent": parent, "child": child, "origin": origin(element),
            "axis": vector(element.find("axis").get("xyz"), "1 0 0")
                    if element.find("axis") is not None else [1.0, 0.0, 0.0],
            "index": None,
        }
        if joint_type != "fixed":
            item["index"] = len(movable_names)
            movable_names.append(item["name"])
        if child in parent_by_child:
            raise ValueError(f"multiple parents for {child}")
        parent_by_child[child] = item
        joints.append(item)
    if movable_names != result["joint_order"]:
        raise ValueError("URDF movable-joint order differs from result joint_order")

    roots = set(link_elements) - set(parent_by_child)
    if len(roots) != 1:
        raise ValueError("URDF must have one root link")
    field = result["field_from_base"]
    if len(field) != 16:
        raise ValueError("invalid field_from_base")
    transforms = {roots.pop(): [field[row * 4:(row + 1) * 4] for row in range(4)]}
    pending = list(joints)
    while pending:
        next_pending = []
        progress = False
        for joint in pending:
            if joint["parent"] not in transforms:
                next_pending.append(joint)
                continue
            value = 0.0 if joint["index"] is None else result["configuration"][joint["index"]]
            transforms[joint["child"]] = multiply(
                multiply(transforms[joint["parent"]], joint["origin"]),
                joint_motion(joint["type"], joint["axis"], value))
            progress = True
        if not progress:
            raise ValueError("disconnected or cyclic URDF")
        pending = next_pending

    expected_spheres = []
    for link in robot.findall("link"):
        for collision in link.findall("collision"):
            expected_spheres.append((link.get("name"), apply(origin(collision), [0, 0, 0])))
    if len(expected_spheres) != len(result["spheres"]):
        raise ValueError("URDF/result sphere count mismatch")

    group_by_source = {}
    for group in result["links"]:
        index = group["joint_configuration_index"]
        if not 0 <= index < len(movable_names) or group["joint"] != movable_names[index]:
            raise ValueError(f"invalid result joint index for {group['name']}")
        for source in group["source_links"]:
            if source in group_by_source:
                raise ValueError(f"source link appears in two groups: {source}")
            group_by_source[source] = group

    maximum_urdf_error = maximum_group_error = maximum_field_error = 0.0
    grouped = 0
    fixed_descendant_spheres = 0
    group_counts = {group["name"]: 0 for group in result["links"]}
    for sphere_id, ((owner, center_owner), actual) in enumerate(
            zip(expected_spheres, result["spheres"])):
        if actual["id"] != sphere_id or actual["link"] != owner:
            raise ValueError(f"sphere order/owner mismatch at {sphere_id}")
        maximum_urdf_error = max(maximum_urdf_error,
                                 error(center_owner, actual["center_urdf_link_m"]))
        expected_field = apply(transforms[owner], center_owner)
        maximum_field_error = max(maximum_field_error,
                                  error(expected_field, actual["center_field_m"]))

        center_group = center_owner
        current = owner
        expected_joint = expected_group_link = None
        expected_index = None
        while current in parent_by_child:
            joint = parent_by_child[current]
            if joint["index"] is not None:
                expected_joint = joint["name"]
                expected_group_link = joint["child"]
                expected_index = joint["index"]
                break
            center_group = apply(joint["origin"], center_group)
            current = joint["parent"]
        if expected_joint is None:
            if any(actual[key] is not None for key in
                   ("joint_group", "group_link", "joint_configuration_index",
                    "center_group_link_m")):
                raise ValueError(f"static sphere {sphere_id} was assigned to a joint")
            if owner in group_by_source:
                raise ValueError(f"static source link {owner} appears in a distance group")
            continue

        grouped += 1
        if owner != expected_group_link:
            fixed_descendant_spheres += 1
        if (actual["joint_group"] != expected_joint or
                actual["group_link"] != expected_group_link or
                actual["joint_configuration_index"] != expected_index):
            raise ValueError(f"wrong joint frame for sphere {sphere_id} on {owner}")
        group = group_by_source.get(owner)
        if group is None or group["name"] != expected_group_link or group["joint"] != expected_joint:
            raise ValueError(f"wrong distance group for sphere {sphere_id} on {owner}")
        group_counts[group["name"]] += 1
        maximum_group_error = max(maximum_group_error,
                                  error(center_group, actual["center_group_link_m"]))

    for group in result["links"]:
        if group_counts[group["name"]] != group["sphere_count"]:
            raise ValueError(f"sphere count mismatch for {group['name']}")
    maximum = max(maximum_urdf_error, maximum_group_error, maximum_field_error)
    if maximum > tolerance:
        raise ValueError(f"sphere frame error {maximum} exceeds tolerance {tolerance}")
    return {
        "spheres": len(result["spheres"]),
        "grouped_spheres": grouped,
        "static_spheres": len(result["spheres"]) - grouped,
        "fixed_descendant_spheres": fixed_descendant_spheres,
        "joint_groups": len(result["links"]),
        "group_sphere_counts": group_counts,
        "max_urdf_link_center_error_m": maximum_urdf_error,
        "max_group_link_center_error_m": maximum_group_error,
        "max_field_center_error_m": maximum_field_error,
        "tolerance_m": tolerance,
        "passed": True,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Validate every sphere's URDF, joint-group, and field frames")
    parser.add_argument("result", type=Path, help="link_distance_node JSON")
    parser.add_argument("--urdf", type=Path, help="override URDF path")
    parser.add_argument("--tolerance", type=float, default=1e-10)
    args = parser.parse_args()
    result = json.loads(args.result.read_text())
    urdf = args.urdf or Path(result["urdf"])
    if not urdf.is_file():
        candidate = args.result.parent / urdf
        if not candidate.is_file():
            parser.error(f"URDF does not exist: {urdf}; pass --urdf")
        urdf = candidate
    if not math.isfinite(args.tolerance) or args.tolerance <= 0:
        parser.error("tolerance must be finite and positive")
    print(json.dumps(audit(args.result, urdf, args.tolerance), indent=2))


if __name__ == "__main__":
    main()
