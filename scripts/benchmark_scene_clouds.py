"""Deterministic meter-scale surface scenes for Phase 4; no sensor claims."""
import itertools
import math
import random


SCENES = {
    "cube": {"boxes": [((-.15, -.15, .50), (.15, .15, .80))]},
    "l": {"boxes": [((-.15, -.15, .50), (.15, -.06, .65)),
                     ((-.15, -.06, .50), (-.06, .15, .65))],
          "gap": ((-.03, -.03, .515), (.12, .12, .635))},
    "u": {"boxes": [((-.15, -.15, .50), (.15, -.06, .65)),
                     ((-.15, -.06, .50), (-.06, .15, .65)),
                     ((.06, -.06, .50), (.15, .15, .65))],
          "gap": ((-.04, -.02, .515), (.04, .12, .635))},
    "thin_wall": {"boxes": [((-.003, -.15, .50), (.003, .15, .80))]},
    "ellipsoid": {"center": (0., 0., .75), "radii": (.18, .13, .24)},
}


def inside_box(point, box):
    return all(lo <= p <= hi for p, lo, hi in zip(point, *box))


def inside_scene(point, scene):
    if "boxes" in scene:
        return any(inside_box(point, box) for box in scene["boxes"])
    return sum(((p-c)/r)**2 for p, c, r in zip(point, scene["center"], scene["radii"])) <= 1


def bounds(scene):
    if "boxes" in scene:
        return ([min(b[0][k] for b in scene["boxes"]) for k in range(3)],
                [max(b[1][k] for b in scene["boxes"]) for k in range(3)])
    return ([c-r for c, r in zip(scene["center"], scene["radii"])],
            [c+r for c, r in zip(scene["center"], scene["radii"])])


def box_volume(box):
    return math.prod(max(0., hi-lo) for lo, hi in zip(*box))


def solid_volume(scene):
    if "boxes" not in scene:
        return 4*math.pi*math.prod(scene["radii"])/3
    # Inclusion-exclusion also works if future scene boxes overlap.
    volume = 0.
    for n in range(1, len(scene["boxes"])+1):
        for boxes in itertools.combinations(scene["boxes"], n):
            intersection = ([max(b[0][k] for b in boxes) for k in range(3)],
                            [min(b[1][k] for b in boxes) for k in range(3)])
            volume += (-1)**(n+1)*box_volume(intersection)
    return volume


def sample_box(box, rng):
    return tuple(rng.uniform(lo, hi) for lo, hi in zip(*box))


def surface_points(scene, count, seed):
    rng = random.Random(seed)
    if "boxes" not in scene:
        points = []
        for _ in range(count):
            z = rng.uniform(-1, 1)
            phi = rng.uniform(0, 2*math.pi)
            radial = math.sqrt(max(0, 1-z*z))
            # Uniform directions scaled to ellipsoid: not uniform surface area.
            direction = (radial*math.cos(phi), radial*math.sin(phi), z)
            points.append(tuple(c+r*d for c, r, d in zip(scene["center"], scene["radii"], direction)))
        return points
    faces, weights = [], []
    for box in scene["boxes"]:
        for axis, side in itertools.product(range(3), (0, 1)):
            faces.append((box, axis, side))
            weights.append(math.prod(box[1][k]-box[0][k] for k in range(3) if k != axis))
    points = []
    while len(points) < count:
        box, axis, side = rng.choices(faces, weights=weights, k=1)[0]
        point = list(sample_box(box, rng))
        point[axis] = box[side][axis]
        outward = point.copy()
        outward[axis] += 1e-9 if side else -1e-9
        if inside_scene(outward, scene):
            continue  # hidden interface of the box union, not an object surface
        points.append(tuple(point))
    return points


def noisy_points(points, halfwidth, seed):
    rng = random.Random(seed)
    return [tuple(p+rng.uniform(-halfwidth, halfwidth) for p in point) for point in points]


def interior_points(scene, count, seed):
    rng = random.Random(seed)
    box = bounds(scene)
    points = []
    while len(points) < count:
        point = sample_box(box, rng)
        if inside_scene(point, scene):
            points.append(point)
    return points
