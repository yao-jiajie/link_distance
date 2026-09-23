# Voxel input representation and construction optimization

The current representation is a tested baseline, not a proven optimum for
runtime, tree size, or geometric accuracy.

## Exact input to Alpha Wrap

For each occupied grid index `(i,j,k)` and size `s`, emit the eight positions
`s*(i+a,j+b,k+c)` for `a,b,c` in `{0,1}`. Deduplicate integer corner keys
globally and order them lexicographically. `preprocessPointCloud()` moves only
the `points` vector into its output. `AlphaWrapVolumeAdapter` calls
`Point_set_oracle::add_point_set(points)`.

Alpha Wrap receives point primitives in meters, with no occupied-box indices,
face connectivity, box dimensions, weights or radii. Voxel size controls
global parameter selection; the oracle itself sees only points. The exact
passed points are saved as `alpha_wrap_input_points.xyz`.

For the supplied camera cloud with `s=0.03 m`:

```text
1655 raw points -> 543 occupied boxes -> 4344 corner occurrences
  -> 1253 unique points -> Point_set_oracle
```

This camera input was already sparsified with 18 mm voxels in the simulator,
so it is not a dense input on which corner voxelization necessarily removes
most samples.

## Geometry semantics

Raw points belong to the occupied box union. The discrete corner set does
not encode those solid boxes in Alpha Wrap, so the resulting wrap still needs
original-point checks and offset recovery. Voxelization does not fit a smooth
surface; it fills occupied grid cells. Small movements near a grid boundary
can change occupancy.

For the 3 cm example, the distances from the 1253 corners to their nearest
raw sample have mean 22.08 mm, 95th percentile 39.63 mm, and maximum 48.58 mm.
These are unsigned distances to the discrete input, not outward error against
the true surface. A box corner can be almost `sqrt(3)*s` from its raw sample.
Voxel size therefore does not bound final geometry error; global alpha/offset
also affect the result.

## Alternative representation counts

Exposed faces exclude shared faces of adjacent occupied voxels, but include
faces towards enclosed empty space. Each exposed square below is split into
two triangles, without coplanar merging or geometric simplification.

| Size (m) | All unique corners | Exposed-face corners | Exposed squares | Triangles |
| --- | --- | --- | --- | --- |
| 0.01 | 7990 | 7990 | 8040 | 16080 |
| 0.02 | 2628 | 2626 | 3144 | 6288 |
| 0.03 | 1253 | 1250 | 1396 | 2792 |

Exposed-face corners save only three points at 3 cm. A triangle-soup oracle
would carry continuous boundary geometry, but this direct construction uses
2792 triangle primitives. Primitive counts alone do not predict runtime;
triangle wrapping and coplanar patch merging have not been benchmarked here.
A single center/farthest point is not a conservative substitute for a box
without an additional conservative geometric construction.

At matched alpha 0.09 m and offset 0.03 m, the earlier raw-point baseline has
473 tree nodes while corners produce 626. The choice should be evaluated
against a specified geometric error budget and full-pipeline time/tree size.

## Same-output implementation optimization

Occupied-cell and corner `std::set` builders are replaced with reserved
vectors and `sort`/`unique`. This avoids per-node allocations while preserving
coordinates, ordering, numerical indexing, and containment checks. Temporary
corner storage holds eight records per occupied voxel before compaction.

```bash
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
build/voxelization_benchmark data/realistic_sparse_camera.xyz
```

The benchmark compares the old ordered-set builder with production in the
same Release process. It alternates execution order, excludes three warm-ups,
and reports medians of 30 runs. Every run compares occupied cells, corner
keys, coordinates/order, and containment outcomes.

| Size (m) | Set construction (ms) | Vector construction (ms) | Set with containment (ms) | Vector with containment (ms) |
| --- | --- | --- | --- | --- |
| 0.01 | 2.29928 | 1.37529 | 2.50595 | 1.58187 |
| 0.02 | 1.39096 | 1.01952 | 1.59951 | 1.22886 |
| 0.03 | 0.76584 | 0.61912 | 0.94342 | 0.79684 |

At 3 cm construction decreases by about 19%. Compare these paired timings,
not absolute times from older runs under different scheduling/cache conditions.

The full regression in `output/voxel_stage_optimization` passes all checks and
has byte-identical `alpha_wrap_input_points.xyz` to the old 3 cm run, with
118 wrap faces, 161 tetrahedra, 101 clusters, 524 planes, and 626 tree nodes.
Its single core time was 12.12 ms, and 39.32 ms including validation; no
whole-pipeline speedup is claimed from that one measurement.
