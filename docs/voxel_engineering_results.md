# Conservative voxel pipeline: first engineering benchmark

The independent corner-voxelization module and its integration with both
Alpha Wrap entry points are implemented. The existing internal tetrahedral
extraction, adjacency, greedy convex aggregation, support-plane reduction,
and MIN–MAX construction algorithms are reused. Jet is available explicitly
for comparisons. `exact-nef` remains a separate legacy route.

## Reproduction and measurement scope

```bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
python3 scripts/benchmark_voxel_pipeline.py data/realistic_sparse_camera.xyz \
  output/another_voxel_benchmark --repeats 3 --validation-samples 1000 \
  --with-nef --timeout 60
```

Measured on the current machine using the Release build, CGAL 5.6.3, no TBB.
Input: 1655 simulated camera samples in meters. Each voxel/raw/Jet case runs
sequentially three times; reported times are column-wise medians. These
medians need not add up across stages. Every successful full run checks all
original samples in the wrap and final expression tree, independently of its
1000 random validation samples. `exact` plane reduction is used throughout.

The complete latest measurements are in
[`summary.csv`](../output/voxel_engineering_benchmark_v2/summary.csv) and
[`summary.json`](../output/voxel_engineering_benchmark_v2/summary.json).
Per-case directories include commands, logs, first-run geometry and validation
reports, and three repetitions in `trials.json`. The earlier
`voxel_engineering_benchmark` directory is an interrupted reporting attempt;
use `voxel_engineering_benchmark_v2` for the complete comparison.

## Voxel size sweep

| Voxel size (m) | Alpha (m) | Offset (m) | Occupied voxels | Unique corners | Wrap faces | Tetrahedra | Convex clusters | Planes before/after | Tree nodes | Core (ms) | With validation (ms) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.01 | 0.03 | 0.01 | 1562 | 7990 | 742 | 1462 | 930 | 4784 / 4784 | 5715 | 191.84 | 472.58 |
| 0.02 | 0.06 | 0.02 | 994 | 2628 | 224 | 346 | 214 | 1120 / 1120 | 1335 | 27.00 | 108.97 |
| 0.03 | 0.09 | 0.03 | 543 | 1253 | 118 | 161 | 101 | 524 / 524 | 626 | 11.20 | 39.36 |

All raw-point voxel/wrap/tree outside counts are zero. All mesh/tetra,
tetra/cluster, tetra-tree, and raw/simplified-tree sampled mismatch counts
are zero. No offset retry was needed for these voxel-factor defaults.
Construction times exclude the measured validation stages and file I/O;
the final column includes validation but still excludes initial point-file
reading and result-file writing. The script also records whole-process wall
time including I/O. Mandatory containment checks are part of validation, so
core time alone is not the latency of a validated collision representation.

Larger voxels in this sweep also select larger alpha/offset. The reduction
from 5715 to 626 nodes therefore cannot be attributed solely to voxelization.

## Matched requested alpha/offset comparisons

| Scale (m) | Preprocessing | Input points to Wrap | Tree nodes | Core median (ms) |
| --- | --- | --- | --- | --- |
| 0.01 | Voxel corners | 7990 | 5715 | 191.84 |
| 0.01 | None | 1655 | 6180 | 171.57 |
| 0.01 | Jet | 1655 | 6029 | 197.74 |
| 0.02 | Voxel corners | 2628 | 1335 | 27.00 |
| 0.02 | None | 1655 | 1442 | 33.39 |
| 0.02 | Jet | 1655 | 1433 | 65.11 |
| 0.03 | Voxel corners | 1253 | 626 | 11.20 |
| 0.03 | None | 1655 | 473 | 12.68 |
| 0.03 | Jet | 1655 | 675 | 49.30 |

At 1 cm the corner representation increases input size by 4.83x. At 3 cm it
reduces input size by 24.3%, but the tree is larger than the unsmoothed
same-scale baseline. Exact coplanar reduction did not remove planes on these
cases. The implemented conservative near-coplanar option has containment
regression coverage; its performance is not measured in the table above.

## Legacy Nef comparison on the same voxel wrap

| Voxel size (m) | Convex parts | Unique planes | Tree nodes | OFF-to-tree including validation (ms) | Result |
| --- | --- | --- | --- | --- | --- |
| 0.01 | 465 | 1516 | 4550 | 30961.35 | Passed |
| 0.02 | 85 | 368 | 991 | 4749.75 | Passed |
| 0.03 | — | — | — | — | Failed: `Facet normal cannot be normalized` |

Nef runs are single repetitions and begin at OFF input. These are not full
point-cloud-to-tree timings and should not be compared directly to alpha-tet
core medians. `shared_prefix_plus_nef_ms` in the summary is explicitly a stage
sum, excluding OFF writing and oracle setup. The 3 cm Nef failure is preserved
in its log; no failed result is reported as a successful collision tree. The
legacy `l_prism.off` two-part decomposition regression still passes.

## Validation coverage and limits

Tests cover negative coordinates, adjacent voxel corner sharing, exact and
neighboring floating-point grid boundaries, closed upper-boundary membership,
independent brute-force box containment, input permutation and duplication,
invalid sizes/coordinates/grid ranges, three parameter modes and conflicts,
voxel-factor overrides, the `alpha_tet` alias, bbox sweeps, successful offset
recovery, exhausted/disabled recovery, and conservative plane reduction.

The adversarial containment test wraps isolated corners too tightly and
verifies that original points between the corners cause a failed exit even
when random validation has zero samples. Occupied boxes containing raw points
is not sufficient to certify a wrap of discrete corners. Successful output
means all supplied samples pass the geometric checks (tree classification uses
the configured numerical tolerance); it does not prove containment of an
unknown continuous physical surface between those samples.
