# Data units

All coordinates in this project are meters (`m`).

- `cube_points.xyz` samples the surface of a `1 m x 1 m x 1 m` cube.
- `l_prism.off` uses the same meter convention.
- `realistic_sparse_camera.xyz` is a deterministic, noisy, multi-view RGB-D
  simulation of a `0.36 m x 0.26 m x 0.48 m` ellipsoid. Its matching ideal
  samples are in `realistic_sparse_camera_ground_truth.xyz`; acquisition and
  noise parameters are recorded in `realistic_sparse_camera_metadata.json`.
- `c_bracket_surface.xyz` is an additional visualization example: an extruded
  nonconvex C-bracket whose inner/outer cylindrical walls are curved and whose
  top, bottom, and two cut faces are planar. `c_bracket_surface.json` records
  the analytic dimensions and per-surface sample counts. Regenerate it with
  `python3 scripts/generate_c_bracket_surface.py`.
- `panda_disconnected_irregular_scene.xyz` combines three disconnected
  irregular surfaces around the fixed Panda example pose: a nonconvex C-bracket,
  a wavy torus, and a capped variable-radius bent tube. Regenerate it with
  `python3 scripts/generate_panda_disconnected_scene.py`; the matching JSON
  records component parameters and the recommended node arguments.
- `point_cloud_cases/` contains 20 reproducible XYZ inputs: the original camera
  cloud plus box/concave/curved/hollow/thin/disconnected surfaces and ellipsoid
  density, noise, partial-view, duplicate and outlier variants. See its
  [manifest](point_cloud_cases/manifest.json) for dimensions and sampling rules,
  [point-cloud preview](point_cloud_cases/preview.svg), and
  [test protocol](../docs/point_cloud_input_tests.md). Synthetic samples have no
  sensor/visibility model except for the copied camera baseline. They sample
  surfaces; the voxel pipeline does not fill the underlying object interiors.

XYZ and OFF do not reliably carry physical-unit metadata. External datasets
must therefore be converted to meters before being passed to the programs.
