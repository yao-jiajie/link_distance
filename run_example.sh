#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_dir="${project_dir}/build"
output_dir="${project_dir}/output"

# Project-wide geometry unit: meter (m). The sample cube spans [0,1]^3 and is
# therefore a 1 m x 1 m x 1 m cube.

cmake -S "${project_dir}" -B "${build_dir}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCGAL_DIR="${project_dir}/_deps/cgal/lib/cmake/CGAL" \
  -DBOOST_ROOT="${project_dir}/_deps/sysroot/usr" \
  -DBoost_NO_SYSTEM_PATHS=ON \
  -DGMP_INCLUDE_DIR="${project_dir}/_deps/sysroot/usr/include/x86_64-linux-gnu" \
  -DGMP_LIBRARIES="/lib/x86_64-linux-gnu/libgmp.so.10" \
  -DMPFR_INCLUDE_DIR="${project_dir}/_deps/sysroot/usr/include" \
  -DMPFR_LIBRARIES="/lib/x86_64-linux-gnu/libmpfr.so.6"
cmake --build "${build_dir}" --parallel 2

mkdir -p "${output_dir}"
"${build_dir}/alpha_wrap_test" \
  "${project_dir}/data/cube_points.xyz" \
  "${output_dir}/cube_wrap.off" \
  --alpha-ratio 2 --offset-ratio 20

"${build_dir}/build_halfspace_tree" \
  "${project_dir}/data/cube_points.xyz" \
  "${output_dir}/cube_alpha_tet" \
  --pipeline alpha-tet \
  --alpha-ratio 2 \
  --offset-ratio 20 \
  --plane-reduction exact \
  --validation-samples 1000

# Explicit Jet benchmark; the engineering main workflow below uses voxels.
"${build_dir}/build_halfspace_tree" \
  "${project_dir}/data/cube_points.xyz" \
  "${output_dir}/cube_sparse_jet" \
  --pipeline alpha-tet \
  --preprocess jet \
  --smooth-neighbors auto \
  --smooth-iterations 1 \
  --alpha-ratio 2 \
  --offset-ratio 20 \
  --plane-reduction exact \
  --validation-samples 1000

"${build_dir}/simulate_depth_camera_cloud" \
  "${project_dir}/data/realistic_sparse_camera.xyz"

"${build_dir}/voxelize_point_cloud" \
  "${project_dir}/data/realistic_sparse_camera.xyz" \
  "${output_dir}/realistic_sparse_camera_voxels" --voxel-size 0.03

"${build_dir}/build_halfspace_tree" \
  "${project_dir}/data/realistic_sparse_camera.xyz" \
  "${output_dir}/realistic_sparse_camera_voxel" \
  --pipeline alpha-tet \
  --voxel-size 0.03 \
  --alpha-voxel-factor 3 \
  --offset-voxel-factor 1 \
  --plane-reduction exact \
  --validation-samples 1000

"${build_dir}/build_halfspace_tree" \
  "${output_dir}/cube_alpha_tet/wrap.off" \
  "${output_dir}/cube_exact_nef" \
  --pipeline exact-nef \
  --point-count 26 \
  --alpha 0.8660254037844386 \
  --offset 0.086602540378443865 \
  --validation-samples 1000 \
  --probe 0.5 0.5 0.5 inside \
  --probe 2 2 2 outside
