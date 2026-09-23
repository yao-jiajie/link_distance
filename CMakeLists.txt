cmake_minimum_required(VERSION 3.16)

project(cgal_halfspace_tools LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(CGAL 5.6.3 EXACT REQUIRED)
find_package(Eigen3 3.2 REQUIRED NO_MODULE)
find_package(Boost REQUIRED)
find_package(Threads REQUIRED)
include(CGAL_Eigen3_support)

option(ROKAE_USE_ABSEIL_BOUNDARY_HASH "Use Abseil flat_hash_set for occupied-voxel neighbor lookup" OFF)
option(ROKAE_USE_ABSEIL_DIRECT_HASH "Use Abseil flat_hash_map for direct-tree interning" OFF)
if(ROKAE_USE_ABSEIL_BOUNDARY_HASH OR ROKAE_USE_ABSEIL_DIRECT_HASH)
  find_package(absl CONFIG REQUIRED)
endif()

add_library(conservative_voxelization src/conservative_voxelization.cpp)
target_include_directories(conservative_voxelization PUBLIC include)
target_link_libraries(conservative_voxelization PUBLIC CGAL::CGAL)
target_compile_options(conservative_voxelization PRIVATE -Wall -Wextra -Wpedantic)

add_library(convex_cluster src/convex_cluster.cpp)
target_include_directories(convex_cluster PUBLIC include)
target_compile_options(convex_cluster PRIVATE -Wall -Wextra -Wpedantic)

add_library(octree_occupancy src/octree_occupancy.cpp src/octree_rectangular_pruning.cpp src/octree_boundary.cpp src/sparse_voxel_occupancy.cpp)
target_link_libraries(octree_occupancy PUBLIC conservative_voxelization convex_cluster)
target_compile_options(octree_occupancy PRIVATE -Wall -Wextra -Wpedantic)
if(ROKAE_USE_ABSEIL_BOUNDARY_HASH)
  target_link_libraries(octree_occupancy PRIVATE absl::flat_hash_set)
  target_compile_definitions(octree_occupancy PRIVATE ROKAE_USE_ABSEIL_BOUNDARY_HASH=1)
endif()

add_library(octree_convex_approximation src/octree_convex_approximation.cpp)
target_link_libraries(octree_convex_approximation PUBLIC octree_occupancy)
target_compile_options(octree_convex_approximation PRIVATE -Wall -Wextra -Wpedantic)

add_library(octree_boundary_tree src/octree_boundary_tree.cpp src/boundary_local_tree.cpp src/boundary_direct_tree.cpp src/boundary_smooth_distance.cpp src/boundary_direct_frame_cache.cpp)
target_link_libraries(octree_boundary_tree PUBLIC octree_occupancy Threads::Threads)
target_compile_options(octree_boundary_tree PRIVATE -Wall -Wextra -Wpedantic)
if(ROKAE_USE_ABSEIL_DIRECT_HASH)
  target_link_libraries(octree_boundary_tree PRIVATE absl::flat_hash_map)
  target_compile_definitions(octree_boundary_tree PRIVATE ROKAE_USE_ABSEIL_DIRECT_HASH=1)
endif()

add_library(spherical_robot_model src/spherical_robot_model.cpp)
target_include_directories(spherical_robot_model PUBLIC include PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(spherical_robot_model PUBLIC Eigen3::Eigen)
target_compile_options(spherical_robot_model PRIVATE -Wall -Wextra -Wpedantic)

add_library(link_distance_evaluator src/link_distance_evaluator.cpp)
target_include_directories(link_distance_evaluator PUBLIC include)
target_link_libraries(link_distance_evaluator PUBLIC spherical_robot_model Threads::Threads)
target_compile_options(link_distance_evaluator PRIVATE -Wall -Wextra -Wpedantic)

add_executable(voxelize_point_cloud src/voxelize_point_cloud.cpp)
target_link_libraries(voxelize_point_cloud PRIVATE conservative_voxelization)
target_compile_options(voxelize_point_cloud PRIVATE -Wall -Wextra -Wpedantic)

include(CTest)
if(BUILD_TESTING)
  add_executable(halfspace_smooth_tree_test tests/halfspace_smooth_tree_test.cpp)
  target_include_directories(halfspace_smooth_tree_test PRIVATE include)
  target_link_libraries(halfspace_smooth_tree_test PRIVATE Threads::Threads)
  target_compile_options(halfspace_smooth_tree_test PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME halfspace_smooth_tree COMMAND halfspace_smooth_tree_test)
  add_executable(spherical_robot_model_test tests/spherical_robot_model_test.cpp)
  target_link_libraries(spherical_robot_model_test PRIVATE link_distance_evaluator)
  target_compile_options(spherical_robot_model_test PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME spherical_robot_model COMMAND spherical_robot_model_test
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/spherical_robot_test.urdf"
    "${CMAKE_CURRENT_SOURCE_DIR}/resources/vamp/panda/panda_spherized.urdf")
  add_executable(boundary_smooth_distance_test tests/boundary_smooth_distance_test.cpp)
  target_link_libraries(boundary_smooth_distance_test PRIVATE octree_boundary_tree)
  target_compile_options(boundary_smooth_distance_test PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME boundary_smooth_distance COMMAND boundary_smooth_distance_test)
  set_tests_properties(boundary_smooth_distance PROPERTIES TIMEOUT 120)
  add_executable(halfspace_compiled_tree_test tests/halfspace_compiled_tree_test.cpp)
  target_include_directories(halfspace_compiled_tree_test PRIVATE include)
  target_link_libraries(halfspace_compiled_tree_test PRIVATE Threads::Threads)
  target_compile_options(halfspace_compiled_tree_test PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME halfspace_compiled_tree COMMAND halfspace_compiled_tree_test)
  set_tests_properties(halfspace_compiled_tree PROPERTIES TIMEOUT 120)
  add_executable(boundary_direct_tree_test tests/boundary_direct_tree_test.cpp)
  target_link_libraries(boundary_direct_tree_test PRIVATE octree_boundary_tree)
  target_compile_options(boundary_direct_tree_test PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME boundary_direct_tree COMMAND boundary_direct_tree_test)
  set_tests_properties(boundary_direct_tree PROPERTIES TIMEOUT 120)
  add_executable(boundary_direct_frame_cache_test tests/boundary_direct_frame_cache_test.cpp)
  target_link_libraries(boundary_direct_frame_cache_test PRIVATE octree_boundary_tree)
  target_compile_options(boundary_direct_frame_cache_test PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME boundary_direct_frame_cache COMMAND boundary_direct_frame_cache_test)
  set_tests_properties(boundary_direct_frame_cache PROPERTIES TIMEOUT 120)
  add_executable(sparse_voxel_occupancy_test tests/sparse_voxel_occupancy_test.cpp)
  target_link_libraries(sparse_voxel_occupancy_test PRIVATE octree_boundary_tree)
  target_compile_options(sparse_voxel_occupancy_test PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME sparse_voxel_occupancy COMMAND sparse_voxel_occupancy_test)
  set_tests_properties(sparse_voxel_occupancy PROPERTIES TIMEOUT 120)
  add_executable(boundary_local_tree_test tests/boundary_local_tree_test.cpp)
  target_link_libraries(boundary_local_tree_test PRIVATE octree_boundary_tree)
  target_compile_options(boundary_local_tree_test PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME boundary_local_tree COMMAND boundary_local_tree_test)
  set_tests_properties(boundary_local_tree PROPERTIES TIMEOUT 120)
  add_executable(octree_boundary_tree_test tests/octree_boundary_tree_test.cpp)
  target_link_libraries(octree_boundary_tree_test PRIVATE octree_boundary_tree)
  target_compile_options(octree_boundary_tree_test PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME octree_boundary_tree COMMAND octree_boundary_tree_test)
  set_tests_properties(octree_boundary_tree PROPERTIES TIMEOUT 120)
  add_executable(octree_convex_approximation_test tests/octree_convex_approximation_test.cpp)
  target_link_libraries(octree_convex_approximation_test PRIVATE octree_convex_approximation)
  target_compile_options(octree_convex_approximation_test PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME octree_convex_approximation COMMAND octree_convex_approximation_test)
  set_tests_properties(octree_convex_approximation PROPERTIES TIMEOUT 120)
  add_executable(conservative_voxelization_test tests/conservative_voxelization_test.cpp)
  target_link_libraries(conservative_voxelization_test PRIVATE conservative_voxelization)
  add_test(NAME conservative_voxelization COMMAND conservative_voxelization_test)
  add_executable(octree_occupancy_test tests/octree_occupancy_test.cpp)
  target_link_libraries(octree_occupancy_test PRIVATE octree_occupancy)
  target_compile_options(octree_occupancy_test PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME octree_occupancy COMMAND octree_occupancy_test)
  add_executable(octree_pruning_test tests/octree_pruning_test.cpp)
  target_link_libraries(octree_pruning_test PRIVATE octree_occupancy)
  target_compile_options(octree_pruning_test PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME octree_pruning COMMAND octree_pruning_test)
  add_executable(octree_rectangular_pruning_test tests/octree_rectangular_pruning_test.cpp)
  target_link_libraries(octree_rectangular_pruning_test PRIVATE octree_occupancy)
  target_compile_options(octree_rectangular_pruning_test PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME octree_rectangular_pruning COMMAND octree_rectangular_pruning_test)
  add_executable(octree_boundary_test tests/octree_boundary_test.cpp)
  target_link_libraries(octree_boundary_test PRIVATE octree_occupancy)
  target_compile_options(octree_boundary_test PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME octree_boundary COMMAND octree_boundary_test)
  add_executable(boundary_frontend_test tests/boundary_frontend_test.cpp)
  target_link_libraries(boundary_frontend_test PRIVATE octree_boundary_tree)
  target_compile_options(boundary_frontend_test PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME boundary_frontend COMMAND boundary_frontend_test)
  set_tests_properties(boundary_frontend PROPERTIES TIMEOUT 120)
  add_executable(boundary_frontend_benchmark tests/boundary_frontend_benchmark.cpp)
  target_link_libraries(boundary_frontend_benchmark PRIVATE octree_boundary_tree)
  target_compile_options(boundary_frontend_benchmark PRIVATE -Wall -Wextra -Wpedantic)
  add_executable(voxelization_benchmark tests/voxelization_benchmark.cpp)
  target_link_libraries(voxelization_benchmark PRIVATE conservative_voxelization)
  add_executable(halfspace_tree_benchmark tests/halfspace_tree_benchmark.cpp)
  target_include_directories(halfspace_tree_benchmark PRIVATE include ${Boost_INCLUDE_DIRS})
  target_compile_options(halfspace_tree_benchmark PRIVATE -Wall -Wextra -Wpedantic)
  add_executable(boundary_direct_frame_benchmark tests/boundary_direct_frame_benchmark.cpp)
  target_link_libraries(boundary_direct_frame_benchmark PRIVATE octree_boundary_tree)
  target_compile_options(boundary_direct_frame_benchmark PRIVATE -Wall -Wextra -Wpedantic)
  add_executable(halfspace_batch_benchmark tests/halfspace_batch_benchmark.cpp)
  target_link_libraries(halfspace_batch_benchmark PRIVATE octree_boundary_tree)
  target_compile_options(halfspace_batch_benchmark PRIVATE -Wall -Wextra -Wpedantic)
  add_executable(link_distance_kernel_benchmark tests/link_distance_kernel_benchmark.cpp)
  target_link_libraries(link_distance_kernel_benchmark PRIVATE link_distance_evaluator octree_boundary_tree)
  target_compile_options(link_distance_kernel_benchmark PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME link_distance_kernel_benchmark_smoke
    COMMAND link_distance_kernel_benchmark
      "${CMAKE_CURRENT_SOURCE_DIR}/resources/vamp/panda/panda_spherized.urdf"
      "${CMAKE_CURRENT_SOURCE_DIR}/data/realistic_sparse_camera.xyz"
      0.02 0.3 7)
  set_tests_properties(link_distance_kernel_benchmark_smoke PROPERTIES TIMEOUT 30)
  add_executable(link_distance_frame_benchmark tests/link_distance_frame_benchmark.cpp)
  target_link_libraries(link_distance_frame_benchmark PRIVATE link_distance_evaluator octree_boundary_tree)
  target_compile_options(link_distance_frame_benchmark PRIVATE -Wall -Wextra -Wpedantic)
  add_test(NAME link_distance_frame_benchmark_smoke
    COMMAND link_distance_frame_benchmark
      "${CMAKE_CURRENT_SOURCE_DIR}/resources/vamp/panda/panda_spherized.urdf"
      "${CMAKE_CURRENT_SOURCE_DIR}/data/panda_complex_scene.xyz"
      0.02 2)
  set_tests_properties(link_distance_frame_benchmark_smoke PROPERTIES TIMEOUT 30)
  find_package(Python3 COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    add_test(NAME point_cloud_input_fixtures
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/point_cloud_input_fixtures_test.py")
    set_tests_properties(point_cloud_input_fixtures PROPERTIES TIMEOUT 60)
    add_test(NAME link_distance_node_smoke
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/link_distance_node_test.py"
        $<TARGET_FILE:link_distance_node>
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/spherical_robot_test.urdf"
        "${CMAKE_CURRENT_SOURCE_DIR}/data/cube_points.xyz"
        "${CMAKE_CURRENT_SOURCE_DIR}/scripts/visualize_link_distances.py"
        "${CMAKE_CURRENT_SOURCE_DIR}/scripts/visualize_link_distances_3d.py"
        "${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_sphere_joint_frames.py")
    set_tests_properties(link_distance_node_smoke PROPERTIES TIMEOUT 60)
    add_test(NAME link_distance_stream
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/link_distance_stream_test.py"
        $<TARGET_FILE:link_distance_node>
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/spherical_robot_test.urdf"
        "${CMAKE_CURRENT_SOURCE_DIR}/data/cube_points.xyz")
    set_tests_properties(link_distance_stream PROPERTIES TIMEOUT 60)
    add_test(NAME link_distance_panda_001
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/link_distance_panda_001_test.py"
        $<TARGET_FILE:link_distance_node>
        "${CMAKE_CURRENT_SOURCE_DIR}/resources/vamp/panda/panda_spherized.urdf"
        "${CMAKE_CURRENT_SOURCE_DIR}/data/realistic_sparse_camera.xyz")
    set_tests_properties(link_distance_panda_001 PROPERTIES TIMEOUT 30)
    add_test(NAME panda_complex_scene
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/panda_complex_scene_test.py"
        $<TARGET_FILE:link_distance_node>
        "${CMAKE_CURRENT_SOURCE_DIR}/resources/vamp/panda/panda_spherized.urdf"
        "${CMAKE_CURRENT_SOURCE_DIR}/data/panda_complex_scene.xyz")
    set_tests_properties(panda_complex_scene PROPERTIES TIMEOUT 90)
    add_test(NAME point_cloud_inputs
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/point_cloud_inputs_test.py"
        $<TARGET_FILE_DIR:build_halfspace_tree>)
    set_tests_properties(point_cloud_inputs PROPERTIES TIMEOUT 180)
    add_test(NAME boundary_frontend_benchmark
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/boundary_frontend_benchmark_test.py"
        $<TARGET_FILE_DIR:boundary_frontend_benchmark>)
    set_tests_properties(boundary_frontend_benchmark PROPERTIES TIMEOUT 120)
    add_test(NAME lse_kernel_benchmark
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/lse_kernel_benchmark_test.py"
        $<TARGET_FILE_DIR:build_halfspace_tree>)
    set_tests_properties(lse_kernel_benchmark PROPERTIES TIMEOUT 120)
    add_test(NAME execution_sharing_benchmark
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/execution_sharing_benchmark_test.py"
        $<TARGET_FILE_DIR:build_halfspace_tree>)
    set_tests_properties(execution_sharing_benchmark PROPERTIES TIMEOUT 120)
    add_test(NAME validation_acceleration_benchmark
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/validation_acceleration_benchmark_test.py"
        $<TARGET_FILE_DIR:build_halfspace_tree>)
    set_tests_properties(validation_acceleration_benchmark PROPERTIES TIMEOUT 120)
    add_test(NAME multi_frame_reuse_benchmark
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/multi_frame_reuse_benchmark_test.py"
        $<TARGET_FILE_DIR:boundary_direct_frame_benchmark>)
    set_tests_properties(multi_frame_reuse_benchmark PROPERTIES TIMEOUT 120)
    add_test(NAME halfspace_batch_benchmark_smoke
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/halfspace_batch_benchmark_test.py"
        $<TARGET_FILE:halfspace_batch_benchmark>)
    set_tests_properties(halfspace_batch_benchmark_smoke PROPERTIES TIMEOUT 120)
    add_test(NAME smooth_distance_pipeline
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/smooth_distance_pipeline_test.py"
        $<TARGET_FILE:build_halfspace_tree>)
    set_tests_properties(smooth_distance_pipeline PROPERTIES TIMEOUT 120)
    add_test(NAME smooth_distance_benchmark
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/smooth_distance_benchmark_test.py"
        $<TARGET_FILE_DIR:build_halfspace_tree>)
    set_tests_properties(smooth_distance_benchmark PROPERTIES TIMEOUT 120)
    add_test(NAME boundary_direct_query_benchmark
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/boundary_direct_query_benchmark_test.py"
        $<TARGET_FILE_DIR:build_halfspace_tree>)
    set_tests_properties(boundary_direct_query_benchmark PROPERTIES TIMEOUT 120)
    add_test(NAME boundary_direct_benchmark
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/boundary_direct_benchmark_test.py"
        $<TARGET_FILE_DIR:build_halfspace_tree>)
    set_tests_properties(boundary_direct_benchmark PROPERTIES TIMEOUT 120)
    add_test(NAME boundary_direct_pipeline
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/boundary_direct_pipeline_test.py"
        $<TARGET_FILE:build_halfspace_tree>)
    set_tests_properties(boundary_direct_pipeline PROPERTIES TIMEOUT 120)
    add_test(NAME boundary_voxel_benchmark
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/boundary_voxel_benchmark_test.py"
        $<TARGET_FILE_DIR:build_halfspace_tree>)
    set_tests_properties(boundary_voxel_benchmark PROPERTIES TIMEOUT 120)
    add_test(NAME boundary_voxel_pipeline
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/boundary_voxel_pipeline_test.py"
        $<TARGET_FILE:build_halfspace_tree>)
    set_tests_properties(boundary_voxel_pipeline PROPERTIES TIMEOUT 120)
    add_test(NAME boundary_local_pipeline
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/boundary_local_pipeline_test.py"
        $<TARGET_FILE:build_halfspace_tree>)
    set_tests_properties(boundary_local_pipeline PROPERTIES TIMEOUT 120)
    add_test(NAME boundary_local_benchmark
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/boundary_local_benchmark_test.py"
        $<TARGET_FILE_DIR:build_halfspace_tree>)
    set_tests_properties(boundary_local_benchmark PROPERTIES TIMEOUT 120)
    add_test(NAME octree_boundary_tree_pipeline
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/octree_boundary_tree_pipeline_test.py"
        $<TARGET_FILE:build_halfspace_tree>)
    set_tests_properties(octree_boundary_tree_pipeline PROPERTIES TIMEOUT 120)
    add_test(NAME boundary_xyz_benchmark
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/boundary_benchmark_test.py"
        $<TARGET_FILE_DIR:build_halfspace_tree>)
    set_tests_properties(boundary_xyz_benchmark PROPERTIES TIMEOUT 120)
    add_test(NAME octree_convex_pipeline
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/octree_convex_pipeline_test.py"
        $<TARGET_FILE:build_halfspace_tree>)
    set_tests_properties(octree_convex_pipeline PROPERTIES TIMEOUT 120)
    add_test(NAME phase4_benchmark
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/phase4_benchmark_test.py"
        $<TARGET_FILE_DIR:build_halfspace_tree>)
    set_tests_properties(phase4_benchmark PROPERTIES TIMEOUT 120)
    add_test(NAME octree_pipeline
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/octree_pipeline_test.py"
        $<TARGET_FILE:build_halfspace_tree>)
    set_tests_properties(octree_pipeline PROPERTIES TIMEOUT 120)
    add_test(NAME octree_boundary_pipeline
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/octree_boundary_pipeline_test.py"
        $<TARGET_FILE:build_halfspace_tree>)
    set_tests_properties(octree_boundary_pipeline PROPERTIES TIMEOUT 120)
    add_test(NAME voxel_pipeline
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/voxel_pipeline_test.py"
        $<TARGET_FILE:build_halfspace_tree> $<TARGET_FILE:alpha_wrap_test>
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/voxel_containment.xyz")
    set_tests_properties(voxel_pipeline PROPERTIES TIMEOUT 180)
  endif()
endif()

add_executable(alpha_wrap_test
  src/alpha_wrap_test.cpp
  src/point_cloud_preprocessing.cpp)
target_include_directories(alpha_wrap_test PRIVATE include)
target_link_libraries(alpha_wrap_test PRIVATE CGAL::CGAL CGAL::Eigen3_support conservative_voxelization)

add_executable(build_halfspace_tree
  src/build_halfspace_tree.cpp
  src/alpha_tet_pipeline.cpp
  src/octree_pipeline.cpp
  src/octree_convex_pipeline.cpp
  src/octree_boundary_tree_pipeline.cpp
  src/boundary_direct_pipeline.cpp
  src/point_cloud_preprocessing.cpp)
target_include_directories(build_halfspace_tree PRIVATE include)
target_link_libraries(build_halfspace_tree PRIVATE CGAL::CGAL CGAL::Eigen3_support conservative_voxelization convex_cluster octree_occupancy octree_convex_approximation octree_boundary_tree)

add_executable(simulate_depth_camera_cloud
  src/simulate_depth_camera_cloud.cpp)

add_executable(link_distance_node src/link_distance_node.cpp)
target_link_libraries(link_distance_node PRIVATE link_distance_evaluator octree_boundary_tree)
target_compile_options(link_distance_node PRIVATE -Wall -Wextra -Wpedantic)

foreach(target alpha_wrap_test build_halfspace_tree simulate_depth_camera_cloud)
  target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
endforeach()
