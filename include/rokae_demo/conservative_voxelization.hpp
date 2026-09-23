#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace rokae_demo
{

using VoxelPoint = CGAL::Exact_predicates_inexact_constructions_kernel::Point_3;
using VoxelIndex = std::array<std::int64_t, 3>;

struct VoxelizationStats
{
  double voxel_size = 0.0;  // meters; grid origin is always (0, 0, 0)
  std::size_t raw_points = 0;
  std::size_t occupied_voxels = 0;
  std::size_t voxel_points = 0;
  double voxel_point_ratio = 0.0;
  std::size_t raw_points_outside_voxels = 0;
  double construction_ms = 0.0;
  double validation_ms = 0.0;
  double total_ms = 0.0;
};

struct ConservativeVoxelization
{
  double voxel_size = 0.0;
  // Lexicographically sorted, unique integer keys; independent of input order.
  std::vector<VoxelIndex> occupied_voxels;
  std::vector<VoxelIndex> corner_indices;
  std::vector<VoxelPoint> points;  // one point per corner_indices entry
  VoxelizationStats stats;

  // Tests membership in the CLOSED union of occupied boxes, not in the
  // discrete corner set. This does not certify a subsequent Alpha Wrap.
  bool contains(const VoxelPoint& point) const;
  std::size_t countOutside(const std::vector<VoxelPoint>& raw_points) const;
};

// Half-open indexing, including negative coordinates. Grid planes are rounded
// double(i * voxel_size); boundary comparisons use those same exported planes.
// Rejects non-finite inputs and unrepresentable grid ranges rather than clamping.
VoxelIndex voxelIndex(const VoxelPoint& point, double voxel_size);

ConservativeVoxelization conservativeVoxelize(
    const std::vector<VoxelPoint>& raw_points, double voxel_size);

}  // namespace rokae_demo
