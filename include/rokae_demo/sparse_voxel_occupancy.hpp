#pragma once
#include <rokae_demo/octree_occupancy.hpp>

namespace rokae_demo
{
// Fixed-resolution occupancy only: no octree nodes, corners, or point-ID lists.
struct SparseVoxelOccupancy
{
  double voxel_size = 0;
  std::vector<VoxelIndex> occupied_voxels;
  VoxelIndex lower{}, upper{}; // inclusive occupied key bounds
  // Optional lexicographic dense bitset, built only for a compact key span.
  // It avoids sorting duplicate raw keys and accelerates closed-cell queries.
  std::array<std::size_t,3> dense_extent{};
  std::vector<std::uint64_t> dense_bits;
  // Optional boundary-plane sign masks produced while enumerating dense bits.
  // Axis vector index is the plane offset from lower[axis].
  std::array<std::vector<unsigned char>,3> boundary_plane_signs;
  OctreeBox box, sampling_box; // sampling_box matches the legacy root extent
  unsigned required_depth = 0; // informational; only fallback needs max_depth
  bool contains(const HS::Vec3&) const;
};
SparseVoxelOccupancy buildSparseVoxelOccupancy(const std::vector<VoxelPoint>&, double voxel_size);
void validateSparseVoxelOccupancy(const SparseVoxelOccupancy&, const std::vector<VoxelPoint>&);
} // namespace rokae_demo
