#pragma once

#include <rokae_demo/conservative_voxelization.hpp>
#include <rokae_demo/convex_cluster.hpp>

namespace rokae_demo
{

struct OctreeBox
{
  std::array<double, 3> lower{}, upper{};
  bool contains(const HS::Vec3& point) const;
  double volume() const;
};

// Common rounded-grid bounds for cubic and rectangular occupied cells.
OctreeBox octreeGridBox(const VoxelIndex& lower, const VoxelIndex& extent, double size);

struct OctreeNode
{
  VoxelIndex lower{};
  std::int64_t width = 1;  // integer fine-voxel units
  unsigned depth = 0;
  OctreeBox box;
  std::array<int, 8> children{{-1, -1, -1, -1, -1, -1, -1, -1}};
  std::size_t occupied_count = 0;
  std::size_t point_count = 0;
  // Point indices are stored only at leaves, never copied at every level.
  std::vector<std::size_t> point_indices;
};

struct OctreeOccupancy
{
  double voxel_size = 0;
  unsigned required_depth = 0;
  std::vector<VoxelIndex> occupied_voxels;  // sorted unique keys
  std::vector<OctreeNode> nodes;           // root is nodes[0]
  std::vector<int> leaf_nodes;             // same order as occupied_voxels
  bool contains(const HS::Vec3& point) const;
};

// Phase 1: sparse storage but fixed fine resolution, WITHOUT parent pruning.
OctreeOccupancy buildOctreeOccupancy(const std::vector<VoxelPoint>& raw_points,
                                     double voxel_size, unsigned max_depth = 20);

// A cut through the original hierarchy. Descendants of selected cells are
// inactive and omitted from the exported hierarchy; original occupancy remains
// available as an independent validation/provenance reference.
struct OctreePruning
{
  bool exact = true;
  std::vector<int> cell_nodes;  // ordered by first covered fine voxel (lexicographic)
  std::vector<std::size_t> voxel_cell_ids;  // fine voxel -> output cluster ID
  std::vector<unsigned char> active_nodes;
  std::size_t merge_operations = 0;  // each accepted eight-child -> parent recovery
};

// O(nodes + occupied voxels), no pairwise merges, no floating volume thresholds.
// exact=false retains the Phase 1 cut, one cell per fine occupied voxel.
OctreePruning pruneOctree(const OctreeOccupancy& octree, bool exact = true);

// Requires a valid source occupancy (validateOctreeOccupancy). Checks full
// coverage, no overlap, full parents only, maximal exact recovery, and provenance.
void validateOctreePruning(const OctreeOccupancy& octree, const OctreePruning& pruning);

// Analytic six-plane extraction. Each AABB is already irredundant, so exact
// local reduction is the identity; the shared backend deduplicates plane IDs.
std::vector<ConvexCluster> octreeAabbClusters(const OctreeOccupancy& octree,
                                             std::vector<HS::Vec3>& vertices);
std::vector<ConvexCluster> octreeAabbClusters(const OctreeOccupancy& octree,
                                             const OctreePruning& pruning,
                                             std::vector<HS::Vec3>& vertices);
std::vector<ConvexCluster> aabbClusters(const std::vector<OctreeBox>& boxes,
                                       const std::vector<std::size_t>& voxel_cell_ids,
                                       std::vector<HS::Vec3>& vertices);

// Structural certificate: one box per unique occupied key, all child partitions
// and point assignments valid. Throws on failure; no randomized proof shortcut.
void validateOctreeOccupancy(const OctreeOccupancy& octree,
                            const std::vector<VoxelPoint>& raw_points);

}  // namespace rokae_demo
