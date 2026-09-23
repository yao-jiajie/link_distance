#pragma once

#include <rokae_demo/octree_occupancy.hpp>

namespace rokae_demo
{
struct OctreeRectangularCell
{
  int source_node = -1;   // existing full cube; -1 for a virtual rectangular cell
  int parent_node = -1;   // width-2 sibling parent for a virtual cell
  unsigned sibling_mask = 0;
  VoxelIndex lower{}, extent{};
  OctreeBox box;
  std::size_t occupied_count = 0;
};

// Exact cubic cut remains the spatial hierarchy and independent reference.
// Rectangular output cells are an overlay, not fictitious cubic octree nodes.
struct OctreeCellPartition
{
  OctreePruning base;
  bool rectangles = false;
  std::vector<OctreeRectangularCell> cells;
  std::vector<std::size_t> voxel_cell_ids;
  std::vector<std::size_t> base_cell_ids;  // cubic-cut cell -> final convex cell
  std::size_t rectangular_groups = 0;
  std::size_t rectangular_merge_count = 0;  // number of non-cubic output rectangles
};

// Minimum disjoint rectangular cover of the given 2x2x2 occupancy mask.
// Bits use x + 2*y + 4*z. Fixed 27 candidates / 256-state lookup, stable ties.
std::vector<unsigned> rectangularSiblingPartition(unsigned occupied_mask);

OctreeCellPartition partitionOctreeCells(const OctreeOccupancy& octree,
                                         OctreePruning base, bool rectangles);
void validateOctreeCellPartition(const OctreeOccupancy& octree,
                                const OctreeCellPartition& partition);
std::vector<ConvexCluster> octreeAabbClusters(const OctreeCellPartition& partition,
                                             std::vector<HS::Vec3>& vertices);
}  // namespace rokae_demo
