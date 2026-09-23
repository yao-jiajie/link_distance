#pragma once

#include <rokae_demo/octree_rectangular_pruning.hpp>

namespace rokae_demo
{
struct OctreeConvexOptions
{
  // Both are explicit CLI choices. The ratio bounds BOTH each candidate's
  // added volume / source volume and the complete final union's inflation.
  double max_volume_inflation = 0;
  double max_fill_distance = 0; // meters, relative to occupied voxel union
  std::size_t max_candidate_voxels = 128;
  std::size_t max_certificate_boxes = 4096;
};

struct OctreeConvexCell
{
  int source_node = -1;
  int baseline_cell = -1; // >=0 for retained exact rect cell; -1 for a hull
  bool non_aabb = false;
  double volume = 0, fill_distance_bound = 0;
  // Visualization of the actual exported halfspace intersection. Double mesh
  // vertices are display data; planes + exact recomputation define geometry.
  std::vector<HS::Vec3> mesh_vertices;
  std::vector<std::array<std::size_t,3>> mesh_triangles;
};

struct OctreeConvexResult
{
  std::vector<ConvexCluster> clusters;
  std::vector<HS::Vec3> source_vertices;
  std::vector<OctreeConvexCell> cells;
  std::vector<std::size_t> voxel_cluster_ids;
  std::size_t candidate_nodes = 0, candidate_hulls = 0, accepted_merges = 0;
  std::size_t rejected_volume = 0, rejected_distance = 0, rejected_gain = 0, rejected_work_limit = 0;
  std::size_t certificate_boxes = 0;
  double occupied_volume = 0, output_volume = 0, volume_inflation = 0, maximum_fill_distance_bound = 0;
};

// Bottom-up within the existing hierarchy, initialized from its exact rect
// partition. No pairwise global merge, Nef, or traditional convex decomposition.
OctreeConvexResult approximateOctreeConvex(const OctreeOccupancy&, const OctreeCellPartition&,
                                          const OctreeConvexOptions&);

// Reconstruct actual halfspace intersections with exact arithmetic, recompute
// volume/distance certificates and check complete voxel/owner coverage.
// Throws on any failed certificate; random probes are additional, not proofs.
void validateOctreeConvex(const OctreeOccupancy&, const OctreeCellPartition&,
                         const OctreeConvexOptions&, const OctreeConvexResult&);

int runOctreeConvexPipeline(int argc, char** argv);
}  // namespace rokae_demo
