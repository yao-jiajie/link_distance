#pragma once
#include <rokae_demo/octree_boundary.hpp>
#include <rokae_demo/sparse_voxel_occupancy.hpp>
#include <functional>

namespace rokae_demo
{
struct BoundarySupportPlane
{
  unsigned axis = 0;
  std::int64_t coordinate = 0;
  unsigned outward_signs = 0; // bit 0: negative, bit 1: positive
  std::vector<std::size_t> patches;
};
struct BoundaryPlaneRegistry
{
  std::array<std::vector<std::int64_t>,3> coordinates;
  std::vector<BoundarySupportPlane> supports;
  bool patch_provenance_complete = true;
};
// Geometry-only registry: no sweep, free rectangles, prisms or tree construction.
BoundaryPlaneRegistry buildBoundaryPlaneRegistry(const OctreeBoundary&);
// Core distance-field path: derive the identical ordered support planes and
// outward signs directly from occupied voxels, without materializing/merging
// rectangular boundary patches. Patch provenance remains intentionally empty.
BoundaryPlaneRegistry buildBoundaryPlaneRegistry(const SparseVoxelOccupancy&);
struct FreePrism3D { VoxelIndex lower{}, upper{}; };
struct BoundaryTreeOptions
{
  unsigned sweep_axis = 0; // 0=X, 1=Y, 2=Z; local transverse axes are cyclic
  std::size_t max_compressed_cells = 16384;
  std::size_t max_event_updates = 1048576;
  std::size_t max_free_prisms = 4096;
  unsigned max_merge_passes = 4;
};
struct BoundaryTreeResult
{
  unsigned sweep_axis = 0;
  std::array<std::vector<std::int64_t>,3> coordinates;
  std::vector<BoundarySupportPlane> supports;
  std::vector<FreePrism3D> prisms;
  std::vector<std::pair<std::size_t,int>> leaf_sources; // support ID + coefficient sign
  TreeRepresentation tree;
  bool fallback = false, merge_fixed_point = false;
  std::string fallback_reason;
  std::size_t compressed_cells = 0, event_updates = 0, prisms_before_merge = 0;
  unsigned merge_passes = 0;
  double registry_ms = 0, sweep_ms = 0, partition_ms = 0, merge_ms = 0, tree_ms = 0;
};

struct BoundaryTreeBatch
{
  std::vector<BoundaryTreeResult> candidates;
  std::vector<double> candidate_build_ms;
  std::size_t selected = 0;
  double all_candidates_ms = 0, selection_ms = 0;
};

// Single-axis sweep. Work-cap fallback keeps the original obstacle rect tree.
BoundaryTreeResult buildBoundaryTree(const OctreeOccupancy&, const OctreeCellPartition&,
                                    const OctreeBoundary&, const BoundaryTreeOptions& = {});
// best_axis=true builds all X/Y/Z candidates. Successful boundary candidates
// outrank cap fallbacks; then minimize (nodes, plane refs, regions, fixed axis).
// If all axes hit a work cap, select the original rect fallback from X.
BoundaryTreeBatch buildBoundaryTreeBatch(const OctreeOccupancy&, const OctreeCellPartition&,
                                        const OctreeBoundary&, const BoundaryTreeOptions&, bool best_axis);
// Shares the identical sweep/tree implementation. Invokes fallback_tree ONLY
// after a sweep work cap, allowing the caller to lazily construct legacy rect.
BoundaryTreeBatch buildVoxelBoundaryTreeBatch(double voxel_size, const OctreeBoundary&,
                                            const BoundaryTreeOptions&, bool best_axis,
                                            const std::function<TreeRepresentation()>& fallback_tree);
// Nonfallback voxel path only. Fallback is validated through the legacy API.
void validateVoxelBoundaryTree(const SparseVoxelOccupancy&, const OctreeBoundary&, const BoundaryTreeResult&);
bool boundaryStrictFree(const BoundaryTreeResult&, const SparseVoxelOccupancy&, const HS::Vec3&,
                        bool* used_zero_guard = nullptr);
// Complete integer compressed-cell coverage proof, event/boundary provenance,
// and tree coefficient/structure check. Throws on certificate failure.
void validateBoundaryTree(const OctreeOccupancy&, const OctreeCellPartition&,
                         const OctreeBoundary&, const BoundaryTreeResult&);
// Strict R^3 \ V. NOT equivalent to evaluateInside(tree)! The nonfallback tree
// represents closure(R^3 \ V); zero values require the closed occupancy oracle.
// The tree and occupancy argument must originate from the SAME cloud snapshot.
bool boundaryStrictFree(const BoundaryTreeResult&, const OctreeOccupancy&, const HS::Vec3&,
                        bool* used_zero_guard = nullptr);
int runOctreeBoundaryTreePipeline(int argc, char** argv);
} // namespace rokae_demo
