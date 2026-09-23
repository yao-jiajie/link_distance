#pragma once
#include <rokae_demo/boundary_local_tree.hpp>
#include <memory>

namespace rokae_demo
{
struct BoundaryDirectOptions
{
  std::size_t max_canonical_cells=16384; // includes unbounded outer intervals
  std::size_t max_split_checks=1048576;
  std::size_t max_nodes=16384; // allocated expression nodes, including temporary ones
  std::size_t max_expanded_nodes=262144;
  unsigned max_depth=128; // both partition recursion and expression depth
  bool intern_subexpressions=false; // exact ordered sharing during construction
  bool materialize_expression=true; // diagnostics/legacy only; serialized DAG is executable
  bool retain_partition_diagnostics=true; // split list only; executable DAG is unaffected
  BoundaryLocalOptions signs;
};
struct BoundaryDirectSplit
{
  std::array<std::size_t,3> lower{}, upper{}; // half-open canonical cell index range
  std::size_t support=0, cut=0, left_occupied=0, right_occupied=0;
};
struct BoundaryDirectResult
{
  BoundaryPlaneRegistry registry;
  std::array<std::size_t,3> dimensions{};
  std::vector<unsigned char> occupied; // boolean canonical cells, outer cells FREE
  std::vector<std::pair<std::size_t,int>> leaf_sources;
  TreeRepresentation tree;
  // Postorder snapshot produced once during direct construction. Consumers
  // compiling the same frame must not traverse ExpressionPtr again.
  HS::SerializedTree serialized;
  std::vector<BoundaryDirectSplit> splits;
  std::size_t split_checks=0, allocated_nodes=0, constructed_nodes=0;
  std::size_t expanded_nodes=0, expanded_plane_refs=0, free_terminals=0, occupied_terminals=0;
  std::size_t phi_full=0, phi_folded=0;
  bool construction_interned=false;
  unsigned recursion_depth=0, expression_depth=0;
  double registry_ms=0, canonical_ms=0, tree_ms=0, prefix_ms=0;
  double recursion_ms=0, serialization_ms=0, core_ms=0;
};

// No sweep/prisms/exterior branches or fallback. Work caps throw; no partial tree.
BoundaryDirectResult buildBoundaryDirectTree(const SparseVoxelOccupancy&,const OctreeBoundary&,
                                             const BoundaryDirectOptions& = {});
// Core-only overload. It bypasses boundary face materialization/rectangle
// merging and derives the same support-coordinate arrangement from occupancy.
BoundaryDirectResult buildBoundaryDirectTree(const SparseVoxelOccupancy&,
                                             const BoundaryDirectOptions& = {});// Boundary provenance + canonical occupancy + raw/simplified all-strata signs.
// A skipped or failed certificate is NOT permission to publish the direct tree.
BoundarySignDiagnostics validateBoundaryDirectTree(const SparseVoxelOccupancy&,const OctreeBoundary&,
                                                   const BoundaryDirectResult&,const BoundaryDirectOptions& = {});
// Expected geometry sign at finite world coordinates, after the above validation.
int boundaryDirectExpectedSign(const BoundaryDirectResult&,double voxel_size,const HS::Vec3&);
void writeBoundaryDirectDiagnostics(const std::filesystem::path&,const BoundaryDirectResult&,const BoundarySignDiagnostics&);
int runBoundaryDirectPipeline(int argc,char** argv);
} // namespace rokae_demo
