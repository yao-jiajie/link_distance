#pragma once
#include <rokae_demo/octree_boundary_tree.hpp>

namespace rokae_demo
{
struct BoundarySignWitness
{
  std::array<std::size_t,3> strata{};
  unsigned dimension = 0;
  int expected = 0, actual = 0;
};
struct BoundarySignDiagnostics
{
  bool applied = false, closure_set_certified = false, strict_sign_certified = false;
  std::string reason;
  std::size_t strata_count = 0, unsafe_free_count = 0, missed_free_count = 0;
  std::size_t boundary_sign_errors = 0, closure_set_errors = 0, raw_simplified_sign_errors = 0;
  // Indexed by dimension: vertex, relative-open edge, relative-open face, open volume.
  std::array<std::size_t,4> free_zeros{}, occupied_zeros{};
  std::vector<BoundarySignWitness> witnesses;
  double time_ms = 0;
};
struct BoundaryLocalOptions
{
  std::size_t max_expanded_nodes = 16384;
  std::size_t max_attempts = 4096;
  std::size_t max_adjacency_checks = 1048576;
  std::size_t max_strata = 200000;
  std::size_t max_sign_operations = 200000000;
  bool packed_sign_propagation = true;
};
// Independent sign certificate for any MIN/MAX tree on this boundary registry.
// Caller must validate raw occupancy and boundary/registry provenance first.
// Coefficients are checked against registered planes, NOT a flat reference tree.
BoundarySignDiagnostics diagnoseBoundaryPlaneSigns(const std::vector<VoxelIndex>& occupied, double voxel_size,
  const BoundaryPlaneRegistry&, const std::vector<std::pair<std::size_t,int>>& leaf_sources,
  const TreeRepresentation&, const BoundaryLocalOptions& = {});
struct BoundaryFaceAdjacency
{
  std::size_t a = 0, b = 0, support = 0;
  int a_sign = 1;
  std::int64_t u0 = 0, u1 = 0, v0 = 0, v1 = 0;
  bool exterior = false;
};
struct BoundaryLocalResult
{
  TreeRepresentation tree;
  std::vector<BoundaryFaceAdjacency> adjacency;
  BoundarySignDiagnostics diagnostics;
  bool obstacle_fallback = false;
  std::string stop_reason;
  std::size_t adjacency_checks = 0, cancelled_interface_pairs = 0, logical_merges = 0;
  std::size_t factoring_failed = 0, attempts = 0, constructed_nodes = 0, expanded_nodes = 0;
  double adjacency_ms = 0, construction_ms = 0, simplify_ms = 0, core_ms = 0;
};

// Finite, symbolic sign certificate on ALL strata of the existing axis-plane
// arrangement, including unbounded intervals. No floating-point sampling or guard.
// Preconditions: validateOctreeOccupancy and validateBoundaryTree have accepted
// the SAME occupancy/reference snapshot. Symbolic proof covers the real-valued
// MIN/MAX signs of the registered planes, not arbitrary floating-point errors.
BoundarySignDiagnostics diagnoseBoundarySigns(const OctreeOccupancy&, const BoundaryTreeResult& reference,
                                              const TreeRepresentation&, const BoundaryLocalOptions& = {});
BoundarySignDiagnostics diagnoseBoundarySigns(const SparseVoxelOccupancy&, const BoundaryTreeResult& reference,
                                              const TreeRepresentation&, const BoundaryLocalOptions& = {});
// Restricted construction templates only. An incomplete sign repair is a valid
// outcome, explicitly reported by diagnostics; no global Boolean optimizer.
// The reference must first pass validateBoundaryTree. The builder does NOT set
// diagnostics: callers must run diagnoseBoundarySigns before deploying it.
BoundaryLocalResult buildBoundaryLocalTree(const BoundaryTreeResult&, const BoundaryLocalOptions& = {});
// Call only after accepting closure/sign-safety diagnostics, or for an original
// obstacle fallback. Zero guard remains necessary unless strict_sign_certified.
bool boundaryLocalStrictFree(const BoundaryLocalResult&, const OctreeOccupancy&, const HS::Vec3&,
                             bool* used_zero_guard = nullptr);
void writeBoundaryLocalReport(const std::filesystem::path&, const BoundaryTreeResult&,
                             const BoundarySignDiagnostics& flat, const BoundaryLocalResult* local);
} // namespace rokae_demo
