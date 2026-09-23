#pragma once
#include <rokae_demo/boundary_direct_tree.hpp>
#include <rokae_demo/boundary_distance_reference.hpp>
#include <rokae_demo/halfspace_smooth_tree.hpp>
#include <rokae_demo/halfspace_validation_cache.hpp>

namespace rokae_demo
{
struct SmoothDistanceAuditOptions
{
  std::size_t max_work=100000000, gradient_samples=128, zero_rays=64;
  bool share_subexpressions=false;
  bool cache_probes=false; // validation only; preserve all occurrences and work limits
};
struct SmoothDistanceProbe
{
  HS::Vec3 point,gradient;
  double hard=0,smooth=0,reference=0;
};
struct SmoothZeroRay
{
  std::size_t patch=0; bool bracketed=false,blocked=false;
  double offset=0,reference_distance=0;
};
struct SmoothDistanceAudit
{
  std::map<std::string,double> metrics;
  std::vector<SmoothDistanceProbe> probes;
  std::vector<SmoothZeroRay> rays;
};
// Requires the existing exact occupancy and boundary-direct certificates first.
// The reference geometry NEVER enters the smooth field or gradient calculation.
SmoothDistanceAudit auditSmoothDistance(const SparseVoxelOccupancy&,const OctreeBoundary&,
  const std::vector<VoxelPoint>& raw,const std::vector<HS::Vec3>& random_points,const HS::SmoothTree&,
  const BoundaryDistanceReference&,const SmoothDistanceAuditOptions& = {},HS::ValidationValueCache* = nullptr);
// Compatibility overload: independently compile and audit the supplied tree,
// retaining the original contract if it differs from the smooth snapshot.
SmoothDistanceAudit auditSmoothDistance(const SparseVoxelOccupancy&,const OctreeBoundary&,const TreeRepresentation&,
  const std::vector<VoxelPoint>& raw,const std::vector<HS::Vec3>& random_points,const HS::SmoothTree&,
  const BoundaryDistanceReference&,const SmoothDistanceAuditOptions& = {});
void writeSmoothDistance(const std::filesystem::path&,const HS::SmoothTree&,const SmoothDistanceAudit&);
} // namespace rokae_demo
