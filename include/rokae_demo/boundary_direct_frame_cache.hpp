#pragma once

#include <rokae_demo/boundary_smooth_distance.hpp>
#include <cstddef>
#include <list>
#include <memory>
#include <optional>
#include <vector>

namespace rokae_demo
{
// Configuration that changes a reusable geometry/execution snapshot. Values
// are compared exactly; changing a work limit also forces recertification.
struct BoundaryDirectFrameOptions
{
  BoundaryDirectOptions direct;
  bool share_subexpressions = true;
  bool smooth = false;
  HS::SmoothTreeOptions smoothing;
  bool cache_validation_values = false;
};

// Immutable geometry and program state plus caller-owned reusable workspaces.
// One snapshot is intended for one session/worker. Share compiledProgram()
// across threads only when each thread supplies separate workspaces.
class BoundaryDirectFrameSnapshot
{
  friend class BoundaryDirectFrameCache;
  SparseVoxelOccupancy occupancy_;
  OctreeBoundary boundary_;
  BoundaryDirectResult direct_;
  BoundarySignDiagnostics diagnostics_;
  std::shared_ptr<const HS::CompiledTree> program_;
  std::optional<HS::SmoothTree> smooth_;
  std::optional<BoundaryDistanceReference> distance_reference_;
  std::optional<HS::ValidationValueCache> validation_cache_;
  std::vector<double> query_workspace_,smooth_workspace_;
  std::vector<HS::ValueGradient> gradient_workspace_;
  BoundaryDirectFrameOptions options_;
  std::size_t generation_ = 0;

public:
  const SparseVoxelOccupancy& occupancy() const { return occupancy_; }
  const OctreeBoundary& boundary() const { return boundary_; }
  const BoundaryDirectResult& direct() const { return direct_; }
  const BoundarySignDiagnostics& diagnostics() const { return diagnostics_; }
  const std::shared_ptr<const HS::CompiledTree>& compiledProgram() const { return program_; }
  const HS::SmoothTree* smoothTree() const { return smooth_ ? &*smooth_ : nullptr; }
  const BoundaryDistanceReference* distanceReference() const
  { return distance_reference_ ? &*distance_reference_ : nullptr; }
  HS::ValidationValueCache* validationCache()
  { return validation_cache_ ? &*validation_cache_ : nullptr; }
  const BoundaryDirectFrameOptions& options() const { return options_; }
  std::size_t generation() const { return generation_; }

  std::vector<double>& queryWorkspace() { return query_workspace_; }
  std::vector<double>& smoothWorkspace() { return smooth_workspace_; }
  std::vector<HS::ValueGradient>& gradientWorkspace() { return gradient_workspace_; }
};

struct BoundaryDirectPreparedFrame
{
  std::shared_ptr<BoundaryDirectFrameSnapshot> snapshot;
  bool cache_hit = false;
  double occupancy_ms = 0;
  double lookup_ms = 0;
  double raw_validation_ms = 0; // cache-hit raw-point containment only
  double boundary_ms = 0;
  double direct_ms = 0;
  double snapshot_build_ms = 0;
};

// Exact-key LRU for a sequence of frames. prepare() always voxelizes the new
// raw frame and checks every raw point against the selected closed occupancy.
// A retained shared_ptr keeps an evicted snapshot and its program alive.
class BoundaryDirectFrameCache
{
  std::size_t capacity_,hits_ = 0,misses_ = 0,evictions_ = 0,generation_ = 0;
  std::list<std::shared_ptr<BoundaryDirectFrameSnapshot>> entries_;

public:
  explicit BoundaryDirectFrameCache(std::size_t capacity = 1);
  BoundaryDirectPreparedFrame prepare(const std::vector<VoxelPoint>& raw,double voxel_size,
                                      const BoundaryDirectFrameOptions& = {});
  void clear() { entries_.clear(); }
  std::size_t capacity() const { return capacity_; }
  std::size_t size() const { return entries_.size(); }
  std::size_t hits() const { return hits_; }
  std::size_t misses() const { return misses_; }
  std::size_t evictions() const { return evictions_; }
};
} // namespace rokae_demo
