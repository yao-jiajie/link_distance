#include <rokae_demo/boundary_direct_frame_cache.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>

namespace rokae_demo
{
namespace
{
using Clock=std::chrono::steady_clock;
double ms(Clock::time_point start)
{ return std::chrono::duration<double,std::milli>(Clock::now()-start).count(); }

bool sameDouble(double a,double b)
{
  static_assert(sizeof(double)==sizeof(std::uint64_t));
  std::uint64_t x=0,y=0;
  std::memcpy(&x,&a,sizeof(a)); std::memcpy(&y,&b,sizeof(b));
  return x==y;
}

bool sameOptions(const BoundaryDirectFrameOptions& a,const BoundaryDirectFrameOptions& b)
{
  const auto& x=a.direct; const auto& y=b.direct;
  const auto& xs=x.signs; const auto& ys=y.signs;
  return x.max_canonical_cells==y.max_canonical_cells && x.max_split_checks==y.max_split_checks &&
    x.max_nodes==y.max_nodes && x.max_expanded_nodes==y.max_expanded_nodes && x.max_depth==y.max_depth &&
    x.intern_subexpressions==y.intern_subexpressions &&
    x.materialize_expression==y.materialize_expression &&
    x.retain_partition_diagnostics==y.retain_partition_diagnostics &&
    xs.max_expanded_nodes==ys.max_expanded_nodes && xs.max_attempts==ys.max_attempts &&
    xs.max_adjacency_checks==ys.max_adjacency_checks && xs.max_strata==ys.max_strata &&
    xs.max_sign_operations==ys.max_sign_operations && xs.packed_sign_propagation==ys.packed_sign_propagation &&
    a.share_subexpressions==b.share_subexpressions && a.smooth==b.smooth &&
    sameDouble(a.smoothing.beta,b.smoothing.beta) && sameDouble(a.smoothing.max_error,b.smoothing.max_error) &&
    a.smoothing.share_subexpressions==b.smoothing.share_subexpressions &&
    a.smoothing.binary_kernel==b.smoothing.binary_kernel &&
    a.smoothing.predecode_kernel==b.smoothing.predecode_kernel &&
    a.smoothing.paired_batch==b.smoothing.paired_batch &&
    a.cache_validation_values==b.cache_validation_values;
}

bool sameGeometry(const SparseVoxelOccupancy& occupancy,double h,const std::vector<VoxelIndex>& keys)
{ return sameDouble(occupancy.voxel_size,h) && occupancy.occupied_voxels==keys; }

void validateRawPoints(const SparseVoxelOccupancy& occupancy,const std::vector<VoxelPoint>& raw)
{
  for(const auto& point:raw)
    if(!occupancy.contains({point.x(),point.y(),point.z()}))
      throw std::runtime_error("Raw point outside reusable closed voxel union");
}
} // namespace

BoundaryDirectFrameCache::BoundaryDirectFrameCache(std::size_t capacity):capacity_(capacity)
{
  if(!capacity_) throw std::invalid_argument("Frame cache capacity must be positive");
}

BoundaryDirectPreparedFrame BoundaryDirectFrameCache::prepare(const std::vector<VoxelPoint>& raw,double voxel_size,
                                                              const BoundaryDirectFrameOptions& options)
{
  if(!std::isfinite(voxel_size) || voxel_size<=0)
    throw std::invalid_argument("Frame cache requires a positive finite voxel size");
  if(options.smooth && options.smoothing.share_subexpressions!=options.share_subexpressions)
    throw std::invalid_argument("Frame smoothing and compiled sharing options must match");

  BoundaryDirectPreparedFrame prepared; auto start=Clock::now();
  auto occupancy=buildSparseVoxelOccupancy(raw,voxel_size); prepared.occupancy_ms=ms(start);
  start=Clock::now();
  const auto found=std::find_if(entries_.begin(),entries_.end(),[&](const auto& entry) {
    return sameOptions(entry->options_,options) && sameGeometry(entry->occupancy_,voxel_size,occupancy.occupied_voxels);
  });
  prepared.lookup_ms=ms(start);
  if(found!=entries_.end())
  {
    prepared.snapshot=*found; entries_.splice(entries_.begin(),entries_,found);
    start=Clock::now();
    validateRawPoints(prepared.snapshot->occupancy_,raw);
    prepared.raw_validation_ms=ms(start);
    prepared.cache_hit=true; ++hits_; return prepared;
  }

  ++misses_; start=Clock::now();
  auto snapshot=std::make_shared<BoundaryDirectFrameSnapshot>();
  snapshot->occupancy_=std::move(occupancy); snapshot->options_=options; snapshot->generation_=++generation_;
  validateSparseVoxelOccupancy(snapshot->occupancy_,raw);
  const auto boundary_start=Clock::now();
  snapshot->boundary_=extractSparseVoxelBoundary(snapshot->occupancy_);
  prepared.boundary_ms=ms(boundary_start);
  const auto direct_start=Clock::now();
  snapshot->direct_=buildBoundaryDirectTree(snapshot->occupancy_,snapshot->boundary_,options.direct);
  prepared.direct_ms=ms(direct_start);
  snapshot->diagnostics_=validateBoundaryDirectTree(snapshot->occupancy_,snapshot->boundary_,snapshot->direct_,options.direct);
  if(snapshot->diagnostics_.applied && snapshot->diagnostics_.strict_sign_certified)
  {
    if(options.share_subexpressions && options.direct.intern_subexpressions &&
       snapshot->direct_.construction_interned)
      snapshot->program_=std::make_shared<const HS::CompiledTree>(snapshot->direct_.serialized,
        snapshot->direct_.tree.leaves,HS::preinterned_program);
    else
      snapshot->program_=std::make_shared<const HS::CompiledTree>(snapshot->direct_.serialized,
        snapshot->direct_.tree.leaves,options.share_subexpressions);
    snapshot->query_workspace_=snapshot->program_->makeWorkspace();
    if(options.smooth)
    {
      snapshot->smooth_.emplace(snapshot->program_,options.smoothing);
      snapshot->distance_reference_.emplace(voxel_size,snapshot->boundary_);
      snapshot->smooth_workspace_=snapshot->smooth_->makeValueWorkspace();
      snapshot->gradient_workspace_=snapshot->smooth_->makeGradientWorkspace();
    }
    if(options.cache_validation_values)
      snapshot->validation_cache_.emplace(snapshot->program_);
  }
  prepared.snapshot=snapshot; prepared.snapshot_build_ms=ms(start);
  if(!snapshot->diagnostics_.applied || !snapshot->diagnostics_.strict_sign_certified)
    return prepared;

  if(entries_.size()==capacity_)
  { entries_.pop_back(); ++evictions_; }
  entries_.push_front(std::move(snapshot));
  return prepared;
}
} // namespace rokae_demo
