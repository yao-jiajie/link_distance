#include <rokae_demo/conservative_voxelization.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>

namespace rokae_demo
{
namespace
{

using Clock = std::chrono::steady_clock;

double milliseconds(Clock::time_point start)
{
  return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

void checkSize(double size)
{
  if(!std::isfinite(size) || size <= 0.0)
    throw std::invalid_argument("voxel_size must be finite and positive (meters)");
}

double gridPlane(std::int64_t index, double size)
{
  const double value = static_cast<double>(index) * size;
  if(!std::isfinite(value))
    throw std::out_of_range("Voxel corner is not representable in double");
  return value;
}

std::int64_t axisIndex(double value, double size)
{
  if(!std::isfinite(value))
    throw std::invalid_argument("Point coordinates must be finite");
  const long double quotient = static_cast<long double>(value) / size;
  // Leave precision headroom for distinct adjacent double grid planes and for
  // index + 1. Never convert an unchecked floating-point value to an integer.
  constexpr std::int64_t limit = std::int64_t{1} << 50;
  if(!std::isfinite(quotient) || quotient <= -limit || quotient >= limit)
    throw std::out_of_range("Coordinate/voxel_size exceeds supported grid range");
  auto index = static_cast<std::int64_t>(std::floor(quotient));
  // Quotient and product rounding may disagree at a decimal grid boundary.
  // Assign using the actual exported planes so no epsilon can hide a missed point.
  for(int attempt = 0; attempt < 4; ++attempt)
  {
    const double lower = gridPlane(index, size);
    const double upper = gridPlane(index + 1, size);
    if(!(lower < upper))
      throw std::out_of_range("Adjacent voxel planes are not distinguishable");
    if(value < lower)
      --index;
    else if(value >= upper)
      ++index;
    else
      return index;
  }
  throw std::out_of_range("Cannot consistently index point on voxel grid");
}

VoxelPoint corner(const VoxelIndex& index, double size)
{
  return {gridPlane(index[0], size), gridPlane(index[1], size),
          gridPlane(index[2], size)};
}

}  // namespace

VoxelIndex voxelIndex(const VoxelPoint& point, double voxel_size)
{
  checkSize(voxel_size);
  return {axisIndex(point.x(), voxel_size), axisIndex(point.y(), voxel_size),
          axisIndex(point.z(), voxel_size)};
}

bool ConservativeVoxelization::contains(const VoxelPoint& point) const
{
  const VoxelIndex primary = voxelIndex(point, voxel_size);
  const std::array<double, 3> coordinates{point.x(), point.y(), point.z()};
  // A grid-plane point can also belong to the cell below it, since the
  // represented occupied boxes are closed even though indexing is half-open.
  for(unsigned int mask = 0; mask < 8; ++mask)
  {
    VoxelIndex candidate = primary;
    bool within = true;
    for(std::size_t axis = 0; axis < 3; ++axis)
    {
      if((mask & (1u << axis)) != 0)
        --candidate[axis];
      within &= coordinates[axis] >= gridPlane(candidate[axis], voxel_size) &&
                coordinates[axis] <= gridPlane(candidate[axis] + 1, voxel_size);
    }
    if(within && std::binary_search(occupied_voxels.begin(), occupied_voxels.end(),
                                   candidate))
      return true;
  }
  return false;
}

std::size_t ConservativeVoxelization::countOutside(
    const std::vector<VoxelPoint>& raw_points) const
{
  return static_cast<std::size_t>(std::count_if(
      raw_points.begin(), raw_points.end(),
      [&](const VoxelPoint& point) { return !contains(point); }));
}

ConservativeVoxelization conservativeVoxelize(
    const std::vector<VoxelPoint>& raw_points, double voxel_size)
{
  checkSize(voxel_size);
  if(raw_points.empty())
    throw std::invalid_argument("Cannot voxelize an empty point cloud");
  const auto start = Clock::now();
  ConservativeVoxelization result;
  result.voxel_size = voxel_size;
  result.stats.voxel_size = voxel_size;
  result.stats.raw_points = raw_points.size();

  // Contiguous storage avoids one tree-node allocation per voxel/corner.
  // Sorting preserves the same lexicographic order as the original std::set.
  auto& occupied = result.occupied_voxels;
  occupied.reserve(raw_points.size());
  for(const VoxelPoint& point : raw_points)
    occupied.push_back(voxelIndex(point, voxel_size));
  std::sort(occupied.begin(), occupied.end());
  occupied.erase(std::unique(occupied.begin(), occupied.end()), occupied.end());

  auto& corners = result.corner_indices;
  if(occupied.size() > corners.max_size() / 8)
    throw std::length_error("Too many occupied voxels to enumerate corners");
  corners.reserve(8 * occupied.size());
  for(const VoxelIndex& voxel : occupied)
  {
    for(int x = 0; x < 2; ++x)
      for(int y = 0; y < 2; ++y)
        for(int z = 0; z < 2; ++z)
          corners.push_back({voxel[0] + x, voxel[1] + y, voxel[2] + z});
  }
  std::sort(corners.begin(), corners.end());
  corners.erase(std::unique(corners.begin(), corners.end()), corners.end());
  result.points.reserve(corners.size());
  for(const VoxelIndex& index : corners)
    result.points.push_back(corner(index, voxel_size));
  result.stats.occupied_voxels = result.occupied_voxels.size();
  result.stats.voxel_points = result.points.size();
  result.stats.voxel_point_ratio =
      static_cast<double>(result.points.size()) / raw_points.size();
  result.stats.construction_ms = milliseconds(start);

  const auto validation_start = Clock::now();
  result.stats.raw_points_outside_voxels = result.countOutside(raw_points);
  result.stats.validation_ms = milliseconds(validation_start);
  result.stats.total_ms = milliseconds(start);
  if(result.stats.raw_points_outside_voxels != 0)
    throw std::runtime_error("Raw points lie outside occupied voxel boxes");
  return result;
}

}  // namespace rokae_demo
