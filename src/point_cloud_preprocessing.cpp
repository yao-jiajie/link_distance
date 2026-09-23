#include <rokae_demo/point_cloud_preprocessing.hpp>

#include <CGAL/compute_average_spacing.h>
#include <CGAL/estimate_scale.h>
#include <CGAL/jet_smooth_point_set.h>
#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/tags.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace rokae_demo
{
namespace
{

double elapsedMilliseconds(const std::chrono::steady_clock::time_point& start)
{
  return std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - start)
      .count();
}

void computeDisplacementStatistics(
    const std::vector<PointCloudPoint>& input,
    const std::vector<PointCloudPoint>& output,
    PointCloudPreprocessStats& stats)
{
  std::vector<double> displacements;
  displacements.reserve(input.size());
  double squared_sum = 0.0;
  for(std::size_t i = 0; i < input.size(); ++i)
  {
    const double displacement = std::sqrt(CGAL::to_double(
        CGAL::squared_distance(input[i], output[i])));
    if(!std::isfinite(displacement))
      throw std::runtime_error("Point smoothing produced a non-finite point");
    displacements.push_back(displacement);
    squared_sum += displacement * displacement;
  }

  stats.mean_displacement =
      std::accumulate(displacements.begin(), displacements.end(), 0.0) /
      static_cast<double>(displacements.size());
  stats.rms_displacement =
      std::sqrt(squared_sum / static_cast<double>(displacements.size()));
  stats.maximum_displacement =
      *std::max_element(displacements.begin(), displacements.end());
  std::sort(displacements.begin(), displacements.end());
  const std::size_t percentile_index = static_cast<std::size_t>(
      std::ceil(0.95 * static_cast<double>(displacements.size())) - 1.0);
  stats.percentile95_displacement = displacements[percentile_index];
}

}  // namespace

const char* pointCloudPreprocessModeName(PointCloudPreprocessMode mode)
{
  switch(mode)
  {
    case PointCloudPreprocessMode::NONE:
      return "none";
    case PointCloudPreprocessMode::JET:
      return "jet";
    case PointCloudPreprocessMode::VOXEL:
      return "voxel";
  }
  return "unknown";
}

PointCloudPreprocessResult preprocessPointCloud(
    const std::vector<PointCloudPoint>& input,
    const PointCloudPreprocessOptions& options)
{
  if(input.empty())
    throw std::runtime_error("Cannot preprocess an empty point cloud");
  if(options.iterations == 0)
    throw std::runtime_error("Smoothing iterations must be positive");

  const auto total_start = std::chrono::steady_clock::now();
  PointCloudPreprocessResult result;
  result.stats.mode = options.mode;
  result.stats.input_points = input.size();
  result.stats.output_points = input.size();
  result.stats.requested_neighbors = options.neighbors;
  result.stats.iterations =
      options.mode == PointCloudPreprocessMode::JET ? options.iterations : 0;

  if(options.mode == PointCloudPreprocessMode::VOXEL)
  {
    auto voxels = conservativeVoxelize(input, options.voxel_size);
    result.points = std::move(voxels.points);
    result.stats.voxelization = voxels.stats;
    result.stats.output_points = result.points.size();
    result.stats.total_ms = elapsedMilliseconds(total_start);
    return result;
  }
  result.points = input;

  if(options.mode == PointCloudPreprocessMode::NONE)
  {
    result.stats.total_ms = elapsedMilliseconds(total_start);
    return result;
  }

  // A quadratic jet in two parameters has six coefficients. Fewer samples
  // cannot define CGAL's default degree-2 fit reliably.
  if(input.size() < 6)
    throw std::runtime_error(
        "Jet smoothing requires at least six input points");

  auto start = std::chrono::steady_clock::now();
  const unsigned int spacing_neighbors = static_cast<unsigned int>(
      std::min<std::size_t>(6, input.size()));
  result.stats.average_spacing = CGAL::to_double(
      CGAL::compute_average_spacing<CGAL::Parallel_if_available_tag>(
          input, spacing_neighbors));

  std::size_t effective_neighbors = options.neighbors;
  if(effective_neighbors == 0)
    effective_neighbors = CGAL::estimate_global_k_neighbor_scale(input);
  effective_neighbors =
      std::max<std::size_t>(6, std::min(effective_neighbors, input.size()));
  result.stats.effective_neighbors = effective_neighbors;
  result.stats.scale_estimation_ms = elapsedMilliseconds(start);

  start = std::chrono::steady_clock::now();
  for(std::size_t iteration = 0; iteration < options.iterations; ++iteration)
  {
    CGAL::jet_smooth_point_set<CGAL::Parallel_if_available_tag>(
        result.points, static_cast<unsigned int>(effective_neighbors));
  }
  result.stats.smoothing_ms = elapsedMilliseconds(start);
  computeDisplacementStatistics(input, result.points, result.stats);
  result.stats.total_ms = elapsedMilliseconds(total_start);
  return result;
}

}  // namespace rokae_demo
