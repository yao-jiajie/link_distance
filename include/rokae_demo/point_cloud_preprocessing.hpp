#pragma once

#include <rokae_demo/conservative_voxelization.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <cstddef>
#include <string>
#include <vector>

namespace rokae_demo
{

using PointCloudKernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using PointCloudPoint = PointCloudKernel::Point_3;

enum class PointCloudPreprocessMode
{
  NONE,
  JET,
  VOXEL
};

struct PointCloudPreprocessOptions
{
  PointCloudPreprocessMode mode = PointCloudPreprocessMode::NONE;
  // Zero selects CGAL's global K-neighbor scale estimator.
  std::size_t neighbors = 0;
  std::size_t iterations = 1;
  double voxel_size = 0.0;
};

struct PointCloudPreprocessStats
{
  PointCloudPreprocessMode mode = PointCloudPreprocessMode::NONE;
  std::size_t input_points = 0;
  std::size_t output_points = 0;
  std::size_t requested_neighbors = 0;
  std::size_t effective_neighbors = 0;
  std::size_t iterations = 0;
  double average_spacing = 0.0;
  double scale_estimation_ms = 0.0;
  double smoothing_ms = 0.0;
  double total_ms = 0.0;
  double mean_displacement = 0.0;
  double rms_displacement = 0.0;
  double percentile95_displacement = 0.0;
  double maximum_displacement = 0.0;
  VoxelizationStats voxelization;
};

struct PointCloudPreprocessResult
{
  std::vector<PointCloudPoint> points;
  PointCloudPreprocessStats stats;
};

const char* pointCloudPreprocessModeName(PointCloudPreprocessMode mode);

PointCloudPreprocessResult preprocessPointCloud(
    const std::vector<PointCloudPoint>& input,
    const PointCloudPreprocessOptions& options);

}  // namespace rokae_demo
