// Same-process benchmark against the previous ordered-set implementation.
// Also verifies identical occupied cells, corner coordinates/order, and checks.
#include <rokae_demo/conservative_voxelization.hpp>

#include <CGAL/IO/read_points.h>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <set>
#include <stdexcept>

namespace
{
using namespace rokae_demo;
using Clock = std::chrono::steady_clock;

double milliseconds(Clock::time_point start)
{
  return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

ConservativeVoxelization orderedSetReference(const std::vector<VoxelPoint>& raw,
                                            double size)
{
  const auto start = Clock::now();
  ConservativeVoxelization result;
  result.voxel_size = size;
  result.stats.voxel_size = size;
  result.stats.raw_points = raw.size();
  std::set<VoxelIndex> occupied;
  for(const auto& point : raw)
    occupied.insert(voxelIndex(point, size));
  result.occupied_voxels.assign(occupied.begin(), occupied.end());
  std::set<VoxelIndex> corners;
  for(const auto& cell : occupied)
    for(int x = 0; x < 2; ++x)
      for(int y = 0; y < 2; ++y)
        for(int z = 0; z < 2; ++z)
          corners.insert({cell[0] + x, cell[1] + y, cell[2] + z});
  result.corner_indices.assign(corners.begin(), corners.end());
  result.points.reserve(corners.size());
  for(const auto& index : corners)
    result.points.emplace_back(index[0] * size, index[1] * size, index[2] * size);
  result.stats.occupied_voxels = occupied.size();
  result.stats.voxel_points = corners.size();
  result.stats.voxel_point_ratio = static_cast<double>(corners.size()) / raw.size();
  result.stats.construction_ms = milliseconds(start);
  const auto validation_start = Clock::now();
  result.stats.raw_points_outside_voxels = result.countOutside(raw);
  result.stats.validation_ms = milliseconds(validation_start);
  result.stats.total_ms = milliseconds(start);
  return result;
}

double median(std::vector<double> values)
{
  std::sort(values.begin(), values.end());
  return (values[(values.size() - 1) / 2] + values[values.size() / 2]) / 2.0;
}

void benchmark(const std::vector<VoxelPoint>& raw, double size)
{
  std::vector<double> old_times, new_times, old_total, new_total;
  ConservativeVoxelization reference, current;
  for(int repeat = 0; repeat < 33; ++repeat)
  {
    // Alternate execution order, and exclude three warm-up iterations.
    if(repeat % 2 == 0)
    {
      reference = orderedSetReference(raw, size);
      current = conservativeVoxelize(raw, size);
    }
    else
    {
      current = conservativeVoxelize(raw, size);
      reference = orderedSetReference(raw, size);
    }
    if(reference.occupied_voxels != current.occupied_voxels ||
       reference.corner_indices != current.corner_indices || reference.points != current.points ||
       reference.stats.raw_points_outside_voxels != 0 || current.stats.raw_points_outside_voxels != 0)
      throw std::runtime_error("Voxel optimization changed geometry or containment");
    if(repeat >= 3)
    {
      old_times.push_back(reference.stats.construction_ms);
      new_times.push_back(current.stats.construction_ms);
      old_total.push_back(reference.stats.total_ms);
      new_total.push_back(current.stats.total_ms);
    }
  }
  std::cout << size << ',' << raw.size() << ',' << current.occupied_voxels.size()
            << ',' << current.points.size() << ',' << median(old_times) << ','
            << median(new_times) << ',' << median(old_total) << ',' << median(new_total)
            << ',' << median(old_times) / median(new_times) << ",true,0\n";
}
}  // namespace

int main(int argc, char** argv)
{
  if(argc != 2)
  {
    std::cerr << "Usage: voxelization_benchmark INPUT.xyz\n";
    return 1;
  }
  try
  {
    std::vector<VoxelPoint> raw;
    if(!CGAL::IO::read_points(argv[1], std::back_inserter(raw)) || raw.empty())
      throw std::runtime_error("Cannot read non-empty benchmark point set");
    std::cout << "voxel_size,raw_points,occupied_voxels,voxel_points,set_construction_ms,"
                 "vector_construction_ms,set_total_ms,vector_total_ms,speedup,"
                 "identical_output,raw_points_outside_voxels\n";
    for(double size : {0.01, 0.02, 0.03})
      benchmark(raw, size);
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
