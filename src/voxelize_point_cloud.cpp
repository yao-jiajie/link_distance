#include <rokae_demo/conservative_voxelization.hpp>
#include <rokae_demo/units.hpp>

#include <CGAL/IO/read_points.h>
#include <CGAL/IO/write_points.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv)
{
  if(argc != 5 || std::string(argv[3]) != "--voxel-size")
  {
    std::cerr << "Usage: voxelize_point_cloud INPUT.xyz OUTPUT_DIR --voxel-size METERS\n";
    return 1;
  }
  try
  {
    std::size_t consumed = 0;
    const std::string argument = argv[4];
    const double size = std::stod(argument, &consumed);
    if(consumed != argument.size() || !std::isfinite(size) || size <= 0.0)
      throw std::invalid_argument("voxel_size must be finite and positive (meters)");
    std::vector<rokae_demo::VoxelPoint> raw;
    if(!CGAL::IO::read_points(argv[1], std::back_inserter(raw)))
      throw std::runtime_error("Cannot read point cloud");
    const auto result = rokae_demo::conservativeVoxelize(raw, size);
    const auto& stats = result.stats;
    const std::filesystem::path output(argv[2]);
    std::filesystem::create_directories(output);
    if(!CGAL::IO::write_points((output / "voxel_corners.xyz").string(), result.points,
                               CGAL::parameters::stream_precision(17)))
      throw std::runtime_error("Cannot write voxel corners");
    std::ofstream occupied(output / "occupied_voxels.csv");
    std::ofstream csv(output / "benchmark.csv");
    std::ofstream json(output / "benchmark.json");
    occupied.exceptions(std::ios::failbit | std::ios::badbit);
    csv.exceptions(std::ios::failbit | std::ios::badbit);
    json.exceptions(std::ios::failbit | std::ios::badbit);
    occupied << "ix,iy,iz\n";
    for(const auto& index : result.occupied_voxels)
      occupied << index[0] << ',' << index[1] << ',' << index[2] << '\n';
    csv << "length_unit,voxel_size,raw_points,occupied_voxels,voxel_points,"
           "voxel_point_ratio,raw_points_outside_voxels,voxelization_ms,validation_ms,total_ms\n"
        << std::setprecision(17) << rokae_demo::kLengthUnit << ',' << size << ','
        << stats.raw_points << ',' << stats.occupied_voxels << ',' << stats.voxel_points << ','
        << stats.voxel_point_ratio << ',' << stats.raw_points_outside_voxels << ','
        << stats.construction_ms << ',' << stats.validation_ms << ',' << stats.total_ms << '\n';
    json << std::setprecision(17)
         << "{\n  \"stage\": \"voxelization\",\n  \"length_unit\": \"m\",\n"
         << "  \"grid_origin\": [0, 0, 0],\n  \"voxel_size\": " << size
         << ",\n  \"raw_points\": " << stats.raw_points
         << ",\n  \"occupied_voxels\": " << stats.occupied_voxels
         << ",\n  \"voxel_points\": " << stats.voxel_points
         << ",\n  \"voxel_point_ratio\": " << stats.voxel_point_ratio
         << ",\n  \"raw_points_outside_voxels\": " << stats.raw_points_outside_voxels
         << ",\n  \"voxelization_ms\": " << stats.construction_ms
         << ",\n  \"validation_ms\": " << stats.validation_ms
         << ",\n  \"total_ms\": " << stats.total_ms << "\n}\n";
    occupied.close();
    csv.close();
    json.close();
    std::cout << "voxel_size=" << size << " raw_points=" << stats.raw_points
              << " occupied_voxels=" << stats.occupied_voxels
              << " voxel_points=" << stats.voxel_points
              << " voxel_point_ratio=" << stats.voxel_point_ratio
              << " raw_points_outside_voxels=" << stats.raw_points_outside_voxels
              << " total_ms=" << stats.total_ms << '\n';
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr << "Voxelization failed: " << error.what() << '\n';
    return 1;
  }
}
