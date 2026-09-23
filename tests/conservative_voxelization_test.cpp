#include <rokae_demo/conservative_voxelization.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>

namespace
{
using namespace rokae_demo;

void require(bool condition, const char* message)
{
  // Do not use assert: the project is built in Release mode.
  if(!condition)
    throw std::runtime_error(message);
}

template<class Function> void rejects(Function function)
{
  bool rejected = false;
  try { function(); }
  catch(const std::exception&) { rejected = true; }
  require(rejected, "Expected invalid geometry to be rejected");
}

void runTests()
{
  auto single = conservativeVoxelize({{0.1, 0.2, 0.3}, {0.9, 0.8, 0.7}}, 1.0);
  require(single.occupied_voxels.size() == 1 && single.points.size() == 8,
          "One occupied voxel must yield eight corners");
  require(single.contains({0, 0, 0}) && single.contains({1, 1, 1}),
          "Both ends of a closed voxel must be contained");
  require(!single.contains({1.01, 0.5, 0.5}), "Outside query must be rejected");
  auto adjacent = conservativeVoxelize({{0.5, 0.5, 0.5}, {1.5, 0.5, 0.5}}, 1.0);
  require(adjacent.points.size() == 12, "Adjacent voxels must share four corners");
  auto negative = conservativeVoxelize({{-0.1, -1.0, -1.1}}, 1.0);
  require(negative.occupied_voxels.front() == VoxelIndex{-1, -1, -2},
          "Negative coordinates must use floor, not truncation");
  require(negative.contains({0, 0, -1}), "Closed negative upper corner is inside");

  const double inf = std::numeric_limits<double>::infinity();
  for(double size : {0.01, 0.03, 1.0})
  {
    for(std::int64_t i : {-100, -3, -1, 0, 1, 3, 100})
    {
      const double boundary = static_cast<double>(i) * size;
      require(voxelIndex({boundary, 0, 0}, size)[0] == i,
              "Exported boundary must index to upper voxel");
      require(voxelIndex({std::nextafter(boundary, -inf), 0, 0}, size)[0] == i - 1,
              "Point just below grid plane must index to lower voxel");
      require(voxelIndex({std::nextafter(boundary, inf), 0, 0}, size)[0] == i,
              "Point just above grid plane must index to upper voxel");
    }
  }

  std::mt19937_64 rng(17);
  std::uniform_real_distribution<double> random(-0.5, 0.5);
  std::vector<VoxelPoint> raw;
  for(int i = 0; i < 1000; ++i)
    raw.emplace_back(random(rng), random(rng), random(rng));
  const auto first = conservativeVoxelize(raw, 0.03);
  require(first.countOutside(raw) == 0, "Raw cloud must be fully contained");
  // Independently check actual corner boxes instead of reusing index lookup.
  for(const VoxelPoint& point : raw)
  {
    bool inside = false;
    for(const VoxelIndex& box : first.occupied_voxels)
    {
      const std::array<double, 3> p{point.x(), point.y(), point.z()};
      bool in_box = true;
      for(std::size_t axis = 0; axis < 3; ++axis)
        in_box &= p[axis] >= box[axis] * 0.03 && p[axis] <= (box[axis] + 1) * 0.03;
      inside |= in_box;
    }
    require(inside, "Exported box union must contain each raw coordinate");
  }
  std::shuffle(raw.begin(), raw.end(), rng);
  raw.push_back(raw.front());
  const auto second = conservativeVoxelize(raw, 0.03);
  require(first.occupied_voxels == second.occupied_voxels &&
          first.corner_indices == second.corner_indices && first.points == second.points,
          "Shuffling or duplicating input must not change output geometry");
  single.occupied_voxels.clear();
  require(single.countOutside({{0.1, 0.2, 0.3}}) == 1,
          "Containment check must detect a missing occupied cell");
  rejects([] { conservativeVoxelize({}, 1.0); });
  for(double size : {0.0, -1.0, inf, std::numeric_limits<double>::quiet_NaN()})
    rejects([&] { conservativeVoxelize({{0, 0, 0}}, size); });
  rejects([&] { conservativeVoxelize({{inf, 0, 0}}, 1.0); });
  rejects([] { conservativeVoxelize({{1e30, 0, 0}}, 0.01); });
  rejects([] { conservativeVoxelize({{1, 0, 0}}, 1e-300); });
}
}  // namespace

int main()
{
  try
  {
    runTests();
    std::cout << "Conservative voxelization tests passed\n";
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
