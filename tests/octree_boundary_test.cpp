#include <rokae_demo/octree_boundary.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>

using namespace rokae_demo;

void check(bool condition, const char* message)
{
  if(!condition) throw std::runtime_error(message);
}

template<class F> void rejects(F function)
{
  bool caught = false;
  try { function(); } catch(const std::exception&) { caught = true; }
  check(caught, "Expected boundary validation to reject mutated data");
}

OctreeBoundary exercise(const std::vector<VoxelPoint>& points, double size = 1)
{
  const auto octree = buildOctreeOccupancy(points, size, 52);
  const auto partition = partitionOctreeCells(octree, pruneOctree(octree), true);
  validateOctreeCellPartition(octree, partition);
  const auto boundary = extractOctreeBoundary(octree, partition);
  validateOctreeBoundary(octree, partition, boundary);
  check(boundary.voxel_faces_before == boundary.faces.size() + boundary.internal_faces_removed,
        "Face accounting failed");
  check(boundary.internal_faces_removed % 2 == 0, "Shared face must remove both directions");
  long double signed_volume3 = 0;
  for(const auto& p : boundary.patches)
  {
    const auto corners = octreeBoundaryCorners(p);
    std::array<long double,3> a{}, b{}, cross{};
    for(unsigned k = 0; k < 3; ++k)
    {
      a[k] = corners[1][k] - corners[0][k];
      b[k] = corners[3][k] - corners[0][k];
    }
    for(unsigned k = 0; k < 3; ++k) cross[k] = a[(k+1)%3]*b[(k+2)%3] - a[(k+2)%3]*b[(k+1)%3];
    check(cross[p.axis]*p.sign == p.source_face_ids.size(), "Incorrect polygon orientation/area");
    check(cross[(p.axis+1)%3] == 0 && cross[(p.axis+2)%3] == 0, "Polygon is not coplanar");
    signed_volume3 += cross[p.axis] * p.plane;
  }
  check(signed_volume3 == 3 * octree.occupied_voxels.size(), "Boundary changed occupied volume");
  // Removing faces must not depend on the selected convex partition.
  const auto fine = partitionOctreeCells(octree, pruneOctree(octree, false), false);
  const auto fine_boundary = extractOctreeBoundary(octree, fine);
  validateOctreeBoundary(octree, fine, fine_boundary);
  check(fine_boundary.patches.size() == boundary.patches.size() && fine_boundary.faces.size() == boundary.faces.size(),
        "Boundary depends on pruning");
  for(std::size_t i = 0; i < boundary.patches.size(); ++i)
    check(octreeBoundaryCorners(boundary.patches[i]) == octreeBoundaryCorners(fine_boundary.patches[i]) &&
          boundary.patches[i].source_face_ids == fine_boundary.patches[i].source_face_ids,
          "Boundary ordering depends on pruning");
  auto reverse = points;
  std::reverse(reverse.begin(), reverse.end());
  const auto reversed = buildOctreeOccupancy(reverse, size, 52);
  const auto reversed_partition = partitionOctreeCells(reversed, pruneOctree(reversed), true);
  const auto other = extractOctreeBoundary(reversed, reversed_partition);
  check(other.patches.size() == boundary.patches.size(), "Input permutation changed boundary");
  for(std::size_t i = 0; i < boundary.patches.size(); ++i)
    check(octreeBoundaryCorners(other.patches[i]) == octreeBoundaryCorners(boundary.patches[i]) &&
          other.patches[i].source_face_ids == boundary.patches[i].source_face_ids, "Unstable boundary patches");

  auto bad = boundary;
  bad.faces[0].sign *= -1;
  rejects([&] { validateOctreeBoundary(octree, partition, bad); });
  bad = boundary; bad.faces[0].convex_cell_id = partition.cells.size();
  rejects([&] { validateOctreeBoundary(octree, partition, bad); });
  bad = boundary; bad.faces.back() = bad.faces.front();
  rejects([&] { validateOctreeBoundary(octree, partition, bad); });
  bad = boundary; bad.patches[0].u1 += 1;
  rejects([&] { validateOctreeBoundary(octree, partition, bad); });
  bad = boundary; bad.patches[0].source_face_ids.push_back(bad.patches[0].source_face_ids.front());
  rejects([&] { validateOctreeBoundary(octree, partition, bad); });
  bad = boundary; bad.patches.pop_back();
  rejects([&] { validateOctreeBoundary(octree, partition, bad); });
  bad = boundary; bad.patches[0].u0 = std::numeric_limits<std::int64_t>::min();
  bad.patches[0].u1 = std::numeric_limits<std::int64_t>::max();
  rejects([&] { validateOctreeBoundary(octree, partition, bad); });
  return boundary;
}

int main()
{
  try
  {
    check(exercise({{-.5,-.5,-.5}}).patches.size() == 6, "Singleton must have six faces");
    for(unsigned mask = 1; mask < 256; ++mask)
    {
      std::vector<VoxelPoint> points;
      for(unsigned bit = 0; bit < 8; ++bit) if(mask & (1u << bit))
        points.emplace_back((bit & 1) - .5, ((bit >> 1) & 1) - .5, ((bit >> 2) & 1) - .5);
      exercise(points);
    }
    std::vector<VoxelPoint> cube, shell, frame;
    for(int x = 0; x < 3; ++x) for(int y = 0; y < 3; ++y) for(int z = 0; z < 3; ++z)
    {
      cube.emplace_back(x+.5,y+.5,z+.5);
      if(x != 1 || y != 1 || z != 1) shell.emplace_back(x+.5,y+.5,z+.5);
      if(z == 0 && (x != 1 || y != 1)) frame.emplace_back(x+.5,y+.5,z+.5);
    }
    const auto solid = exercise(cube);
    check(solid.faces.size() == 54 && solid.patches.size() == 6, "Solid cube should merge to six polygons");
    const auto hollow = exercise(shell);
    check(hollow.faces.size() == 60 && hollow.patches.size() == 12, "Cavity walls lost");
    check(exercise(frame).faces.size() == 32, "Planar hole filled");
    // Edge-only/vertex-only contact must not erase any area of boundary.
    check(exercise({{.5,.5,.5},{1.5,1.5,.5}}).faces.size() == 12, "Edge contact incorrectly culled");
    check(exercise({{.5,.5,.5},{1.5,1.5,1.5}}).faces.size() == 12, "Vertex contact incorrectly culled");
    exercise({{-.025,.015,.015},{.005,.015,.015},{.005,.015,.015}}, .03);
    exercise({{0,0,0},{.03,.03,.03},{-.03,-.03,-.03}}, .03);
    // Sparse, large coordinate gaps must not trigger a dense grid scan.
    check(exercise({{-1000000.5,0,0},{1000000.5,0,0}}).patches.size() == 12, "Disconnected components joined");
    std::mt19937 generator(38142);
    for(unsigned trial = 0; trial < 40; ++trial)
    {
      std::vector<VoxelPoint> points{{.5,.5,.5}};
      for(int x = -2; x <= 2; ++x) for(int y = -2; y <= 2; ++y) for(int z = -2; z <= 2; ++z)
        if(generator()%3 == 0) points.emplace_back(x+.5,y+.5,z+.5);
      exercise(points);
    }
    std::cout << "Octree boundary exact coverage, holes, orientation, determinism and mutation tests passed\n";
  }
  catch(const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
