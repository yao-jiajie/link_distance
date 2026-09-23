#include <rokae_demo/octree_rectangular_pruning.hpp>

#include <algorithm>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>

using namespace rokae_demo;

void check(bool value, const char* message)
{
  if(!value) throw std::runtime_error(message);
}

template<class F> void rejects(F function)
{
  bool caught = false;
  try { function(); } catch(const std::exception&) { caught = true; }
  check(caught, "Expected rectangular validation failure");
}

unsigned population(unsigned mask)
{
  unsigned count = 0;
  for(; mask; mask >>= 1) count += mask & 1u;
  return count;
}

// Independent enumeration of rectangular bit subsets: a nonempty subset is a
// rectangle iff its bounding integer box has exactly the subset's cardinality.
bool isRectangle(unsigned mask)
{
  if(!mask) return false;
  std::array<unsigned,3> lo{2,2,2}, hi{0,0,0};
  for(unsigned bit = 0; bit < 8; ++bit) if(mask & (1u << bit))
    for(unsigned axis = 0; axis < 3; ++axis)
    {
      const auto position = (bit >> axis) & 1u;
      lo[axis] = std::min(lo[axis], position);
      hi[axis] = std::max(hi[axis], position);
    }
  return (hi[0]-lo[0]+1)*(hi[1]-lo[1]+1)*(hi[2]-lo[2]+1) == population(mask);
}

std::array<unsigned,256> independentCosts()
{
  std::array<unsigned,256> cost{};
  cost.fill(9); cost[0] = 0;
  for(unsigned mask = 1; mask < 256; ++mask)
    for(unsigned sub = mask; sub; sub = (sub - 1) & mask)
      if(isRectangle(sub)) cost[mask] = std::min(cost[mask], 1 + cost[mask ^ sub]);
  return cost;
}

OctreeCellPartition compare(const std::vector<VoxelPoint>& points, double size)
{
  const auto octree = buildOctreeOccupancy(points, size);
  validateOctreeOccupancy(octree, points);
  const auto base = pruneOctree(octree);
  const auto partition = partitionOctreeCells(octree, base, true);
  validateOctreeCellPartition(octree, partition);
  check(partition.cells.size() <= base.cell_nodes.size(), "Rectangles increased cell count");
  std::vector<HS::Vec3> vertices;
  auto clusters = octreeAabbClusters(partition, vertices);
  ClusterValidation validation;
  validateClusters(clusters, vertices, 0, validation);
  check(validation.failures.empty(), "Rectangular support planes are invalid");
  auto tree = buildClusterTree(clusters, 0, "octree_aabb");
  tree.simplified = HS::simplifyToFixedPoint(tree.raw);
  const auto probe = [&](const HS::Vec3& p) {
    check(HS::evaluateInside(tree.simplified, tree.leaves, p, 0) == octree.contains(p),
          "Rectangular tree changed occupied union");
    check(HS::evaluate(tree.raw, tree.leaves, p) == HS::evaluate(tree.simplified, tree.leaves, p),
          "Simplification changed scalar evaluation");
  };
  for(const auto& p : points) probe({p.x(),p.y(),p.z()});
  for(const auto& p : vertices) probe(p);
  const auto& box = octree.nodes[0].box;
  for(unsigned x = 0; x < 7; ++x) for(unsigned y = 0; y < 7; ++y) for(unsigned z = 0; z < 7; ++z)
    probe({box.lower[0] + (static_cast<double>(x)-1)*(box.upper[0]-box.lower[0])/4,
           box.lower[1] + (static_cast<double>(y)-1)*(box.upper[1]-box.lower[1])/4,
           box.lower[2] + (static_cast<double>(z)-1)*(box.upper[2]-box.lower[2])/4});
  auto reversed = points;
  std::reverse(reversed.begin(), reversed.end());
  const auto other = buildOctreeOccupancy(reversed, size);
  const auto reverse_cut = partitionOctreeCells(other, pruneOctree(other), true);
  std::vector<HS::Vec3> other_vertices;
  auto other_clusters = octreeAabbClusters(reverse_cut, other_vertices);
  const auto other_tree = buildClusterTree(other_clusters, 0, "octree_aabb");
  check(partition.voxel_cell_ids == reverse_cut.voxel_cell_ids &&
        HS::structuralKey(tree.raw) == HS::structuralKey(other_tree.raw), "Unstable rectangle ordering");
  check(tree.leaves.size() == other_tree.leaves.size(), "Unstable rectangular plane count");
  for(std::size_t i = 0; i < tree.leaves.size(); ++i)
  {
    const auto& a = tree.leaves[i]; const auto& b = other_tree.leaves[i];
    check(a.offset == b.offset && a.normal.x == b.normal.x && a.normal.y == b.normal.y && a.normal.z == b.normal.z,
          "Input permutation changed rectangular geometry");
  }
  return partition;
}

int main()
{
  try
  {
    const auto optimum = independentCosts();
    for(unsigned mask = 0; mask < 256; ++mask)
    {
      const auto pieces = rectangularSiblingPartition(mask);
      check(pieces.size() == optimum[mask], "Lookup partition is not minimal");
      unsigned covered = 0;
      for(const auto piece : pieces)
      {
        check(isRectangle(piece) && !(covered & piece) && (piece & mask) == piece,
              "Invalid or overlapping rectangular pattern");
        covered |= piece;
      }
      check(covered == mask, "Sibling pattern dropped or filled a voxel");
      if(!mask) continue;
      std::vector<VoxelPoint> points;
      for(unsigned bit = 0; bit < 8; ++bit) if(mask & (1u << bit))
        points.emplace_back(-1.5 + (bit & 1u), -1.5 + ((bit >> 1) & 1u), -1.5 + ((bit >> 2) & 1u));
      check(compare(points, 1).cells.size() == optimum[mask], "Pipeline disagrees with optimal pattern cover");
    }
    rejects([] { rectangularSiblingPartition(256); });
    std::vector<VoxelPoint> full, plane, shell;
    for(int x = 0; x < 4; ++x) for(int y = 0; y < 4; ++y) for(int z = 0; z < 4; ++z)
    {
      VoxelPoint p(x-1.5,y-1.5,z-1.5); full.push_back(p);
      if(z == 0) plane.push_back(p);
      if(x == 0 || y == 0 || z == 0 || x == 3 || y == 3 || z == 3) shell.push_back(p);
    }
    check(compare(full,1).cells.size() == 1, "Lost larger exact cubic recovery");
    check(compare(plane,1).cells.size() == 4, "Local slabs should not merge across sibling parents");
    compare(shell,1);
    std::mt19937_64 rng(29);
    for(int trial = 0; trial < 16; ++trial)
    {
      std::vector<VoxelPoint> subset;
      for(const auto& p : full) if(rng() % 10 < 7) subset.push_back(p);
      compare(subset,1);
    }
    const std::vector<VoxelPoint> pair{{-0.045,-0.015,0.015},{-0.015,-0.015,0.015}};
    auto duplicate = pair; duplicate.insert(duplicate.end(),pair.begin(),pair.end());
    check(compare(duplicate,0.03).cells.size() == 1, "Decimal/negative-coordinate pair failed to merge");
    const auto octree = buildOctreeOccupancy(pair,0.03);
    const auto partition = partitionOctreeCells(octree,pruneOctree(octree),true);
    auto bad = partition; bad.voxel_cell_ids[0] = 42;
    rejects([&] { validateOctreeCellPartition(octree,bad); });
    bad = partition; bad.base_cell_ids[0] = 42;
    rejects([&] { validateOctreeCellPartition(octree,bad); });
    bad = partition; bad.cells[0].occupied_count = 8;
    rejects([&] { validateOctreeCellPartition(octree,bad); });
    bad = partition; bad.cells[0].sibling_mask = 7;
    rejects([&] { validateOctreeCellPartition(octree,bad); });
    bad = partition; bad.cells[0].box.upper[1] += 0.03;
    rejects([&] { validateOctreeCellPartition(octree,bad); });
    bad = partition; bad.rectangular_groups = 0;
    rejects([&] { validateOctreeCellPartition(octree,bad); });
    bad = partition; bad.rectangular_merge_count = 0;
    rejects([&] { validateOctreeCellPartition(octree,bad); });
    bad = partition; bad.cells[0].lower[0] = std::numeric_limits<std::int64_t>::max();
    rejects([&] { validateOctreeCellPartition(octree,bad); });
    rejects([&] { partitionOctreeCells(octree,pruneOctree(octree,false),true); });
    const auto wide = buildOctreeOccupancy({VoxelPoint(-1e9,0,0),VoxelPoint(1e9,0,0)},1,32);
    const auto sparse = partitionOctreeCells(wide,pruneOctree(wide),true);
    validateOctreeCellPartition(wide,sparse);
    check(sparse.cells.size() == 2 && sparse.rectangular_merge_count == 0, "Sparse root erroneously grouped");
    std::cout << "All 256 masks, rectangle coverage, source mapping and shared tree checks passed\n";
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
