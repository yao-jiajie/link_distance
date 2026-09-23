#include <rokae_demo/octree_occupancy.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>

using namespace rokae_demo;

void check(bool condition, const char* message)
{
  if(!condition) throw std::runtime_error(message);
}

template<class F> void rejects(F function)
{
  bool rejected = false;
  try { function(); } catch(const std::exception&) { rejected = true; }
  check(rejected, "Expected rejection of invalid pruning cut");
}

std::vector<VoxelPoint> block(int width)
{
  std::vector<VoxelPoint> points;
  for(int x = 0; x < width; ++x)
    for(int y = 0; y < width; ++y)
      for(int z = 0; z < width; ++z)
        points.emplace_back(x - 2 + 0.5, y - 2 + 0.5, z - 2 + 0.5);
  return points;
}

void compare(const std::vector<VoxelPoint>& points, double size)
{
  const auto octree = buildOctreeOccupancy(points, size);
  validateOctreeOccupancy(octree, points);
  const auto pruned = pruneOctree(octree);
  validateOctreePruning(octree, pruned);
  const auto none = pruneOctree(octree, false);
  validateOctreePruning(octree, none);
  check(pruned.cell_nodes.size() <= none.cell_nodes.size(), "Pruning increased cell count");
  std::vector<HS::Vec3> before_vertices, after_vertices;
  auto before_clusters = octreeAabbClusters(octree, none, before_vertices);
  auto after_clusters = octreeAabbClusters(octree, pruned, after_vertices);
  ClusterValidation validation;
  validateClusters(after_clusters, after_vertices, 0, validation);
  check(validation.failures.empty(), "Merged box halfspaces failed support validation");
  auto before = buildClusterTree(before_clusters, 0, "octree_aabb");
  auto after = buildClusterTree(after_clusters, 0, "octree_aabb");
  after.simplified = HS::simplifyToFixedPoint(after.raw);
  const auto probe = [&](const HS::Vec3& p) {
    const bool reference = octree.contains(p);
    check(HS::evaluateInside(before.raw, before.leaves, p, 0) == reference, "Fine tree mismatch");
    check(HS::evaluateInside(after.simplified, after.leaves, p, 0) == reference, "Pruned tree mismatch");
  };
  for(const auto& p : before_vertices) probe(p);  // Includes removed internal faces
  for(const auto& p : after_vertices) probe(p);
  for(const auto& p : points) probe({p.x(), p.y(), p.z()});
  std::mt19937_64 generator(19);
  std::uniform_real_distribution<double> sample(-3, 3);
  for(int i = 0; i < 300; ++i) probe({sample(generator), sample(generator), sample(generator)});
  // Repeated pruning and input permutations have stable cut/plane/tree ordering.
  const auto again = pruneOctree(octree);
  check(pruned.cell_nodes == again.cell_nodes && pruned.voxel_cell_ids == again.voxel_cell_ids,
        "Pruning is not deterministic");
  auto reversed = points;
  std::reverse(reversed.begin(), reversed.end());
  const auto other = buildOctreeOccupancy(reversed, size);
  const auto reverse_cut = pruneOctree(other);
  check(pruned.cell_nodes == reverse_cut.cell_nodes && pruned.voxel_cell_ids == reverse_cut.voxel_cell_ids,
        "Input permutation changed the pruning cut");
}

int main()
{
  try
  {
    const auto points = block(4);
    const auto octree = buildOctreeOccupancy(points, 1);
    const auto pruned = pruneOctree(octree);
    validateOctreePruning(octree, pruned);
    check(pruned.cell_nodes == std::vector<int>{0} && pruned.merge_operations == 9,
          "4x4x4 occupancy must recover root via nine eight-child merges");
    check(std::count(pruned.active_nodes.begin(), pruned.active_nodes.end(), 1) == 1,
          "Descendants of recovered root remain active");
    std::vector<HS::Vec3> vertices;
    auto merged = octreeAabbClusters(octree, pruned, vertices);
    auto fine = octreeAabbClusters(octree, vertices);
    const auto new_tree = buildClusterTree(merged, 0, "octree_aabb");
    const auto old_tree = buildClusterTree(fine, 0, "octree_aabb");
    check(HS::evaluate(old_tree.raw, old_tree.leaves, {0,0,0}) == 0 &&
          HS::evaluate(new_tree.raw, new_tree.leaves, {0,0,0}) == -2,
          "Geometrically exact pruning can change scalar values on internal faces");
    check(HS::serializeTree(new_tree.raw).nodes.size() == 7, "Recovered root tree is not minimal AABB");
    auto bad = pruned;
    bad.cell_nodes.push_back(octree.leaf_nodes.front());
    rejects([&] { validateOctreePruning(octree, bad); });
    bad = pruned; bad.voxel_cell_ids[0] = 1;
    rejects([&] { validateOctreePruning(octree, bad); });
    bad = pruned; bad.active_nodes[1] = 1;
    rejects([&] { validateOctreePruning(octree, bad); });
    bad = pruned; bad.merge_operations = 8;
    rejects([&] { validateOctreePruning(octree, bad); });
    bad = pruneOctree(octree, false); bad.exact = true;
    rejects([&] { validateOctreePruning(octree, bad); });  // not maximally recovered

    compare(points, 1);
    compare(block(2), 1);
    auto missing = points;
    missing.pop_back();
    const auto hole = buildOctreeOccupancy(missing, 1);
    const auto partial = pruneOctree(hole);
    check(partial.cell_nodes.size() == 14 && partial.merge_operations == 7,
          "Missing voxel must prevent ancestor recovery");
    bad = partial; bad.cell_nodes[0] = 0;
    rejects([&] { validateOctreePruning(hole, bad); });
    compare(missing, 1);
    auto shell = points;
    shell.erase(std::remove_if(shell.begin(), shell.end(), [](const VoxelPoint& p) {
      return std::abs(p.x()) < 1 && std::abs(p.y()) < 1 && std::abs(p.z()) < 1;
    }), shell.end());
    compare(shell, 1);  // internal cavity must not be filled
    auto noisy_duplicates = points;
    noisy_duplicates.insert(noisy_duplicates.end(), points.begin(), points.end());
    compare(noisy_duplicates, 1);
    auto decimal = points;
    for(auto& p : decimal) p = VoxelPoint(p.x() * 0.03, p.y() * 0.03, p.z() * 0.03);
    compare(decimal, 0.03);
    compare({VoxelPoint(0, 0, 0)}, 0.03);
    std::mt19937_64 random_subset(71);
    for(int trial = 0; trial < 12; ++trial)
    {
      std::vector<VoxelPoint> subset;
      for(const auto& p : points)
        if(random_subset() % 100 < static_cast<unsigned>(40 + 5 * trial)) subset.push_back(p);
      compare(subset, 1);
    }
    // No width^3 overflow or dense traversal for enormous, mostly empty roots.
    const auto wide = buildOctreeOccupancy({VoxelPoint(-1e9,0,0), VoxelPoint(1e9,0,0)}, 1, 32);
    const auto sparse_cut = pruneOctree(wide);
    validateOctreePruning(wide, sparse_cut);
    check(sparse_cut.merge_operations == 0 && sparse_cut.cell_nodes.size() == 2, "Sparse root incorrectly full");
    std::cout << "Exact octree pruning, coverage, and tree consistency tests passed\n";
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
