#include <rokae_demo/octree_occupancy.hpp>

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
  bool rejected = false;
  try { function(); } catch(const std::exception&) { rejected = true; }
  check(rejected, "Expected explicit rejection");
}

int main()
{
  try
  {
    const std::vector<VoxelPoint> single{{-0.005, 0.002, 0.005}, {-0.005, 0.002, 0.005}};
    auto one = buildOctreeOccupancy(single, 0.01, 0);
    validateOctreeOccupancy(one, single);
    check(one.nodes.size() == 1 && one.leaf_nodes.size() == 1, "Single voxel root");
    check(one.nodes[0].lower == VoxelIndex{-1, 0, 0}, "Negative coordinate floor");
    check(one.nodes[0].point_indices.size() == 2, "Duplicate raw points retained");
    std::vector<HS::Vec3> vertices;
    auto clusters = octreeAabbClusters(one, vertices);
    ClusterValidation validation;
    validateClusters(clusters, vertices, 0, validation);
    check(validation.failures.empty(), "AABB support checks");
    auto tree = buildClusterTree(clusters, 0, "octree_aabb");
    tree.simplified = HS::simplifyToFixedPoint(tree.raw);
    check(tree.leaves.size() == 6 && HS::serializeTree(tree.simplified).nodes.size() == 7,
          "One box has six planes and one MAX");
    // Rebuilding the same clusters must not append stale leaf references.
    tree = buildClusterTree(clusters, 0, "octree_aabb");
    check(clusters[0].leaf_ids.size() == 6, "Repeated tree construction");
    for(const auto& p : vertices)
      check(one.contains(p) && HS::evaluateInside(tree.raw, tree.leaves, p, 0), "Closed box corners");
    check(!one.contains({std::nextafter(0.0, 1.0), 0.005, 0.005}), "No epsilon expansion");

    std::vector<VoxelPoint> block;
    for(int x = -1; x < 1; ++x)
      for(int y = -1; y < 1; ++y)
        for(int z = -1; z < 1; ++z)
          block.emplace_back(x + 0.5, y + 0.5, z + 0.5);
    auto full = buildOctreeOccupancy(block, 1, 1);
    validateOctreeOccupancy(full, block);
    check(full.nodes.size() == 9 && full.leaf_nodes.size() == 8, "Phase 1 must not prune parents");
    rejects([&] { buildOctreeOccupancy(block, 1, 0); });
    auto broken = full;
    broken.nodes[0].children[0] = -1;
    rejects([&] { validateOctreeOccupancy(broken, block); });
    broken = full;
    broken.nodes[broken.leaf_nodes[0]].point_indices[0] = 7;
    rejects([&] { validateOctreeOccupancy(broken, block); });

    const double s = 0.03;
    std::vector<VoxelPoint> cloud;
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> random(-0.2, 0.2);
    for(int i = 0; i < 100; ++i) cloud.emplace_back(random(rng), random(rng), random(rng));
    // Points exactly on rounded decimal planes and immediately on either side.
    for(int i = -5; i <= 5; ++i)
    {
      const double face = i * s;
      cloud.emplace_back(face, face, face);
      cloud.emplace_back(std::nextafter(face, -1.0), 0, 0);
      cloud.emplace_back(std::nextafter(face, 1.0), 0, 0);
    }
    const auto octree = buildOctreeOccupancy(cloud, s);
    const auto voxel = conservativeVoxelize(cloud, s);
    validateOctreeOccupancy(octree, cloud);
    check(octree.occupied_voxels == voxel.occupied_voxels, "Shared voxel indexing semantics");
    auto reversed = cloud;
    std::reverse(reversed.begin(), reversed.end());
    const auto other = buildOctreeOccupancy(reversed, s);
    check(octree.occupied_voxels == other.occupied_voxels && octree.leaf_nodes == other.leaf_nodes,
          "Permutation invariant geometry and node IDs");
    for(const auto& p : voxel.points)
      check(octree.contains({p.x(), p.y(), p.z()}), "All reference corners contained");
    for(int i = 0; i < 1000; ++i)
    {
      const VoxelPoint p(random(rng), random(rng), random(rng));
      check(octree.contains({p.x(), p.y(), p.z()}) == voxel.contains(p), "Independent occupancy comparison");
    }
    const auto high = buildOctreeOccupancy({VoxelPoint(std::ldexp(1.0, 50) - 0.5, 0, 0)}, 1);
    check(high.contains({std::ldexp(1.0, 50), 0, 0}), "Upper grid boundary containment");
    const std::vector<VoxelPoint> distant{{-1e9, 0, 0}, {1e9, 0, 0}};
    const auto sparse = buildOctreeOccupancy(distant, 1, 32);
    validateOctreeOccupancy(sparse, distant);
    check(sparse.leaf_nodes.size() == 2 && sparse.nodes.size() <= 65,
          "Large empty spans must not allocate dense grids");
    rejects([&] { buildOctreeOccupancy({}, 1); });
    rejects([&] { buildOctreeOccupancy(single, 0); });
    rejects([&] { buildOctreeOccupancy(single, -1); });
    rejects([&] { buildOctreeOccupancy(single, std::numeric_limits<double>::quiet_NaN()); });
    rejects([&] { buildOctreeOccupancy({VoxelPoint(INFINITY, 0, 0)}, 1); });
    rejects([&] { buildOctreeOccupancy(single, 1, 53); });
    rejects([&] { buildOctreeOccupancy({VoxelPoint(1e100, 0, 0)}, 0.01); });
    rejects([&] { buildOctreeOccupancy({VoxelPoint(0, 0, 0)}, 1e200); });
    rejects([&] { buildOctreeOccupancy({VoxelPoint(0, 0, 0)}, 1e-200); });
    std::cout << "Octree occupancy and shared halfspace backend checks passed\n";
    return 0;
  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
