#include <rokae_demo/octree_convex_approximation.hpp>

#include <iostream>
#include <random>
#include <stdexcept>

using namespace rokae_demo;
void check(bool condition, const char* message)
{
  if(!condition) throw std::runtime_error(message);
}
template<class F> void rejects(F f)
{
  bool threw = false;
  try { f(); } catch(const std::exception&) { threw = true; }
  check(threw,"Invalid certificate was accepted");
}

OctreeConvexResult run(const std::vector<VoxelPoint>& points, double h, const OctreeConvexOptions& options)
{
  const auto octree = buildOctreeOccupancy(points,h);
  const auto baseline = partitionOctreeCells(octree,pruneOctree(octree),true);
  auto result = approximateOctreeConvex(octree,baseline,options);
  validateOctreeConvex(octree,baseline,options,result);
  check(result.clusters.size() <= baseline.cells.size(),"Cell count increased");
  auto tree = buildClusterTree(result.clusters,0,"hull_test",true);
  tree.simplified = HS::simplifyToFixedPoint(tree.raw);
  ClusterValidation validation;
  validateClusters(result.clusters,result.source_vertices,0,validation);
  check(validation.failures.empty(),"Shared support validation failed");
  const auto test = [&](const HS::Vec3& p, bool source) {
    const auto a = HS::evaluate(tree.raw,tree.leaves,p);
    const auto b = HS::evaluate(tree.simplified,tree.leaves,p);
    check(a == b,"Tree simplification changed scalar");
    check(!source || b <= 0,"Source point missed");
    check((b <= 0) == insideClusterUnion(result.clusters,p,0),"Plane/tree decision mismatch");
    if(options.max_volume_inflation == 0 || options.max_fill_distance == 0)
      check((b <= 0) == octree.contains(p),"Zero budget changed occupancy");
  };
  for(const auto& p : points) test({p.x(),p.y(),p.z()},true);
  for(const auto& p : result.source_vertices) test(p,true);
  std::mt19937_64 random(714);
  std::uniform_real_distribution<double> t(0,1);
  for(const auto id : octree.leaf_nodes)
  {
    const auto& box = octree.nodes[id].box;
    for(int i = 0; i < 30; ++i)
      test({box.lower[0]+t(random)*(box.upper[0]-box.lower[0]),
            box.lower[1]+t(random)*(box.upper[1]-box.lower[1]),
            box.lower[2]+t(random)*(box.upper[2]-box.lower[2])},true);
  }
  const auto& box = octree.nodes[0].box;
  for(int i = 0; i < 200; ++i)
  {
    HS::Vec3 p{box.lower[0]+t(random)*(box.upper[0]-box.lower[0]),
               box.lower[1]+t(random)*(box.upper[1]-box.lower[1]),
               box.lower[2]+t(random)*(box.upper[2]-box.lower[2])};
    test(p,octree.contains(p));
  }
  return result;
}

int main()
{
  try
  {
    const std::vector<VoxelPoint> l{{.5,.5,.5},{1.5,.5,.5},{.5,1.5,.5}};
    OctreeConvexOptions options;
    options.max_volume_inflation = .2; options.max_fill_distance = .8;
    const auto merged = run(l,1,options);
    check(merged.clusters.size() == 1 && merged.cells[0].non_aabb,"L must become a non-AABB convex prism");
    check(merged.clusters[0].reduced_planes.size() == 7,"L prism must have seven support planes");
    check(std::abs(merged.output_volume-3.5) < 1e-10 && merged.volume_inflation <= .2,"L prism volume budget incorrect");
    check(merged.maximum_fill_distance_bound <= .8,"L fill bound incorrect");
    auto reversed = l; std::reverse(reversed.begin(),reversed.end());
    const auto repeat = run(reversed,1,options);
    check(merged.voxel_cluster_ids == repeat.voxel_cluster_ids && merged.output_volume == repeat.output_volume,
          "Input ordering changed approximation");

    auto budget = options; budget.max_volume_inflation = .1;
    check(run(l,1,budget).rejected_volume > 0,"Volume constraint did not reject L fill");
    budget = options; budget.max_fill_distance = .1;
    check(run(l,1,budget).rejected_distance > 0,"Distance constraint did not reject L fill");
    budget = options; budget.max_candidate_voxels = 2;
    check(run(l,1,budget).rejected_work_limit > 0,"Voxel work cap not enforced");
    budget = options; budget.max_certificate_boxes = 1;
    check(run(l,1,budget).rejected_work_limit > 0,"Certificate work cap not enforced");
    budget = options; budget.max_volume_inflation = 0;
    check(run(l,1,budget).clusters.size() == 2,"Zero volume budget filled a concavity");
    budget = options; budget.max_fill_distance = 0;
    check(run(l,1,budget).clusters.size() == 2,"Zero distance budget filled a concavity");
    budget = {}; // Exact longer rectangle, beyond the old 2x2x2 local pruning.
    const auto rectangle = run({{.5,.5,.5},{1.5,.5,.5},{2.5,.5,.5}},1,budget);
    check(rectangle.clusters.size() == 1 && rectangle.volume_inflation == 0 && !rectangle.cells[0].non_aabb,
          "Exact hierarchical rectangle recovery failed");
    run({{0,0,0}},.03,options);
    std::vector<VoxelPoint> decimal;
    for(const auto& p : l) decimal.emplace_back((p.x()-5)*.03,(p.y()-5)*.03,(p.z()-5)*.03);
    budget = options; budget.max_fill_distance *= .03;
    check(run(decimal,.03,budget).cells[0].non_aabb,"Negative decimal coordinates failed");
    std::vector<VoxelPoint> u;
    for(int x = 0; x < 4; ++x) for(int y = 0; y < 4; ++y)
      if(x == 0 || x == 3 || y == 0) u.emplace_back(x+.5,y+.5,.5);
    budget = options; budget.max_volume_inflation = 2; budget.max_fill_distance = .1;
    const auto protected_u = run(u,1,budget);
    check(!insideClusterUnion(protected_u.clusters,{2,3,.5},0),"Distance budget failed to protect U opening");
    std::mt19937 random(741);
    for(int trial = 0; trial < 12; ++trial)
    {
      std::vector<VoxelPoint> points;
      for(unsigned i = 0; i < 8; ++i)
        if(random()%2 || points.empty()) points.emplace_back((i&1)+.5,((i>>1)&1)+.5,((i>>2)&1)+.5);
      run(points,1,options);
    }

    const auto octree = buildOctreeOccupancy(l,1);
    const auto base = partitionOctreeCells(octree,pruneOctree(octree),true);
    auto bad = merged; bad.clusters[0].exact_planes[0].offset -= .1;
    rejects([&] { validateOctreeConvex(octree,base,options,bad); });
    bad = merged; bad.voxel_cluster_ids[0] = 7;
    rejects([&] { validateOctreeConvex(octree,base,options,bad); });
    bad = merged; bad.cells[0].volume += 1;
    rejects([&] { validateOctreeConvex(octree,base,options,bad); });
    bad = merged; bad.maximum_fill_distance_bound = 0;
    rejects([&] { validateOctreeConvex(octree,base,options,bad); });
    bad = merged; bad.cells[0].fill_distance_bound = 0;
    rejects([&] { validateOctreeConvex(octree,base,options,bad); });
    bad = merged; bad.cells[0].non_aabb = false;
    rejects([&] { validateOctreeConvex(octree,base,options,bad); });
    bad = merged; bad.source_vertices[0].x += .1;
    rejects([&] { validateOctreeConvex(octree,base,options,bad); });
    budget = options; budget.max_fill_distance = std::numeric_limits<double>::quiet_NaN();
    rejects([&] { approximateOctreeConvex(octree,base,budget); });

    Plane a,b; a.normal = {1,0,0}; a.offset = 1;
    b = a; b.normal = {1,1e-8,0};
    std::vector<HS::HalfspaceLeaf> leaves;
    findOrAddLeaf(a,0,0,leaves,"strict",true);
    findOrAddLeaf(b,1,0,leaves,"strict",true);
    check(leaves.size() == 2,"Strict leaf reuse aliased distinct normal coefficients");
    findOrAddLeaf(a,2,0,leaves,"strict",true);
    check(leaves.size() == 2,"Identical support plane was not reused");
    std::cout << "Bounded non-AABB hull, source coverage, budgets, and strict tree tests passed\n";
    return 0;
  }
  catch(const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
