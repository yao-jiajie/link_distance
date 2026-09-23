#include <rokae_demo/octree_boundary_tree.hpp>
#include <iostream>
#include <random>
#include <stdexcept>
#include <tuple>

using namespace rokae_demo;
void check(bool v, const char* msg) { if(!v) throw std::runtime_error(msg); }
template<class F> void rejects(F f)
{ bool rejected = false; try { f(); } catch(const std::exception&) { rejected = true; } check(rejected,"Expected certificate rejection"); }

void testAxis(const std::vector<VoxelPoint>& points, double h, BoundaryTreeOptions options)
{
  const auto octree = buildOctreeOccupancy(points,h);
  const auto base = partitionOctreeCells(octree,pruneOctree(octree),true);
  const auto boundary = extractOctreeBoundary(octree,base);
  const auto r = buildBoundaryTree(octree,base,boundary,options);
  check(r.sweep_axis == options.sweep_axis,"Sweep axis was ignored");
  validateBoundaryTree(octree,base,boundary,r);
  // Exhaust all Cartesian strata for these small geometries: open cells,
  // faces, edges and vertices. Extend outside the bbox too.
  std::array<std::vector<double>,3> coordinates;
  for(unsigned k = 0; k < 3; ++k)
  {
    const auto& c = r.coordinates[k];
    coordinates[k].push_back(static_cast<double>(c.front())*h-h);
    for(std::size_t i = 0; i < c.size(); ++i)
    {
      coordinates[k].push_back(static_cast<double>(c[i])*h);
      if(i+1 < c.size()) coordinates[k].push_back(static_cast<double>(c[i])*h +
         (static_cast<double>(c[i+1])*h-static_cast<double>(c[i])*h)/2);
    }
    coordinates[k].push_back(static_cast<double>(c.back())*h+h);
  }
  for(const auto x : coordinates[0]) for(const auto y : coordinates[1]) for(const auto z : coordinates[2])
  {
    const HS::Vec3 p{x,y,z};
    check(boundaryStrictFree(r,octree,p) == !octree.contains(p),"Strict free classification/strata mismatch");
    check(HS::evaluate(r.tree.raw,r.tree.leaves,p) == HS::evaluate(r.tree.simplified,r.tree.leaves,p),"Simplification changed scalar");
  }
  for(const auto& p : points) check(!boundaryStrictFree(r,octree,{p.x(),p.y(),p.z()}),"Raw point classified free");
  auto reverse = points; std::reverse(reverse.begin(),reverse.end());
  const auto other = buildOctreeOccupancy(reverse,h);
  const auto other_base = partitionOctreeCells(other,pruneOctree(other),true);
  const auto again = buildBoundaryTree(other,other_base,extractOctreeBoundary(other,other_base),options);
  check(HS::structuralKey(again.tree.raw) == HS::structuralKey(r.tree.raw),"Input ordering changed tree");
  auto bad = r; bad.tree.leaves[0].offset += .25;
  rejects([&] { validateBoundaryTree(octree,base,boundary,bad); });
  bad = r; bad.supports[0].outward_signs = 0;
  rejects([&] { validateBoundaryTree(octree,base,boundary,bad); });
  if(!r.fallback && !r.prisms.empty())
  {
    bad = r; bad.prisms.push_back(bad.prisms[0]);
    rejects([&] { validateBoundaryTree(octree,base,boundary,bad); });
    bad = r; bad.prisms.pop_back();
    rejects([&] { validateBoundaryTree(octree,base,boundary,bad); });
  }
}

void test(const std::vector<VoxelPoint>& points, double h = 1, BoundaryTreeOptions options = {})
{
  for(unsigned axis = 0; axis < 3; ++axis)
  { options.sweep_axis = axis; testAxis(points,h,options); }
}

BoundaryTreeBatch testBatch(const std::vector<VoxelPoint>& points, BoundaryTreeOptions options = {})
{
  const auto o = buildOctreeOccupancy(points,1);
  const auto b = partitionOctreeCells(o,pruneOctree(o),true);
  const auto boundary = extractOctreeBoundary(o,b);
  const auto batch = buildBoundaryTreeBatch(o,b,boundary,options,true);
  check(batch.candidates.size() == 3 && batch.candidate_build_ms.size() == 3,"Best-axis must build three candidates");
  const auto key = [](const auto& r) {
    return std::make_tuple(r.fallback,HS::serializeTree(r.tree.simplified).nodes.size(),r.tree.leaf_references,r.prisms.size(),r.sweep_axis);
  };
  const auto& winner = batch.candidates.at(batch.selected);
  double total = 0;
  for(unsigned axis = 0; axis < 3; ++axis)
  {
    const auto& r = batch.candidates[axis];
    check(r.sweep_axis == axis && key(winner) <= key(r),"Wrong best-axis selection");
    validateBoundaryTree(o,b,boundary,r);
    options.sweep_axis = axis;
    const auto alone = buildBoundaryTreeBatch(o,b,boundary,options,false);
    check(alone.candidates.size() == 1 && alone.selected == 0,"Single-axis batch has extra candidates");
    check(HS::structuralKey(r.tree.simplified) == HS::structuralKey(alone.candidates[0].tree.simplified),"Batch changed candidate geometry");
    total += batch.candidate_build_ms[axis];
  }
  check(batch.all_candidates_ms >= total && batch.selection_ms >= 0,"Best-axis build cost omitted candidates");
  return batch;
}

int main()
{
  try
  {
    test({{.5,.5,.5}}); test({{0,0,0}},.03);
    test({{.5,.5,.5},{1.5,.5,.5}});
    test({{.5,.5,.5},{1.5,.5,.5},{.5,1.5,.5}});
    test({{.5,.5,.5},{1.5,1.5,.5}}); // edge-only contact
    test({{.5,.5,.5},{1.5,1.5,1.5}}); // vertex-only contact
    test({{-.015,-.015,-.015},{.015,.015,.015}},.03);
    std::vector<VoxelPoint> solid,shell,u,stair,wall,corridor;
    for(int x = 0; x < 3; ++x) for(int y = 0; y < 3; ++y) for(int z = 0; z < 3; ++z)
    {
      const VoxelPoint p(x+.5,y+.5,z+.5); solid.push_back(p);
      if(x != 1 || y != 1 || z != 1) shell.push_back(p);
      if(x != 1 || y == 0) u.push_back(p);
      if(y <= x) stair.push_back(p);
      if(x == 1) wall.push_back(p);
      if(x != 1) corridor.push_back(p);
    }
    for(const auto& cloud : {solid,shell,u,stair,wall,corridor}) test(cloud);
    check(testBatch(solid).selected == 0,"Equal trees must tie-break X then Y then Z");
    const auto all = testBatch(u);
    std::size_t min_events = std::numeric_limits<std::size_t>::max(), max_events = 0;
    for(const auto& c : all.candidates)
    { min_events = std::min(min_events,c.event_updates); max_events = std::max(max_events,c.event_updates); }
    check(min_events < max_events,"Partial-fallback fixture must distinguish axes");
    BoundaryTreeOptions mixed_cap; mixed_cap.max_event_updates = min_events;
    const auto mixed = testBatch(u,mixed_cap);
    check(!mixed.candidates[mixed.selected].fallback,"Fallback selected over successful boundary candidate");
    check(std::any_of(mixed.candidates.begin(),mixed.candidates.end(),[](const auto& c) { return c.fallback; }),"Partial fallback not exercised");
    BoundaryTreeOptions cap; cap.max_compressed_cells = 1; test(shell,1,cap);
    const auto capped = testBatch(shell,cap);
    check(capped.selected == 0 && capped.candidates[0].fallback,"All-fallback batch must use original rect from X");
    cap = {}; cap.max_event_updates = 1; test(shell,1,cap);
    cap = {}; cap.max_free_prisms = 1; test(u,1,cap);
    cap = {}; cap.max_merge_passes = 1; test(u,1,cap);
    std::mt19937 random(731);
    for(int trial = 0; trial < 16; ++trial)
    {
      std::vector<VoxelPoint> cloud;
      for(const auto& p : solid) if(random()%3 == 0) cloud.push_back(p);
      if(cloud.empty()) cloud.push_back(solid[0]);
      test(cloud);
    }
    // Tiny work cap prevents a large compressed Cartesian product.
    std::vector<VoxelPoint> distant;
    for(int i = 0; i < 30; ++i) distant.emplace_back(10*i+.5,13*i+.5,17*i+.5);
    auto o = buildOctreeOccupancy(distant,1); auto b = partitionOctreeCells(o,pruneOctree(o),true);
    auto r = buildBoundaryTree(o,b,extractOctreeBoundary(o,b));
    check(r.fallback && r.fallback_reason == "compressed_cell_limit","Sparse Cartesian work cap failed");
    validateBoundaryTree(o,b,extractOctreeBoundary(o,b),r);
    BoundaryTreeOptions invalid; invalid.sweep_axis = 3;
    rejects([&] { buildBoundaryTree(o,b,extractOctreeBoundary(o,b),invalid); });
    rejects([&] { boundaryStrictFree(r,o,{NAN,0,0}); });
    std::cout << "Boundary XYZ/best-axis, closure certificate, strict strata, deterministic selection and fallback tests passed\n";
    return 0;
  }
  catch(const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
