#include <rokae_demo/boundary_local_tree.hpp>
#include <iostream>
#include <random>
#include <stdexcept>

using namespace rokae_demo;
void check(bool b,const char* why) { if(!b) throw std::runtime_error(why); }
std::size_t zeros(const BoundarySignDiagnostics& d)
{ return d.free_zeros[0]+d.free_zeros[1]+d.free_zeros[2]+d.free_zeros[3]; }
void sameDiagnostics(const BoundarySignDiagnostics& a,const BoundarySignDiagnostics& b)
{
  check(a.applied==b.applied && a.closure_set_certified==b.closure_set_certified &&
        a.strict_sign_certified==b.strict_sign_certified && a.reason==b.reason && a.strata_count==b.strata_count &&
        a.unsafe_free_count==b.unsafe_free_count && a.missed_free_count==b.missed_free_count &&
        a.boundary_sign_errors==b.boundary_sign_errors && a.closure_set_errors==b.closure_set_errors &&
        a.raw_simplified_sign_errors==b.raw_simplified_sign_errors && a.free_zeros==b.free_zeros &&
        a.occupied_zeros==b.occupied_zeros && a.witnesses.size()==b.witnesses.size(),
        "Packed/scalar sign diagnostics differ");
  for(std::size_t i=0;i<a.witnesses.size();++i)
    check(a.witnesses[i].strata==b.witnesses[i].strata && a.witnesses[i].dimension==b.witnesses[i].dimension &&
          a.witnesses[i].expected==b.witnesses[i].expected && a.witnesses[i].actual==b.witnesses[i].actual,
          "Packed/scalar sign witnesses differ");
}

void testIsolatedLowerDimensionalZeros()
{
  // An independently assembled strict complement of two separated cubes.
  // Inject a zero line/point without creating any artificial zero faces: the
  // all-strata diagnostic must reject it even if face-only checks would pass.
  const auto o=buildOctreeOccupancy({{.5,.5,.5},{2.5,2.5,2.5}},1);
  const auto p=partitionOctreeCells(o,pruneOctree(o),true);
  const auto r=buildBoundaryTree(o,p,extractOctreeBoundary(o,p));
  validateBoundaryTree(o,p,extractOctreeBoundary(o,p),r);
  const auto leaf=[&](unsigned axis,std::int64_t coordinate,int sign) {
    for(std::size_t i=0;i<r.leaf_sources.size();++i)
    {
      const auto& link=r.leaf_sources[i]; const auto& s=r.supports[link.first];
      if(s.axis==axis && s.coordinate==coordinate && link.second==sign) return HS::makeLeaf(static_cast<int>(i));
    }
    throw std::runtime_error("Missing isolated-zero test plane: axis="+std::to_string(axis)+" coordinate="+
                             std::to_string(coordinate)+" sign="+std::to_string(sign));
  };
  std::vector<HS::ExpressionPtr> cubes;
  for(const auto low : {0,2})
  {
    std::vector<HS::ExpressionPtr> planes;
    for(unsigned axis=0;axis<3;++axis)
    { planes.push_back(leaf(axis,low,1)); planes.push_back(leaf(axis,low+1,-1)); }
    cubes.push_back(HS::makeNode(HS::TreeOperator::MIN,std::move(planes)));
  }
  auto ideal=r.tree; ideal.raw=HS::makeNode(HS::TreeOperator::MAX,std::move(cubes)); ideal.simplified=ideal.raw;
  check(diagnoseBoundarySigns(o,r,ideal).strict_sign_certified,"Independent cube complement is not strict");
  for(unsigned axes=2;axes<=3;++axes)
  {
    std::vector<HS::ExpressionPtr> absolute_terms;
    for(unsigned axis=0;axis<axes;++axis)
    {
      const auto coordinate=axis==1 ? 2 : 0;
      absolute_terms.push_back(leaf(axis,coordinate,1)); absolute_terms.push_back(leaf(axis,coordinate,-1));
    }
    auto bad=ideal;
    bad.raw=HS::makeNode(HS::TreeOperator::MAX,{ideal.raw,HS::makeNode(HS::TreeOperator::MIN,std::move(absolute_terms))});
    bad.simplified=bad.raw; const auto d=diagnoseBoundarySigns(o,r,bad);
    check(d.closure_set_certified && !d.strict_sign_certified && !d.unsafe_free_count,"Artificial zero test changed closure geometry");
    check(d.free_zeros[2]==0 && d.free_zeros[0]>0,"Diagnostic missed isolated lower-dimensional zeros");
    check((d.free_zeros[1]>0)==(axes==2),"Point/line zero diagnostic did not distinguish dimensions");
  }
}

BoundaryLocalResult test(const std::vector<VoxelPoint>& points,unsigned axis,bool require_strict,double h=1)
{
  const auto o = buildOctreeOccupancy(points,h);
  const auto p = partitionOctreeCells(o,pruneOctree(o),true); const auto b = extractOctreeBoundary(o,p);
  BoundaryTreeOptions opts; opts.sweep_axis = axis;
  const auto flat = buildBoundaryTree(o,p,b,opts); validateBoundaryTree(o,p,b,flat);
  const auto before = diagnoseBoundarySigns(o,flat,flat.tree);
  BoundaryLocalOptions scalar; scalar.packed_sign_propagation=false;
  sameDiagnostics(before,diagnoseBoundarySigns(o,flat,flat.tree,scalar));
  check(before.applied && before.closure_set_certified,"Flat sign certificate failed");
  auto local = buildBoundaryLocalTree(flat);
  local.diagnostics = diagnoseBoundarySigns(o,flat,local.tree);
  const auto& d = local.diagnostics;
  sameDiagnostics(d,diagnoseBoundarySigns(o,flat,local.tree,scalar));
  check(d.applied && d.closure_set_certified && !d.unsafe_free_count && !d.boundary_sign_errors,"Local rewrite changed geometry/boundary");
  check(zeros(d) <= zeros(before),"Rewrite introduced new artificial zeros");
  check(!require_strict || d.strict_sign_certified,"Required strict local fixture still has artificial zeros");
  for(const auto& point : points)
  {
    const HS::Vec3 q{point.x(),point.y(),point.z()};
    check(!boundaryLocalStrictFree(local,o,q),"Raw point classified free");
  }
  // Independent world-coordinate checks at all finite strata in the small fixtures.
  std::array<std::vector<double>,3> positions;
  for(unsigned a = 0; a < 3; ++a)
  {
    const auto& c = flat.coordinates[a]; positions[a].push_back(static_cast<double>(c.front()-1)*h);
    for(std::size_t i = 0; i < c.size(); ++i)
    {
      positions[a].push_back(static_cast<double>(c[i])*h);
      if(i+1 < c.size()) positions[a].push_back(static_cast<double>(c[i])*h+(static_cast<double>(c[i+1])*h-static_cast<double>(c[i])*h)/2);
    }
    positions[a].push_back(static_cast<double>(c.back()+1)*h);
  }
  for(double x : positions[0]) for(double y : positions[1]) for(double z : positions[2])
  {
    HS::Vec3 q{x,y,z};
    check(boundaryLocalStrictFree(local,o,q) == !o.contains(q),"Local world classifier mismatch");
    check(HS::evaluate(local.tree.raw,local.tree.leaves,q) == HS::evaluate(local.tree.simplified,local.tree.leaves,q),"Simplifier changed local scalar");
  }
  auto bad = local.tree; bad.leaves[0].offset += .25; bool rejected = false;
  try { diagnoseBoundarySigns(o,flat,bad); } catch(const std::exception&) { rejected = true; }
  check(rejected,"Changed plane dictionary accepted");
  BoundaryLocalOptions cap; cap.max_strata = 1;
  check(!diagnoseBoundarySigns(o,flat,local.tree,cap).applied,"Sign work cap ignored");
  cap = {}; cap.max_sign_operations = 1;
  check(!diagnoseBoundarySigns(o,flat,local.tree,cap).strict_sign_certified,"Skipped certificate claimed success");
  for(unsigned mode=0;mode<2;++mode)
  {
    cap = {}; if(mode) cap.max_attempts=1; else cap.max_adjacency_checks=1;
    const auto limited=buildBoundaryLocalTree(flat,cap);
    const auto certificate=diagnoseBoundarySigns(o,flat,limited.tree);
    check(certificate.closure_set_certified && !certificate.unsafe_free_count,"Work cap made local geometry unsafe");
  }
  auto reverse = points; std::reverse(reverse.begin(),reverse.end());
  const auto other = buildOctreeOccupancy(reverse,h); const auto cut = partitionOctreeCells(other,pruneOctree(other),true);
  const auto again = buildBoundaryLocalTree(buildBoundaryTree(other,cut,extractOctreeBoundary(other,cut),opts));
  check(HS::structuralKey(local.tree.raw) == HS::structuralKey(again.tree.raw),"Local construction is nondeterministic");
  return local;
}

int main()
{
  try
  {
    testIsolatedLowerDimensionalZeros();
    std::vector<VoxelPoint> l_cavity,u_cavity,stair_cavity,hollow;
    for(int x=-2;x<2;++x) for(int y=-1;y<3;++y) for(int z=-1;z<2;++z)
    {
      const VoxelPoint p(x+.5,y+.5,z+.5);
      if(!(z==0 && ((x==-1 && y==0) || (x==0 && (y==0 || y==1))))) l_cavity.push_back(p);
    }
    for(int x=-1;x<4;++x) for(int y=-1;y<4;++y) for(int z=-1;z<2;++z)
    {
      const VoxelPoint p(x+.5,y+.5,z+.5);
      const bool inner=z==0 && x>=0 && x<3 && y>=0 && y<3;
      if(!(inner && (x!=1 || y==0))) u_cavity.push_back(p);
      if(!(inner && y<=x)) stair_cavity.push_back(p);
      if(!inner) hollow.push_back(p);
    }
    for(unsigned axis=0;axis<3;++axis)
    {
      test({{.5,.5,.5}},axis,true);
      test({{.5,.5,.5},{1.5,.5,.5}},axis,true); // exact prism/bbox case
      const auto l=test(l_cavity,axis,true);
      check(l.logical_merges && !l.diagnostics.free_zeros[2],"L template did not run");
      check(HS::evaluate(l.tree.raw,l.tree.leaves,{0,.5,.5})<0,"L internal seam not repaired");
      check(HS::evaluate(l.tree.raw,l.tree.leaves,{-.5,1.5,.5})>0,"L notch incorrectly filled");
      check(HS::evaluate(l.tree.raw,l.tree.leaves,{0,1.5,.5})==0,"L true boundary removed");
      test(hollow,axis,true); test(u_cavity,axis,false); test(stair_cavity,axis,false);
      // Minimal documented limitation: a finite gap touches four exterior
      // halfspaces. Root-gate matches disappear after the first logical unions.
      const auto gap=test({{.5,.5,.5},{2.5,.5,.5}},axis,false);
      check(gap.adjacency.size()==4 && gap.logical_merges>0,"Exterior interfaces omitted");
      check(std::all_of(gap.adjacency.begin(),gap.adjacency.end(),[](const auto& f) { return f.exterior; }),"Wrong gap graph");
      check(!gap.diagnostics.strict_sign_certified && zeros(gap.diagnostics)>0,"Expected bounded-template counterexample disappeared; update test/proof");
      check(HS::evaluate(gap.tree.raw,gap.tree.leaves,{1.5,.5,0})==0,"Gap witness changed; update counterexample");
      test({{-.015,-.015,-.015},{.015,.015,.015}},axis,false,.03);
    }
    std::mt19937 generator(947);
    for(int trial=0;trial<12;++trial)
    {
      std::vector<VoxelPoint> cloud;
      for(int x=0;x<3;++x) for(int y=0;y<3;++y) for(int z=0;z<3;++z)
        if(generator()%3==0) cloud.emplace_back(x+.5,y+.5,z+.5);
      if(cloud.empty()) cloud.emplace_back(.5,.5,.5);
      test(cloud,static_cast<unsigned>(trial%3),false);
    }
    std::cout << "Local sign strata, L rewrite, exterior graph, safe incomplete results and caps passed\n";
    return 0;
  }
  catch(const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
