#include <rokae_demo/boundary_direct_tree.hpp>
#include <rokae_demo/halfspace_compiled_tree.hpp>
#include <iostream>
#include <random>

using namespace rokae_demo;
void check(bool ok,const char* why) { if(!ok) throw std::runtime_error(why); }
template<class F> void rejects(F f,const char* why)
{ bool rejected=false; try { f(); } catch(const std::exception&) { rejected=true; } check(rejected,why); }
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

void sameSerialized(const HS::SerializedTree& a,const HS::SerializedTree& b)
{
  check(a.root_id==b.root_id && a.nodes.size()==b.nodes.size(),
        "Fast/full direct serialized size differs");
  for(std::size_t i=0;i<a.nodes.size();++i)
    check(a.nodes[i].id==b.nodes[i].id && a.nodes[i].op==b.nodes[i].op &&
          a.nodes[i].leaf_id==b.nodes[i].leaf_id &&
          a.nodes[i].children==b.nodes[i].children,
          "Fast/full direct serialized DAG differs");
}

void test(const std::vector<VoxelPoint>& points,double h=1,bool deep=true)
{
  const auto o=buildSparseVoxelOccupancy(points,h); validateSparseVoxelOccupancy(o,points);
  const auto b=extractVoxelBoundary(o.occupied_voxels); auto r=buildBoundaryDirectTree(o,b);
  const auto d=validateBoundaryDirectTree(o,b,r);
  check(d.applied && d.strict_sign_certified,"Direct tree failed strict all-strata signs");
  const auto fast=buildBoundaryDirectTree(o);
  check(!fast.registry.patch_provenance_complete,
        "Fast direct registry unexpectedly claims patch provenance");
  check(validateBoundaryDirectTree(o,b,fast).strict_sign_certified,
        "Fast direct tree failed strict all-strata signs");
  check(fast.registry.coordinates==r.registry.coordinates &&
        fast.registry.supports.size()==r.registry.supports.size() &&
        fast.leaf_sources==r.leaf_sources,
        "Fast/full direct registry size differs");
  for(std::size_t i=0;i<r.registry.supports.size();++i)
  {
    const auto& actual=fast.registry.supports[i];
    const auto& expected=r.registry.supports[i];
    check(actual.axis==expected.axis && actual.coordinate==expected.coordinate &&
          actual.outward_signs==expected.outward_signs && actual.patches.empty(),
          "Fast/full direct support geometry differs");
  }
  sameSerialized(fast.serialized,r.serialized);
  check(fast.tree.leaves.size()==r.tree.leaves.size(),
        "Fast/full direct leaf count differs");
  for(std::size_t i=0;i<r.tree.leaves.size();++i)
    check(fast.tree.leaves[i].normal.x==r.tree.leaves[i].normal.x &&
          fast.tree.leaves[i].normal.y==r.tree.leaves[i].normal.y &&
          fast.tree.leaves[i].normal.z==r.tree.leaves[i].normal.z &&
          fast.tree.leaves[i].offset==r.tree.leaves[i].offset,
          "Fast/full direct leaf coefficient differs");
  BoundaryDirectOptions scalar_options; scalar_options.signs.packed_sign_propagation=false;
  sameDiagnostics(d,validateBoundaryDirectTree(o,b,r,scalar_options));
  check(r.tree.raw==r.tree.simplified,"Direct unexpectedly invokes simplifier");
  check(HS::treeNodeReferenceCount(r.tree.raw)==r.expanded_nodes,"Expanded reference counter mismatch");
  const auto independently_serialized=HS::serializeTree(r.tree.raw);
  check(r.expanded_plane_refs==r.tree.leaf_references &&
        r.constructed_nodes==r.serialized.nodes.size() &&
        r.serialized.root_id==independently_serialized.root_id &&
        r.serialized.nodes.size()==independently_serialized.nodes.size(),
        "Tree statistics/serialized snapshot mismatch");
  for(std::size_t i=0;i<r.serialized.nodes.size();++i)
  {
    const auto& actual=r.serialized.nodes[i];
    const auto& expected=independently_serialized.nodes[i];
    check(actual.id==expected.id && actual.op==expected.op &&
          actual.leaf_id==expected.leaf_id && actual.children==expected.children,
          "Construction-emitted serialization differs from tree traversal");
  }
  const HS::CompiledTree compiled(r.serialized,r.tree.leaves); auto workspace=compiled.makeWorkspace();
  check(compiled.nodeCount()==r.constructed_nodes,"Compiled DAG node count mismatch");
  BoundaryDirectOptions intern_options; intern_options.intern_subexpressions=true;
  const auto interned=buildBoundaryDirectTree(o,b,intern_options);
  auto serialized_options=intern_options;
  serialized_options.materialize_expression=false;
  const auto serialized_only=buildBoundaryDirectTree(o,b,serialized_options);
  const auto fast_serialized_only=buildBoundaryDirectTree(o,serialized_options);
  sameSerialized(fast_serialized_only.serialized,serialized_only.serialized);
  check(fast_serialized_only.leaf_sources==serialized_only.leaf_sources &&
        validateBoundaryDirectTree(o,b,fast_serialized_only,serialized_options).strict_sign_certified,
        "Fast serialized-only direct result differs or failed validation");
  check(!serialized_only.tree.raw && !serialized_only.tree.simplified,
        "Serialized-only direct build unexpectedly materialized expressions");
  const auto serialized_diagnostics=validateBoundaryDirectTree(o,b,serialized_only,serialized_options);
  check(serialized_diagnostics.applied && serialized_diagnostics.strict_sign_certified,
        "Serialized-only direct tree failed lazy strict validation");
  const auto independently_interned=HS::serializeTree(interned.tree.raw);
  check(interned.serialized.root_id==independently_interned.root_id &&
        interned.serialized.nodes.size()==independently_interned.nodes.size() &&
        serialized_only.serialized.root_id==interned.serialized.root_id &&
        serialized_only.serialized.nodes.size()==interned.serialized.nodes.size() &&
        serialized_only.leaf_sources==interned.leaf_sources,
        "Interned direct snapshot size mismatch");
  for(std::size_t i=0;i<interned.serialized.nodes.size();++i)
  {
    const auto& actual=interned.serialized.nodes[i];
    const auto& expected=independently_interned.nodes[i];
    const auto& serialized=serialized_only.serialized.nodes[i];
    check(actual.id==expected.id && actual.op==expected.op &&
          actual.leaf_id==expected.leaf_id && actual.children==expected.children &&
          serialized.id==actual.id && serialized.op==actual.op &&
          serialized.leaf_id==actual.leaf_id && serialized.children==actual.children,
          "Interned direct snapshot differs from tree traversal");
  }
  const HS::CompiledTree legacy_shared(r.serialized,r.tree.leaves,true);
  const HS::CompiledTree construction_shared(
      interned.serialized,interned.tree.leaves,HS::preinterned_program);
  const HS::CompiledTree serialized_shared(
      serialized_only.serialized,serialized_only.tree.leaves,HS::preinterned_program);
  auto legacy_workspace=legacy_shared.makeWorkspace();
  auto construction_workspace=construction_shared.makeWorkspace();
  auto serialized_workspace=serialized_shared.makeWorkspace();
  check(HS::structuralKey(r.tree.raw)==HS::structuralKey(interned.tree.raw) &&
        legacy_shared.nodeCount()==construction_shared.nodeCount(),
        "Construction-time interning changed the ordered expression DAG");
  for(const auto& split:r.splits)
  {
    const auto& plane=r.registry.supports.at(split.support);
    check(!plane.patches.empty() && split.cut>split.lower[plane.axis] && split.cut<split.upper[plane.axis],"Invalid local split");
    check(r.registry.coordinates[plane.axis][split.cut-1]==plane.coordinate,"Artificial split coordinate");
  }
  const auto query=[&](const HS::Vec3& p) {
    const double q=HS::evaluate(r.tree.raw,r.tree.leaves,p); const int sign=(q>0)-(q<0);
    const double fast=compiled.evaluate(p,workspace);
    check(q==fast && std::signbit(q)==std::signbit(fast),"Compiled boundary scalar/sign bit mismatch");
    const double legacy=legacy_shared.evaluate(p,legacy_workspace);
    const double construction=construction_shared.evaluate(p,construction_workspace);
    const double serialized=serialized_shared.evaluate(p,serialized_workspace);
    check(legacy==construction && std::signbit(legacy)==std::signbit(construction) &&
          serialized==construction && std::signbit(serialized)==std::signbit(construction),
          "Construction-time/legacy interning value mismatch");
    const int expected=boundaryDirectExpectedSign(r,h,p);
    check(sign==expected,"Direct world sign mismatch");
    check((q<0)==!o.contains(p),"Direct world strict classifier mismatch");
  };
  for(const auto& p:points) query({p.x(),p.y(),p.z()});
  std::array<std::vector<double>,3> coordinates;
  for(unsigned a=0;a<3;++a)
  {
    const auto& c=r.registry.coordinates[a]; coordinates[a].push_back(static_cast<double>(c.front())*h-h);
    for(std::size_t i=0;i<c.size();++i)
    {
      const double v=static_cast<double>(c[i])*h; coordinates[a].push_back(v);
      if(i+1<c.size()) coordinates[a].push_back(v+(static_cast<double>(c[i+1])*h-v)/2);
    }
    coordinates[a].push_back(static_cast<double>(c.back())*h+h);
  }
  for(const auto x:coordinates[0]) for(const auto y:coordinates[1]) for(const auto z:coordinates[2]) query({x,y,z});
  if(!deep) return; // all 255 patterns compare DAG/recursive values on every stratum
  for(const auto& k:o.occupied_voxels)
  {
    const auto box=octreeGridBox(k,{1,1,1},h);
    for(unsigned mask=0;mask<8;++mask)
    {
      const double v[3]{mask&1 ? box.upper[0] : box.lower[0],mask&2 ? box.upper[1] : box.lower[1],mask&4 ? box.upper[2] : box.lower[2]};
      for(unsigned axis=0;axis<3;++axis) for(int side:{-1,1})
      {
        auto p=HS::Vec3{v[0],v[1],v[2]}; const auto delta=std::nextafter(v[axis],side<0 ? -INFINITY : INFINITY);
        if(axis==0) p.x=delta; else if(axis==1) p.y=delta; else p.z=delta; query(p);
      }
    }
  }
  auto reversed=points; std::reverse(reversed.begin(),reversed.end());
  const auto oo=buildSparseVoxelOccupancy(reversed,h); const auto rr=buildBoundaryDirectTree(oo,extractVoxelBoundary(oo.occupied_voxels));
  check(HS::structuralKey(r.tree.raw)==HS::structuralKey(rr.tree.raw),"Nondeterministic direct expression");
  auto bad=r; bad.tree.leaves[0].offset+=.125;
  rejects([&] { validateBoundaryDirectTree(o,b,bad); },"Modified coefficient accepted");
  bad=r; bad.registry.supports[0].patches.clear();
  rejects([&] { validateBoundaryDirectTree(o,b,bad); },"Modified provenance accepted");
  bad=r; bad.occupied[0]=1;
  rejects([&] { validateBoundaryDirectTree(o,b,bad); },"Occupied exterior accepted");
  BoundaryDirectOptions cap; cap.signs.max_strata=1;
  check(!validateBoundaryDirectTree(o,b,r,cap).applied,"Skipped certificate claimed success");
  cap={}; cap.signs.max_sign_operations=1;
  check(!validateBoundaryDirectTree(o,b,r,cap).strict_sign_certified,"Sign op cap ignored");
  for(unsigned mode=0;mode<5;++mode)
  {
    cap={};
    if(mode==0) cap.max_canonical_cells=1; else if(mode==1) cap.max_split_checks=1;
    else if(mode==2) cap.max_nodes=1; else if(mode==3) cap.max_expanded_nodes=1; else cap.max_depth=1;
    rejects([&] { buildBoundaryDirectTree(o,b,cap); },"Construction cap ignored");
  }
}
int main()
{
  try
  {
    for(unsigned mask=1;mask<256;++mask)
    {
      std::vector<VoxelPoint> points;
      for(unsigned i=0;i<8;++i) if(mask&(1u<<i)) points.emplace_back((i&1)+.5,((i>>1)&1)+.5,((i>>2)&1)+.5);
      test(points,1,false);
    }
    test({{.5,.5,.5}}); test({{.5,.5,.5},{2.5,.5,.5}});
    test({{.5,.5,.5},{1.5,1.5,.5}}); test({{.5,.5,.5},{1.5,1.5,1.5}});
    std::vector<VoxelPoint> l,u,stair,cavity;
    for(int x=0;x<3;++x) for(int y=0;y<3;++y)
    {
      if(!x||!y) l.emplace_back(x+.5,y+.5,.5);
      if(x!=1||!y) u.emplace_back(x+.5,y+.5,.5);
      if(y<=x) stair.emplace_back(x+.5,y+.5,.5);
      for(int z=0;z<3;++z) if(x!=1||y!=1||z!=1) cavity.emplace_back(x+.5,y+.5,z+.5);
    }
    test(l); test(u); test(stair); test(cavity);
    test({{-.015,-.015,-.015},{.015,.015,.015}},.03);
    // The large empty coordinate gap must take the bounded-memory binary-search
    // fallback rather than allocating a dense canonical-axis lookup.
    test({{.5,.5,.5},{5000.5,.5,.5}},1,false);
    std::mt19937 random(9431);
    for(unsigned trial=0;trial<12;++trial)
    {
      std::vector<VoxelPoint> points;
      for(int x=-1;x<3;++x) for(int y=-1;y<3;++y) for(int z=-1;z<3;++z)
        if(random()%4==0) points.emplace_back((x+.5)*.04,(y+.5)*.04,(z+.5)*.04);
      test(points,.04);
    }
    std::cout << "Boundary direct: 255 patterns, strict all-strata/world signs, provenance and work caps passed\n";
  }
  catch(const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
