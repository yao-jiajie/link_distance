#include <rokae_demo/sparse_voxel_occupancy.hpp>
#include <rokae_demo/octree_boundary_tree.hpp>
#include <rokae_demo/boundary_local_tree.hpp>
#include "boundary_frontend_checks.hpp"
#include <iostream>
#include <limits>
#include <random>

using namespace rokae_demo;
void check(bool ok,const char* why) { if(!ok) throw std::runtime_error(why); }
void test(const std::vector<VoxelPoint>& points,double h)
{
  const auto v=buildSparseVoxelOccupancy(points,h), again=buildSparseVoxelOccupancy(points,h);
  const auto o=buildOctreeOccupancy(points,h,52);
  validateSparseVoxelOccupancy(v,points); validateOctreeOccupancy(o,points);
  if(!v.dense_bits.empty())
  {
    const auto cached=buildBoundaryPlaneRegistry(v);
    auto fallback=v; fallback.boundary_plane_signs={};
    const auto rebuilt=buildBoundaryPlaneRegistry(fallback);
    check(cached.coordinates==rebuilt.coordinates && cached.supports.size()==rebuilt.supports.size(),
          "Cached boundary registry shape mismatch");
    for(std::size_t i=0;i<cached.supports.size();++i)
      check(cached.supports[i].axis==rebuilt.supports[i].axis &&
            cached.supports[i].coordinate==rebuilt.supports[i].coordinate &&
            cached.supports[i].outward_signs==rebuilt.supports[i].outward_signs,
            "Cached boundary registry support mismatch");
  }
  check(v.occupied_voxels==o.occupied_voxels && v.occupied_voxels==again.occupied_voxels,"Key mismatch");
  check(v.sampling_box.lower==o.nodes[0].box.lower && v.sampling_box.upper==o.nodes[0].box.upper,"Sampling region changed");
  const auto brute=[&](const HS::Vec3& p) {
    for(const auto& k:v.occupied_voxels) if(octreeGridBox(k,{1,1,1},h).contains(p)) return true;
    return false;
  };
  const auto query=[&](const HS::Vec3& p) { check(v.contains(p)==brute(p) && v.contains(p)==o.contains(p),"Closed occupancy mismatch"); };
  for(const auto& k:v.occupied_voxels)
  {
    const auto b=octreeGridBox(k,{1,1,1},h);
    for(unsigned x=0;x<3;++x) for(unsigned y=0;y<3;++y) for(unsigned z=0;z<3;++z)
    {
      const unsigned ix[3]{x,y,z}; double q[3]{};
      for(unsigned a=0;a<3;++a) q[a]=ix[a]==0 ? b.lower[a] : ix[a]==2 ? b.upper[a] : b.lower[a]+(b.upper[a]-b.lower[a])/2;
      query({q[0],q[1],q[2]});
      // Probe one representable float on either side of every face/edge/vertex.
      for(unsigned a=0;a<3;++a) for(int sign:{-1,1})
      {
        auto p=HS::Vec3{q[0],q[1],q[2]};
        const auto value=std::nextafter(q[a],sign<0 ? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::infinity());
        if(a==0) p.x=value; else if(a==1) p.y=value; else p.z=value; query(p);
      }
    }
  }
  query({1e300,-1e300,1e300});
  const auto partition=partitionOctreeCells(o,pruneOctree(o),true);
  const auto old_boundary=extractOctreeBoundary(o,partition);
  const auto key_boundary=extractVoxelBoundary(v.occupied_voxels),new_boundary=extractSparseVoxelBoundary(v);
  boundary_frontend_checks::sameBoundary(key_boundary,new_boundary);
  validateVoxelBoundary(v.occupied_voxels,new_boundary);
  for(unsigned axis=0;axis<3;++axis)
  {
    BoundaryTreeOptions opts; opts.sweep_axis=axis;
    const auto old=buildBoundaryTree(o,partition,old_boundary,opts);
    bool fallback=false;
    const auto batch=buildVoxelBoundaryTreeBatch(h,new_boundary,opts,false,[&] { fallback=true; return old.tree; });
    const auto& fresh=batch.candidates[0];
    check(fallback==old.fallback,"Fallback policy changed");
    check(HS::structuralKey(old.tree.raw)==HS::structuralKey(fresh.tree.raw),"Raw flat expression changed");
    check(HS::structuralKey(old.tree.simplified)==HS::structuralKey(fresh.tree.simplified),"Simplified flat expression changed");
    if(!fresh.fallback)
    {
      validateVoxelBoundaryTree(v,new_boundary,fresh);
      const auto a=diagnoseBoundarySigns(o,old,old.tree), b=diagnoseBoundarySigns(v,fresh,fresh.tree);
      check(a.applied==b.applied && a.free_zeros==b.free_zeros && a.occupied_zeros==b.occupied_zeros &&
            a.closure_set_certified==b.closure_set_certified,"Sign strata changed");
    }
  }
}
int main()
{
  try
  {
    for(double h:{1.,.03,.04})
    {
      test({{-.5*h,-.5*h,-.5*h},{.5*h,.5*h,.5*h},{2.5*h,.5*h,.5*h},{2.5*h,.5*h,.5*h}},h);
      test({{0,0,0},{h,h,h}},h);
    }
    const auto limit=(std::int64_t{1}<<50)-1;
    test({{static_cast<double>(limit),0,0}},1);
    test({{static_cast<double>(-limit),0,0}},1);
    const std::vector<VoxelPoint> wide_points{
      {static_cast<double>(limit),static_cast<double>(limit),static_cast<double>(limit)},
      {static_cast<double>(-limit),static_cast<double>(-limit),static_cast<double>(-limit)}};
    const auto wide=buildSparseVoxelOccupancy(wide_points,1);
    validateSparseVoxelOccupancy(wide,wide_points);
    check(wide.occupied_voxels==std::vector<VoxelIndex>{{-limit,-limit,-limit},{limit,limit,limit}},
          "Wide-key fallback changed lexicographic occupancy");
    check(wide.dense_bits.empty(),"Wide-key occupancy unexpectedly retained a dense lookup");
    std::mt19937 rng(20260910);
    for(unsigned trial=0;trial<12;++trial)
    {
      std::vector<VoxelPoint> points;
      for(int x=-2;x<2;++x) for(int y=-2;y<2;++y) for(int z=-2;z<2;++z)
        if(rng()%4==0) points.emplace_back((x+.5)*.03,(y+.5)*.03,(z+.5)*.03);
      test(points,.03);
    }
    const auto o=buildSparseVoxelOccupancy({{.5,.5,.5}},1); bool rejected=false;
    try { o.contains({std::numeric_limits<double>::quiet_NaN(),0,0}); } catch(const std::exception&) { rejected=true; }
    check(rejected,"Nonfinite query accepted");
    auto bad=o; bad.occupied_voxels.push_back({5,5,5}); rejected=false;
    try { validateSparseVoxelOccupancy(bad,{{.5,.5,.5}}); } catch(const std::exception&) { rejected=true; }
    check(rejected,"Corrupt occupancy accepted");
    std::cout << "Sparse voxels: closed planes/edges/vertices, nextafter, index limits, identical flat geometry passed\n";
  }
  catch(const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
