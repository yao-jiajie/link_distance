#include <rokae_demo/boundary_smooth_distance.hpp>
#include <iostream>

using namespace rokae_demo;
void check(bool b,const char* why) { if(!b) throw std::runtime_error(why); }
void sameAudit(const SmoothDistanceAudit& a,const SmoothDistanceAudit& b)
{
  auto left=a.metrics,right=b.metrics;
  for(const auto* key:{"audit_ms","probe_cache_enabled","probe_cache_entries","probe_cache_hits","probe_cache_misses",
                       "probe_evaluations","diagnostic_cache_hits","diagnostic_evaluations","hard_cache_hits","hard_evaluations"})
  { left.erase(key); right.erase(key); }
  check(left==right && a.probes.size()==b.probes.size() && a.rays.size()==b.rays.size(),"Cached audit summary differs");
  for(std::size_t i=0;i<a.probes.size();++i)
  {
    const auto& x=a.probes[i]; const auto& y=b.probes[i];
    check(x.point.x==y.point.x && x.point.y==y.point.y && x.point.z==y.point.z &&
          x.gradient.x==y.gradient.x && x.gradient.y==y.gradient.y && x.gradient.z==y.gradient.z &&
          x.hard==y.hard && x.smooth==y.smooth && x.reference==y.reference,"Cached audit probe differs");
  }
  for(std::size_t i=0;i<a.rays.size();++i)
  {
    const auto& x=a.rays[i]; const auto& y=b.rays[i];
    check(x.patch==y.patch && x.bracketed==y.bracketed && x.blocked==y.blocked &&
          x.offset==y.offset && x.reference_distance==y.reference_distance,"Cached audit ray differs");
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
      const auto o=buildSparseVoxelOccupancy(points,1); validateSparseVoxelOccupancy(o,points);
      const auto b=extractVoxelBoundary(o.occupied_voxels); const auto r=buildBoundaryDirectTree(o,b);
      check(validateBoundaryDirectTree(o,b,r).strict_sign_certified,"Hard prerequisite failed");
      HS::SmoothTree smooth(r.tree.raw,r.tree.leaves,{0,.01}); auto w=smooth.makeValueWorkspace();
      HS::SmoothTree sharing(r.tree.raw,r.tree.leaves,{0,.01,true}); auto sharing_work=sharing.makeValueWorkspace();
      HS::SmoothTree binary(r.tree.raw,r.tree.leaves,{0,.01,true,true}); auto binary_work=binary.makeValueWorkspace();
      check(smooth.beta()==sharing.beta() && smooth.errorBound()==sharing.errorBound(),"Pattern CSE changed error budget");
      const BoundaryDistanceReference ref(1,b);
      for(double x:{-.5,0.,.5,1.,1.5,2.,2.5}) for(double y:{-.5,0.,.5,1.,1.5,2.,2.5}) for(double z:{-.5,0.,.5,1.,1.5,2.,2.5})
      {
        const HS::Vec3 p{x,y,z}; const double d=smooth.evaluate(p,w),hard=-HS::evaluate(r.tree.raw,r.tree.leaves,p),reference=ref.signedDistance(p,o);
        const auto cse=sharing.evaluate(p,sharing_work);
        check(cse==d && std::signbit(cse)==std::signbit(d),"Pattern CSE changed scalar/signed zero");
        const auto fast=binary.evaluate(p,binary_work);
        check(fast==d && std::signbit(fast)==std::signbit(d),"Pattern binary kernel changed scalar/signed zero");
        check(d<=hard && hard-d<=smooth.errorBound()+1e-12,"Pattern LSE bound failed");
        check(!(d>0 && o.contains(p)),"Pattern false free");
        check((reference>0)==!o.contains(p),"Reference sign failed");
      }
    }
    std::vector<VoxelPoint> points{{.5,.5,.5}}; const auto o=buildSparseVoxelOccupancy(points,1);
    const auto b=extractVoxelBoundary(o.occupied_voxels); const auto r=buildBoundaryDirectTree(o,b);
    const BoundaryDistanceReference ref(1,b);
    check(ref.signedDistance({.5,.5,.5},o)==-.5,"Inside distance incorrect");
    check(ref.signedDistance({0,.5,.5},o)==0,"Boundary distance incorrect");
    check(std::abs(ref.signedDistance({-1,-1,.5},o)-std::sqrt(2.))<1e-15,"Finite patch corner distance incorrect");
    check(ref.signedDistance({0,3,.5},o)==2,"Infinite support plane substituted for finite patch");
    HS::SmoothTree smooth(r.tree.raw,r.tree.leaves,{0,.001});
    auto audit=auditSmoothDistance(o,b,points,{{-1,-1,.5},{2,.5,.5}},smooth,ref);
    check(audit.metrics.at("hard_sdf_error_max_m")>.4,"Hard proxy distortion was hidden");
    check(audit.metrics.at("false_free_count")==0 && audit.metrics.at("gradient_checks")>0,"Smooth audit failed");
    check(audit.metrics.at("zero_rays_bracketed")==6 && audit.metrics.at("zero_ray_max_offset_m")>0,"Zero displacement not measured");
    auto legacy=auditSmoothDistance(o,b,r.tree,points,{{-1,-1,.5},{2,.5,.5}},smooth,ref);
    auto expected_metrics=audit.metrics; expected_metrics.erase("audit_ms"); legacy.metrics.erase("audit_ms");
    check(legacy.metrics==expected_metrics,"Reused and legacy audit metrics differ");
    auto program=std::make_shared<const HS::CompiledTree>(r.tree.raw,r.tree.leaves);
    HS::SmoothTree cached_smooth(program,{0,.001}); SmoothDistanceAuditOptions cached_options; cached_options.cache_probes=true;
    HS::ValidationValueCache hard_cache(program); auto hard_work=program->makeWorkspace();
    hard_cache.evaluate({.5,.5,.5},hard_work); hard_cache.evaluate({-1,-1,.5},hard_work);
    const auto cached=auditSmoothDistance(o,b,points,{{-1,-1,.5},{2,.5,.5}},cached_smooth,ref,cached_options,&hard_cache);
    sameAudit(audit,cached);
    check(cached.metrics.at("probe_cache_hits")>0 && cached.metrics.at("probe_evaluations")<cached.metrics.at("distance_probes"),
          "Repeated distance probes were not reused");
    check(cached.metrics.at("hard_cache_hits")>=2 && cached.metrics.at("diagnostic_cache_hits")>0,
          "Earlier hard values or primary smooth values were not reused");
    check(cached.metrics.at("work_used")==audit.metrics.at("work_used"),"Cache changed the conservative work charge");
    auto other_program=std::make_shared<const HS::CompiledTree>(r.tree.raw,r.tree.leaves);
    HS::ValidationValueCache wrong_cache(other_program);
    bool wrong_cache_rejected=false;
    try { auditSmoothDistance(o,b,points,{},cached_smooth,ref,cached_options,&wrong_cache); }
    catch(const std::exception&) { wrong_cache_rejected=true; }
    check(wrong_cache_rejected,"Audit accepted a cache from another program snapshot");
    auto changed=r.tree; for(auto& leaf:changed.leaves) leaf.offset+=10;
    bool changed_rejected=false;
    try { auditSmoothDistance(o,b,changed,points,{},smooth,ref); }
    catch(const std::exception&) { changed_rejected=true; }
    check(changed_rejected,"Legacy audit ignored the supplied tree coefficients");
    bool limited=false; try { SmoothDistanceAuditOptions cap; cap.max_work=1; auditSmoothDistance(o,b,points,{},smooth,ref,cap); }
    catch(const std::exception&) { limited=true; } check(limited,"Smooth work cap ignored");
    bool mismatched=false; try { SmoothDistanceAuditOptions mismatch; mismatch.share_subexpressions=true;
      auditSmoothDistance(o,b,points,{},smooth,ref,mismatch); }
    catch(const std::exception&) { mismatched=true; } check(mismatched,"Smooth audit accepted a mismatched program mode");
    std::cout<<"Boundary smooth distance: 255 patterns, finite faces, signs, gradients, zero rays and work limits passed\n";
  }
  catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
