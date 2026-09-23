#include <rokae_demo/boundary_smooth_distance.hpp>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace rokae_demo
{
namespace
{
using Clock=std::chrono::steady_clock;
void require(bool v,const char* why) { if(!v) throw std::runtime_error(why); }
double ms(Clock::time_point t) { return std::chrono::duration<double,std::milli>(Clock::now()-t).count(); }
std::ofstream file(const std::filesystem::path& p)
{
  std::ofstream out(p); require(bool(out),"Cannot write smooth distance output");
  out.exceptions(std::ios::badbit|std::ios::failbit); out<<std::setprecision(17)<<std::boolalpha; return out;
}
struct Errors
{
  std::vector<double> values;
  void add(double value) { require(std::isfinite(value),"Nonfinite distance error"); values.push_back(std::abs(value)); }
  void finish(const std::string& prefix,std::map<std::string,double>& m)
  {
    m[prefix+"_count"]=values.size();
    long double sum=0,square=0; for(const auto v:values) { sum+=v; square+=static_cast<long double>(v)*v; }
    std::sort(values.begin(),values.end());
    m[prefix+"_max_m"]=values.empty() ? 0 : values.back();
    m[prefix+"_mean_m"]=values.empty() ? 0 : static_cast<double>(sum/values.size());
    m[prefix+"_rmse_m"]=values.empty() ? 0 : static_cast<double>(std::sqrt(square/values.size()));
    m[prefix+"_p95_m"]=values.empty() ? 0 : values[(95*values.size()+99)/100-1];
  }
};
struct CachedProbe
{
  HS::ValueGradient sample;
  double hard=0,smooth=0,reference=0;
  bool filled=false;
};
} // namespace

static SmoothDistanceAudit auditSmoothDistanceImpl(const SparseVoxelOccupancy& occupancy,const OctreeBoundary& boundary,
  const std::vector<VoxelPoint>& raw,const std::vector<HS::Vec3>& random_points,const HS::SmoothTree& smooth,
  const BoundaryDistanceReference& reference,const SmoothDistanceAuditOptions& options,const HS::CompiledTree& hard,
  HS::ValidationValueCache* hard_cache=nullptr)
{
  const auto start=Clock::now(); SmoothDistanceAudit result; auto& m=result.metrics;
  require(options.max_work>0,"Smooth audit work limit must be positive");
  std::size_t remaining=options.max_work;
  const auto charge=[&](std::size_t n) { require(n<=remaining,"Smooth audit work limit"); remaining-=n; };
  require(hard.structuralSharingEnabled()==options.share_subexpressions,
          "Smooth audit sharing option does not match compiled program");
  require(!hard_cache || (options.cache_probes && &hard_cache->program()==&hard),
          "Audit cache must use the same immutable program and caching must be enabled");
  const auto previous_hard_hits=hard_cache ? hard_cache->hits() : 0;
  auto hw=hard.makeWorkspace();
  auto sw=smooth.makeValueWorkspace(); auto gw=smooth.makeGradientWorkspace();
  constexpr std::size_t max_probes=200000;
  std::unordered_map<HS::ValidationPointKey,CachedProbe,HS::ValidationPointHash> cache;
  if(options.cache_probes)
    cache.reserve(std::min(max_probes,std::min(raw.size(),max_probes)+9*std::min(occupancy.occupied_voxels.size(),max_probes)
      +3*std::min(boundary.patches.size(),max_probes)+std::min(random_points.size(),max_probes)+27));
  std::size_t cache_hits=0,probe_evaluations=0,diagnostic_hits=0,diagnostic_evaluations=0;
  Errors smoothing,hard_error,sdf_error,near_error,boundary_bias,outside_error,inside_error;
  m["false_free_count"]=m["raw_points_positive"]=m["upper_order_errors"]=m["bound_errors"]=0;
  m["min_gradient_norm"]=INFINITY; m["max_gradient_norm"]=m["near_boundary_small_gradient_count"]=0;
  const auto eval=[&](const HS::Vec3& p) {
    charge(smooth.workPerQuery()); // Preserve the original conservative work cap even on hits.
    if(options.cache_probes)
    {
      const auto found=cache.find(HS::ValidationPointKey(p));
      if(found!=cache.end()) { ++diagnostic_hits; return found->second.smooth; }
    }
    ++diagnostic_evaluations;
    return smooth.evaluate(p,sw);
  };
  const auto probe=[&](const HS::Vec3& p,bool is_raw=false) {
    require(result.probes.size()<max_probes,"Smooth audit probe count limit");
    charge(hard.workPerQuery()+2*smooth.workPerQuery()+reference.patchCount());
    const HS::ValidationPointKey key(p);
    const auto found=options.cache_probes ? cache.find(key) : cache.end();
    CachedProbe values;
    if(found!=cache.end())
    {
      values=found->second; ++cache_hits;
    }
    else
    {
      ++probe_evaluations;
      values.hard=-(hard_cache ? hard_cache->evaluate(p,hw) : hard.evaluate(p,hw));
      values.smooth=smooth.evaluate(p,sw);
      values.reference=reference.signedDistance(p,occupancy); values.sample=smooth.evaluateWithGradient(p,gw);
      values.filled=occupancy.contains(p);
    }
    const double d0=values.hard, ds=values.smooth, dr=values.reference;
    const auto& sample=values.sample;
    require(ds==sample.value,"Smooth scalar/gradient value disagreement");
    const bool filled=values.filled;
    const double tolerance=64*std::numeric_limits<double>::epsilon()*smooth.sourceNodeCount()*(1+std::max(std::abs(d0),std::abs(ds)));
    m["upper_order_errors"]+=ds>d0; m["bound_errors"]+=d0-ds>smooth.errorBound()+tolerance;
    m["false_free_count"]+=ds>0 && filled; m["raw_points_positive"]+=is_raw && ds>0;
    if(ds>d0 || d0-ds>smooth.errorBound()+tolerance || (ds>0 && filled))
    {
      std::ostringstream why; why<<std::setprecision(17)<<"Smooth conservativeness failure at "<<p.x<<','<<p.y<<','<<p.z
        <<" hard="<<d0<<" smooth="<<ds<<" reference="<<dr;
      throw std::runtime_error(why.str());
    }
    require((dr>0)==!filled && (dr==0 || (d0>0)==(dr>0)),"Boundary distance reference sign mismatch");
    const auto& g=sample.gradient; const double norm=std::hypot(g.x,g.y,g.z);
    require(norm<=smooth.lipschitzBound()+1e-12*(1+smooth.lipschitzBound()),"Smooth gradient violates leaf normal bound");
    m["min_gradient_norm"]=std::min(m["min_gradient_norm"],norm); m["max_gradient_norm"]=std::max(m["max_gradient_norm"],norm);
    smoothing.add(d0-ds); hard_error.add(d0-dr); sdf_error.add(ds-dr);
    if(std::abs(dr)<=occupancy.voxel_size) { near_error.add(ds-dr); m["near_boundary_small_gradient_count"]+=norm<.1; }
    if(dr==0) boundary_bias.add(ds); else if(dr>0) outside_error.add(ds-dr); else inside_error.add(ds-dr);
    if(options.cache_probes && found==cache.end()) cache.emplace(key,values);
    result.probes.push_back({p,g,d0,ds,dr});
  };
  for(const auto& p:raw) probe({p.x(),p.y(),p.z()},true);
  for(const auto& key:occupancy.occupied_voxels)
  {
    const auto b=octreeGridBox(key,{1,1,1},occupancy.voxel_size);
    probe({b.lower[0]+(b.upper[0]-b.lower[0])/2,b.lower[1]+(b.upper[1]-b.lower[1])/2,b.lower[2]+(b.upper[2]-b.lower[2])/2});
    for(unsigned mask=0;mask<8;++mask) probe({mask&1?b.upper[0]:b.lower[0],mask&2?b.upper[1]:b.lower[1],mask&4?b.upper[2]:b.lower[2]});
  }
  const auto patchPoint=[&](const OctreeBoundaryPatch& patch,double offset) {
    const auto c=octreeBoundaryCorners(patch); double p[3]{};
    for(unsigned a=0;a<3;++a) { const double lo=static_cast<double>(c[0][a])*occupancy.voxel_size,hi=static_cast<double>(c[2][a])*occupancy.voxel_size; p[a]=lo+(hi-lo)/2; }
    p[patch.axis]+=patch.sign*offset; return HS::Vec3{p[0],p[1],p[2]};
  };
  for(const auto& p:boundary.patches)
  { probe(patchPoint(p,0)); probe(patchPoint(p,occupancy.voxel_size*.25)); probe(patchPoint(p,-occupancy.voxel_size*.25)); }
  // Deterministic off-face edge/corner probes: a small random sample may miss
  // the Euclidean-vs-plane distance discrepancy even for a single cube.
  double exterior[3][3];
  for(unsigned a=0;a<3;++a)
  {
    const auto& box=occupancy.sampling_box;
    exterior[a][0]=box.lower[a]-occupancy.voxel_size*.25;
    exterior[a][1]=box.lower[a]+(box.upper[a]-box.lower[a])/2;
    exterior[a][2]=box.upper[a]+occupancy.voxel_size*.25;
  }
  for(unsigned x=0;x<3;++x) for(unsigned y=0;y<3;++y) for(unsigned z=0;z<3;++z)
    probe({exterior[0][x],exterior[1][y],exterior[2][z]});
  for(const auto& p:random_points) probe(p);
  m["distance_probes"]=result.probes.size();
  m["probe_cache_enabled"]=options.cache_probes;
  m["probe_cache_entries"]=cache.size(); m["probe_cache_hits"]=cache_hits;
  m["probe_cache_misses"]=options.cache_probes ? probe_evaluations : 0;
  m["probe_evaluations"]=probe_evaluations;
  const double step=std::min(occupancy.voxel_size*1e-4,1e-3/smooth.beta());
  m["gradient_fd_step_m"]=step; m["gradient_checks"]=m["gradient_skipped"]=m["gradient_max_abs_error"]=0;
  const auto count=std::min(options.gradient_samples,result.probes.size());
  for(std::size_t i=0;i<count;++i)
  {
    const auto& s=result.probes[i*result.probes.size()/count]; const double p[3]{s.point.x,s.point.y,s.point.z},g[3]{s.gradient.x,s.gradient.y,s.gradient.z};
    for(unsigned a=0;a<3;++a)
    {
      double lo[3]{p[0],p[1],p[2]},hi[3]{p[0],p[1],p[2]}; lo[a]-=step; hi[a]+=step;
      if(!(step>0 && lo[a]<p[a] && hi[a]>p[a] && std::isfinite(lo[a]) && std::isfinite(hi[a]))) { ++m["gradient_skipped"]; continue; }
      const double l=eval({lo[0],lo[1],lo[2]}),h=eval({hi[0],hi[1],hi[2]});
      const double error=std::abs((h-l)/(hi[a]-lo[a])-g[a]);
      const double tolerance=3e-5*(1+smooth.lipschitzBound())+32*std::numeric_limits<double>::epsilon()*(1+std::abs(l)+std::abs(h))/(hi[a]-lo[a]);
      m["gradient_max_abs_error"]=std::max(m["gradient_max_abs_error"],error); ++m["gradient_checks"];
      require(error<=tolerance,"Smooth analytic gradient finite-difference failure");
    }
  }
  require(!count || m["gradient_checks"]>0,"All requested gradient steps unrepresentable");
  m["zero_rays_bracketed"]=m["zero_rays_blocked"]=m["zero_ray_max_offset_m"]=m["zero_ray_max_reference_distance_m"]=0;
  const auto rays=std::min(options.zero_rays,boundary.patches.size());
  for(std::size_t i=0;i<rays;++i)
  {
    SmoothZeroRay ray; ray.patch=i*boundary.patches.size()/rays; const auto& patch=boundary.patches[ray.patch];
    double lo=0,hi=0; const double initial=eval(patchPoint(patch,0)); require(initial<=0,"Positive smooth distance on true boundary");
    ray.bracketed=initial==0;
    for(unsigned step_id=1;step_id<=64 && !ray.bracketed;++step_id)
    {
      hi=2*occupancy.voxel_size*step_id/64; const auto p=patchPoint(patch,hi);
      if(occupancy.contains(p)) { ray.blocked=true; break; }
      if(eval(p)>0) ray.bracketed=true; else lo=hi;
    }
    if(ray.bracketed)
    {
      for(unsigned k=0;k<40 && hi>lo;++k) { const double mid=lo+(hi-lo)/2; if(eval(patchPoint(patch,mid))>0) hi=mid; else lo=mid; }
      ray.offset=lo+(hi-lo)/2; charge(reference.patchCount()); ray.reference_distance=reference.unsignedDistance(patchPoint(patch,ray.offset));
      ++m["zero_rays_bracketed"]; m["zero_ray_max_offset_m"]=std::max(m["zero_ray_max_offset_m"],ray.offset);
      m["zero_ray_max_reference_distance_m"]=std::max(m["zero_ray_max_reference_distance_m"],ray.reference_distance);
    }
    m["zero_rays_blocked"]+=ray.blocked; result.rays.push_back(ray);
  }
  m["zero_rays_tested"]=rays; m["zero_rays_unbracketed"]=rays-m["zero_rays_bracketed"];
  smoothing.finish("smoothing_error",m); hard_error.finish("hard_sdf_error",m); sdf_error.finish("sdf_error",m);
  near_error.finish("near_sdf_error",m); boundary_bias.finish("true_boundary_value_bias",m);
  outside_error.finish("outside_sdf_error",m); inside_error.finish("inside_sdf_error",m);
  m["diagnostic_cache_hits"]=diagnostic_hits;
  m["diagnostic_evaluations"]=diagnostic_evaluations;
  const auto hard_hits=hard_cache ? hard_cache->hits()-previous_hard_hits : 0;
  m["hard_cache_hits"]=hard_hits; m["hard_evaluations"]=probe_evaluations-hard_hits;
  m["work_used"]=options.max_work-remaining; m["audit_ms"]=ms(start); return result;
}

SmoothDistanceAudit auditSmoothDistance(const SparseVoxelOccupancy& occupancy,const OctreeBoundary& boundary,
  const std::vector<VoxelPoint>& raw,const std::vector<HS::Vec3>& random_points,const HS::SmoothTree& smooth,
  const BoundaryDistanceReference& reference,const SmoothDistanceAuditOptions& options,HS::ValidationValueCache* hard_cache)
{
  return auditSmoothDistanceImpl(occupancy,boundary,raw,random_points,smooth,reference,options,smooth.compiledProgram(),hard_cache);
}

SmoothDistanceAudit auditSmoothDistance(const SparseVoxelOccupancy& occupancy,const OctreeBoundary& boundary,
  const TreeRepresentation& tree,const std::vector<VoxelPoint>& raw,const std::vector<HS::Vec3>& random_points,
  const HS::SmoothTree& smooth,const BoundaryDistanceReference& reference,const SmoothDistanceAuditOptions& options)
{
  const auto start=Clock::now();
  const HS::CompiledTree hard(tree.raw,tree.leaves,options.share_subexpressions);
  auto result=auditSmoothDistanceImpl(occupancy,boundary,raw,random_points,smooth,reference,options,hard);
  result.metrics["audit_ms"]=ms(start);
  return result;
}

void writeSmoothDistance(const std::filesystem::path& output,const HS::SmoothTree& smooth,const SmoothDistanceAudit& audit)
{
  auto field=file(output/"smooth_distance.json");
  field<<"{\"schema\":1,\"tree\":\"free_direct_tree_raw.json\",\"mode\":\"negative_conservative_upper_lse\","
    <<"\"length_unit\":\"m\",\"beta_unit\":\"1/m\",\"distance_sign\":\"positive_free_negative_occupied_approximate_zero\","
    <<"\"beta\":"<<smooth.beta()<<",\"scalar_error_bound_m\":"<<smooth.errorBound()<<",\"path_log_arity_weight\":"<<smooth.errorWeight()
    <<",\"is_exact_euclidean_sdf\":false,\"zero_level_set_exact\":false,\"sdf_accuracy_certified\":false,"
    <<"\"analytic_real_arithmetic_upper_bound\":true,\"floating_point_interval_certified\":false,\"cbf_certified\":false,"
    <<"\"tree_and_geometry_unchanged\":true,\"reference_is_validation_only\":true}\n"; field.close();
  auto report=file(output/"distance_validation.json");
  report<<"{\"passed\":true,\"scope\":\"sampled_numerical_bounds_gradients_and_distance_errors\","
    <<"\"zero_ray_scope\":\"first_sampled_outward_bracket_within_2_voxels_not_Hausdorff\",\"metrics\":{";
  bool first=true; for(const auto& v:audit.metrics) { if(!first) report<<','; first=false; report<<'"'<<v.first<<"\":"<<v.second; }
  report<<"},\"rays\":[";
  for(std::size_t i=0;i<audit.rays.size();++i)
  { const auto& r=audit.rays[i]; report<<(i?",":"")<<"{\"patch\":"<<r.patch<<",\"bracketed\":"<<r.bracketed<<",\"blocked\":"<<r.blocked
      <<",\"offset_m\":"<<r.offset<<",\"reference_distance_m\":"<<r.reference_distance<<'}'; }
  report<<"]}\n"; report.close();
  auto probes=file(output/"distance_samples.csv"); probes<<"x_m,y_m,z_m,hard_proxy_m,smooth_proxy_m,reference_sdf_m,gradient_x,gradient_y,gradient_z\n";
  for(const auto& p:audit.probes) probes<<p.point.x<<','<<p.point.y<<','<<p.point.z<<','<<p.hard<<','<<p.smooth<<','<<p.reference<<','
    <<p.gradient.x<<','<<p.gradient.y<<','<<p.gradient.z<<'\n';
  probes.close();
}
} // namespace rokae_demo
