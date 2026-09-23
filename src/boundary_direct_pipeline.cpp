#include <rokae_demo/boundary_direct_tree.hpp>
#include <rokae_demo/halfspace_compiled_tree.hpp>
#include <rokae_demo/boundary_smooth_distance.hpp>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>

namespace rokae_demo
{
namespace
{
namespace fs=std::filesystem;
using Clock=std::chrono::steady_clock;
using Metrics=std::map<std::string,double>;
enum class ValidationCacheMode { Auto,None,Exact };
double ms(Clock::time_point t) { return std::chrono::duration<double,std::milli>(Clock::now()-t).count(); }
void require(bool ok,const char* why) { if(!ok) throw std::runtime_error(why); }
struct Options
{
  fs::path input,output; double h=0; std::size_t samples=1000;
  std::size_t max_probe_ops=1000000000, max_query_ops=8000000000ULL;
  bool query=false, dag=true, cross_check=false; BoundaryDirectOptions work;
  bool share_subexpressions=true; ValidationCacheMode validation_cache=ValidationCacheMode::Auto;
  bool smooth=false; HS::SmoothTreeOptions smoothing{0,0,false,true}; SmoothDistanceAuditOptions distance_audit;
};
std::size_t integer(const std::string& value)
{
  require(!value.empty() && value.find_first_not_of("0123456789")==std::string::npos,"Expected nonnegative integer");
  const auto n=std::stoull(value); require(n<=std::numeric_limits<std::size_t>::max(),"Integer overflow"); return n;
}
Options parse(int argc,char** argv)
{
  require(argc>=3,"Expected input.xyz and output folder"); Options o; o.input=argv[1]; o.output=argv[2];
  std::set<std::string> seen;
  for(int i=3;i<argc;i+=2)
  {
    const std::string flag=argv[i]; require(i+1<argc && seen.insert(flag).second,"Missing/repeated direct option");
    const std::string v=argv[i+1];
    if(flag=="--pipeline") require(v=="octree","Expected --pipeline octree routing name");
    else if(flag=="--tree-source") require(v=="boundary","Expected --tree-source boundary");
    else if(flag=="--boundary-expression") require(v=="direct","Expected direct expression");
    else if(flag=="--boundary-occupancy") require(v=="voxels","Direct requires sparse voxel occupancy, not Octree");
    else if(flag=="--boundary-diagnostics") require(v=="true","Direct sign certification cannot be disabled");
    else if(flag=="--boundary-query")
    { require(v=="dag" || v=="recursive","boundary-query must be dag/recursive"); o.dag=v=="dag"; }
    else if(flag=="--boundary-query-cross-check")
    { require(v=="true" || v=="false","boundary-query-cross-check must be true/false"); o.cross_check=v=="true"; }
    else if(flag=="--execution-sharing")
    { require(v=="none" || v=="structural","execution-sharing must be none/structural"); o.share_subexpressions=v=="structural"; }
    else if(flag=="--validation-cache")
    {
      require(v=="auto" || v=="none" || v=="exact","validation-cache must be auto/none/exact");
      o.validation_cache=v=="auto" ? ValidationCacheMode::Auto : v=="exact" ? ValidationCacheMode::Exact : ValidationCacheMode::None;
    }
    else if(flag=="--sign-propagation")
    { require(v=="scalar" || v=="packed","sign-propagation must be scalar/packed"); o.work.signs.packed_sign_propagation=v=="packed"; }
    else if(flag=="--distance-field")
    { require(v=="exact" || v=="lse","distance-field must be exact/lse (exact preserves the existing tree only)"); o.smooth=v=="lse"; }
    else if(flag=="--lse-kernel")
    { require(v=="generic" || v=="binary","lse-kernel must be generic/binary"); o.smoothing.binary_kernel=v=="binary"; }
    else if(flag=="--lse-beta" || flag=="--lse-error")
    {
      std::size_t used=0; const double n=std::stod(v,&used);
      require(used==v.size() && std::isfinite(n) && n>0,"LSE parameter must be finite positive");
      if(flag=="--lse-beta") o.smoothing.beta=n; else o.smoothing.max_error=n;
    }
    else if(flag=="--sdf-max-work") o.distance_audit.max_work=integer(v);
    else if(flag=="--sdf-gradient-samples") o.distance_audit.gradient_samples=integer(v);
    else if(flag=="--sdf-zero-rays") o.distance_audit.zero_rays=integer(v);
    else if(flag=="--voxel-size")
    { std::size_t used=0; o.h=std::stod(v,&used); require(used==v.size() && std::isfinite(o.h) && o.h>0,"voxel-size must be finite positive meters"); }
    else if(flag=="--validation-samples") o.samples=integer(v);
    else if(flag=="--boundary-max-cells") o.work.max_canonical_cells=integer(v);
    else if(flag=="--boundary-max-direct-nodes") o.work.max_nodes=integer(v);
    else if(flag=="--boundary-max-expanded-nodes") o.work.max_expanded_nodes=integer(v);
    else if(flag=="--boundary-max-split-checks") o.work.max_split_checks=integer(v);
    else if(flag=="--boundary-max-direct-depth")
    { const auto n=integer(v); require(n<=256,"Direct depth exceeds 256"); o.work.max_depth=static_cast<unsigned>(n); }
    else if(flag=="--boundary-max-strata") o.work.signs.max_strata=integer(v);
    else if(flag=="--boundary-max-sign-ops") o.work.signs.max_sign_operations=integer(v);
    else if(flag=="--boundary-max-probe-ops") o.max_probe_ops=integer(v);
    else if(flag=="--boundary-max-query-ops") o.max_query_ops=integer(v);
    else if(flag=="--query-benchmark")
    { require(v=="true" || v=="false","query-benchmark must be true/false"); o.query=v=="true"; }
    else throw std::invalid_argument("Unsupported direct option (no sweep/prisms/Octree/fallback): "+flag);
  }
  require(o.h>0 && o.max_probe_ops && o.max_query_ops && o.work.signs.max_strata && o.work.signs.max_sign_operations,
          "Explicit positive voxel-size and positive validation work limits required");
  if(o.smooth)
    require((o.smoothing.beta>0)!=(o.smoothing.max_error>0) && o.distance_audit.max_work>0,
            "LSE requires exactly one of --lse-beta / --lse-error and a positive audit work limit");
  else require(!seen.count("--lse-beta") && !seen.count("--lse-error") && !seen.count("--lse-kernel") && !seen.count("--sdf-max-work") &&
               !seen.count("--sdf-gradient-samples") && !seen.count("--sdf-zero-rays") && !seen.count("--validation-cache"),
               "LSE options require --distance-field lse");
  o.smoothing.share_subexpressions=o.distance_audit.share_subexpressions=o.share_subexpressions;
  return o;
}
std::vector<VoxelPoint> readCloud(const fs::path& path)
{
  std::ifstream in(path); require(bool(in),"Cannot read XYZ"); std::vector<VoxelPoint> points; std::string line;
  while(std::getline(in,line))
  {
    const auto comment=line.find('#'); if(comment!=std::string::npos) line.resize(comment);
    std::istringstream row(line); row>>std::ws; if(row.eof()) continue;
    double x,y,z; require(bool(row>>x>>y>>z) && std::isfinite(x) && std::isfinite(y) && std::isfinite(z),"Invalid XYZ point");
    row>>std::ws; require(row.eof(),"XYZ requires exactly three columns"); points.emplace_back(x,y,z);
  }
  require(!in.bad() && !points.empty(),"Empty/unreadable XYZ"); return points;
}
std::ofstream file(const fs::path& path)
{
  std::ofstream out(path); require(bool(out),"Cannot write direct output"); out.exceptions(std::ios::badbit|std::ios::failbit);
  out<<std::setprecision(17)<<std::boolalpha; return out;
}
template<class Range> void array(std::ostream& out,const Range& v)
{ out<<'['; bool first=true; for(const auto& x:v) { if(!first) out<<','; first=false; out<<x; } out<<']'; }
std::vector<HS::Vec3> randomPoints(const SparseVoxelOccupancy& o,std::size_t count)
{
  std::mt19937_64 generator(20260909); std::uniform_real_distribution<double> unit(-.05,1.05);
  const auto& b=o.sampling_box; std::vector<HS::Vec3> points; points.reserve(count);
  for(std::size_t i=0;i<count;++i)
    points.push_back({b.lower[0]+unit(generator)*(b.upper[0]-b.lower[0]),b.lower[1]+unit(generator)*(b.upper[1]-b.lower[1]),
                      b.lower[2]+unit(generator)*(b.upper[2]-b.lower[2])});
  return points;
}
template<class Function> double queryMedian(const std::vector<HS::Vec3>& points,Function f)
{
  std::vector<double> times; volatile double sink=0;
  for(const auto& p:points) sink=f(p);
  for(unsigned repeat=0;repeat<7;++repeat)
  {
    double sum=0; const auto start=Clock::now();
    for(unsigned round=0;round<5;++round) for(const auto& p:points) sum+=f(p);
    times.push_back(ms(start)*1000/(5*points.size())); sink=sum;
  }
  (void)sink; std::sort(times.begin(),times.end()); return times[3];
}
void geometryOutputs(const Options& o,const SparseVoxelOccupancy& occupancy,const OctreeBoundary& boundary,
                     const BoundaryDirectResult& r,const BoundarySignDiagnostics& d,const Metrics& m)
{
  writeTreeJson(o.output/"free_direct_tree_raw.json",r.tree,r.tree.raw,"raw","boundary-direct");
  writeTreeJson(o.output/"free_direct_tree_simplified.json",r.tree,r.tree.simplified,"simplified","boundary-direct");
  writeVoxelBoundary(o.output,o.h,boundary);
  writeBoundaryDirectDiagnostics(o.output/"boundary_sign_diagnostics.json",r,d);
  auto registry=file(o.output/"boundary_supports.json"); registry<<"{\"length_unit\":\"m\",\"supports\":[";
  for(std::size_t i=0;i<r.registry.supports.size();++i)
  {
    const auto& s=r.registry.supports[i]; registry<<(i ? "," : "")<<"{\"id\":"<<i<<",\"axis\":"<<s.axis
      <<",\"coordinate_grid\":"<<s.coordinate<<",\"coordinate_m\":"<<static_cast<double>(s.coordinate)*o.h<<",\"patches\":";
    array(registry,s.patches); registry<<'}';
  }
  registry<<"],\"leaf_sources\":[";
  for(std::size_t i=0;i<r.leaf_sources.size();++i)
    registry<<(i ? "," : "")<<"{\"leaf_id\":"<<i<<",\"support\":"<<r.leaf_sources[i].first<<",\"sign\":"<<r.leaf_sources[i].second<<'}';
  registry<<"]}\n"; registry.close();
  auto splits=file(o.output/"direct_splits.json");
  splits<<"{\"rule\":\"min(mixed_children,max_child_canonical_cells,axis,cut_rank)\",\"includes_unbounded_cells\":true,\"dimensions\":";
  array(splits,r.dimensions); splits<<",\"splits\":[";
  for(std::size_t i=0;i<r.splits.size();++i)
  {
    const auto& s=r.splits[i]; splits<<(i ? "," : "")<<"{\"lower\":"; array(splits,s.lower); splits<<",\"upper\":"; array(splits,s.upper);
    splits<<",\"support\":"<<s.support<<",\"cut_rank\":"<<s.cut<<",\"left_occupied\":"<<s.left_occupied<<",\"right_occupied\":"<<s.right_occupied<<'}';
  }
  splits<<"]}\n"; splits.close();
  auto reference=file(o.output/"occupancy_reference.json");
  reference<<"{\"schema\":2,\"debug_only\":true,\"backend\":\"voxels\",\"voxel_size\":"<<o.h<<",\"occupied_keys\":[";
  for(std::size_t i=0;i<occupancy.occupied_voxels.size();++i) { if(i) reference<<','; array(reference,occupancy.occupied_voxels[i]); }
  reference<<"]}\n"; reference.close();
  auto validation=file(o.output/"validation_report.json");
  validation<<"{\"passed\":true,\"strict_sign_certified\":true,\"standalone_tree_is_strict_complement\":true,"
    <<"\"zero_level_set\":\"boundary_of_closed_voxel_union\",\"certificate\":\"all_support_plane_arrangement_strata_including_unbounded\","
    <<"\"tested_points\":"<<m.at("tested_points")<<",\"raw_points_classified_free\":"<<m.at("raw_points_classified_free")
    <<",\"world_sign_errors\":"<<m.at("world_sign_errors")
    <<",\"scalar_cross_checks\":"<<m.at("scalar_cross_checks")<<",\"scalar_cross_check_errors\":"<<m.at("scalar_cross_check_errors")
    <<",\"query_backend\":\""<<(o.dag ? "dag" : "recursive")<<"\",\"failures\":[]}\n";
  validation.close();
}
} // namespace

int runBoundaryDirectPipeline(int argc,char** argv)
{
  try
  {
    auto o=parse(argc,argv); Metrics m; auto start=Clock::now(); const auto raw=readCloud(o.input); m["input_io_ms"]=ms(start);
    const auto core=Clock::now(); start=Clock::now(); const auto occupancy=buildSparseVoxelOccupancy(raw,o.h); m["occupancy_ms"]=ms(start);
    const auto boundary=extractSparseVoxelBoundary(occupancy);
    const auto r=buildBoundaryDirectTree(occupancy,boundary,o.work);
    m["geometry_core_ms"]=ms(core); start=Clock::now();
    std::shared_ptr<const HS::CompiledTree> program; std::vector<double> workspace;
    if(o.dag || o.cross_check || o.smooth)
    {
      program=std::make_shared<const HS::CompiledTree>(
        r.serialized,r.tree.leaves,o.share_subexpressions);
      if(o.dag || o.cross_check) workspace=program->makeWorkspace();
    }
    m["evaluator_compile_ms"]=program ? ms(start) : 0;
    std::optional<HS::SmoothTree> smooth;
    std::vector<double> smooth_workspace; std::vector<HS::ValueGradient> gradient_workspace;
    start=Clock::now();
    if(o.smooth)
    {
      smooth.emplace(program,o.smoothing);
      require(&smooth->compiledProgram()==program.get(),"Smooth tree did not retain the shared execution program");
      smooth_workspace=smooth->makeValueWorkspace(); gradient_workspace=smooth->makeGradientWorkspace();
    }
    if(smooth) m["smooth_compile_ms"]=ms(start);
    constexpr std::size_t automatic_cache_min_work=256;
    const bool use_validation_cache=smooth && (o.validation_cache==ValidationCacheMode::Exact ||
      (o.validation_cache==ValidationCacheMode::Auto && smooth->workPerQuery()>=automatic_cache_min_work));
    o.distance_audit.cache_probes=use_validation_cache;
    m["validation_cache_requested_auto"]=o.validation_cache==ValidationCacheMode::Auto;
    m["validation_cache_selected"]=use_validation_cache;
    m["validation_cache_auto_min_work"]=automatic_cache_min_work;
    const std::size_t query_work=o.dag ? program->workPerQuery() : r.expanded_nodes;
    require(!o.cross_check || program->workPerQuery()<=std::numeric_limits<std::size_t>::max()-r.expanded_nodes,
            "Direct cross-check work estimate overflow");
    const std::size_t probe_work=o.cross_check ? program->workPerQuery()+r.expanded_nodes : query_work;
    const auto value=[&](const HS::Vec3& p) { return o.dag ? program->evaluate(p,workspace) : HS::evaluate(r.tree.raw,r.tree.leaves,p); };
    m["core_runtime_ms"]=m["core_build_time"]=ms(core);
    start=Clock::now(); validateSparseVoxelOccupancy(occupancy,raw); const auto d=validateBoundaryDirectTree(occupancy,boundary,r,o.work);
    m["certificate_ms"]=ms(start);
    if(!d.applied || !d.strict_sign_certified)
    {
      fs::create_directories(o.output); writeBoundaryDirectDiagnostics(o.output/"boundary_sign_failure.json",r,d);
      throw std::runtime_error("Direct strict sign certificate failed/skipped; no new tree published. See boundary_sign_failure.json");
    }
    start=Clock::now(); std::size_t budget=o.max_probe_ops/probe_work;
    std::optional<HS::ValidationValueCache> validation_cache;
    if(use_validation_cache && o.dag) validation_cache.emplace(program);
    m["tested_points"]=m["world_sign_errors"]=m["raw_points_classified_free"]=0;
    m["scalar_cross_checks"]=m["scalar_cross_check_errors"]=0;
    const auto check=[&](const HS::Vec3& p,bool is_raw=false) {
      require(budget>0,"Direct world probe work limit; no new tree published"); --budget;
      const double q=validation_cache ? validation_cache->evaluate(p,workspace) : value(p);
      const int actual=(q>0)-(q<0);
      if(o.cross_check)
      {
        const double reference=o.dag ? HS::evaluate(r.tree.raw,r.tree.leaves,p) : program->evaluate(p,workspace);
        ++m["scalar_cross_checks"];
        if(q!=reference || std::signbit(q)!=std::signbit(reference))
        { ++m["scalar_cross_check_errors"]; throw std::runtime_error("Direct compiled/recursive scalar or signed-zero mismatch"); }
      }
      const int expected=boundaryDirectExpectedSign(r,o.h,p); const bool occupied=occupancy.contains(p);
      ++m["tested_points"]; m["raw_points_classified_free"]+=is_raw && q<0;
      if(!std::isfinite(q) || actual!=expected || occupied!=(expected>=0) || (is_raw && !occupied))
      {
        ++m["world_sign_errors"]; std::ostringstream why;
        why<<std::setprecision(17)<<"Direct world sign failure at ("<<p.x<<','<<p.y<<','<<p.z<<"): expected "<<expected<<", q="<<q;
        throw std::runtime_error(why.str());
      }
    };
    require(raw.size()<=budget && o.samples<=budget-raw.size(),"Direct requested probe work limit");
    for(const auto& p:raw) check({p.x(),p.y(),p.z()},true);
    for(const auto& key:occupancy.occupied_voxels)
    {
      const auto b=octreeGridBox(key,{1,1,1},o.h);
      check({b.lower[0]+(b.upper[0]-b.lower[0])/2,b.lower[1]+(b.upper[1]-b.lower[1])/2,b.lower[2]+(b.upper[2]-b.lower[2])/2});
      for(unsigned mask=0;mask<8;++mask) check({(mask&1)?b.upper[0]:b.lower[0],(mask&2)?b.upper[1]:b.lower[1],(mask&4)?b.upper[2]:b.lower[2]});
    }
    for(const auto& patch:boundary.patches)
    {
      const auto corners=octreeBoundaryCorners(patch); double p[3]{};
      for(unsigned a=0;a<3;++a) { const double lo=static_cast<double>(corners[0][a])*o.h,hi=static_cast<double>(corners[2][a])*o.h; p[a]=lo+(hi-lo)/2; }
      check({p[0],p[1],p[2]}); p[patch.axis]-=.25*o.h; check({p[0],p[1],p[2]}); p[patch.axis]+=.5*o.h; check({p[0],p[1],p[2]});
    }
    for(const auto& p:randomPoints(occupancy,o.samples)) check(p);
    m["world_probe_cache_hits"]=validation_cache ? validation_cache->hits() : 0;
    m["world_probe_evaluations"]=validation_cache ? validation_cache->evaluations() : m["tested_points"];
    m["probe_validation_ms"]=ms(start); m["validation_ms"]=m["certificate_ms"]+m["probe_validation_ms"]; m["total_ms"]=ms(core);
    std::optional<BoundaryDistanceReference> distance_reference; std::optional<SmoothDistanceAudit> distance_audit;
    if(smooth)
    {
      start=Clock::now(); distance_reference.emplace(o.h,boundary);
      m["sdf_reference_prepare_ms"]=ms(start); start=Clock::now();
      distance_audit=auditSmoothDistance(occupancy,boundary,raw,randomPoints(occupancy,o.samples),*smooth,*distance_reference,
                                        o.distance_audit,validation_cache ? &*validation_cache : nullptr);
      m["sdf_validation_ms"]=ms(start)+m["sdf_reference_prepare_ms"];
      for(const auto& item:distance_audit->metrics) m["distance_"+item.first]=item.second;
      m["exact_validation_ms"]=m["validation_ms"]; m["validation_ms"]+=m["sdf_validation_ms"]; m["total_ms"]=ms(core);
      m["lse_beta"]=smooth->beta(); m["lse_error_bound_m"]=smooth->errorBound(); m["lse_path_weight"]=smooth->errorWeight();
      m["smooth_program_bytes"]=smooth->programBytes();
      m["smooth_execution_nodes"]=smooth->nodeCount(); m["smooth_shared_nodes"]=smooth->sharedNodeCount();
      m["smooth_execution_edges"]=smooth->edgeCount();
      m["smooth_exp_calls_per_query"]=smooth->expCallsPerQuery(); m["smooth_log_calls_per_query"]=smooth->logCallsPerQuery();
      m["smooth_binary_kernel"]=smooth->binaryKernelEnabled(); m["smooth_binary_nodes"]=smooth->binaryNodeCount();
      m["smooth_generic_nodes"]=smooth->logCallsPerQuery()-smooth->binaryNodeCount();
      m["smooth_reuses_compiled_program"]=1;
      m["smooth_value_workspace_bytes"]=smooth_workspace.capacity()*sizeof(double);
      m["smooth_gradient_workspace_bytes"]=gradient_workspace.capacity()*sizeof(HS::ValueGradient);
      m["validation_cache_entries"]=validation_cache ? validation_cache->size() : 0;
    }
    m["raw_points"]=raw.size(); m["occupied_voxels"]=occupancy.occupied_voxels.size(); m["voxel_size"]=o.h;
    m["octree_nodes"]=m["free_rectangles"]=m["free_prisms"]=m["exterior_free_branches"]=m["fallback"]=m["simplifier_applied"]=0;
    m["raw_exposed_faces"]=boundary.faces.size(); m["merged_boundary_patches"]=boundary.patches.size();
    m["boundary_extraction_ms"]=boundary.extraction_ms; m["boundary_merge_ms"]=boundary.merge_ms;
    m["registry_ms"]=r.registry_ms; m["canonical_ms"]=r.canonical_ms; m["tree_ms"]=r.tree_ms;
    m["canonical_cells"]=r.occupied.size(); m["split_count"]=r.splits.size(); m["split_checks"]=r.split_checks;
    m["allocated_tree_nodes"]=r.allocated_nodes;
    m["tree_nodes"]=m["raw_tree_nodes"]=m["constructed_tree_nodes"]=m["simplified_tree_nodes"]=r.constructed_nodes;
    m["expanded_node_references"]=r.expanded_nodes; m["plane_references"]=r.expanded_plane_refs;
    m["query_node_visits"]=o.dag ? program->nodeCount() : r.expanded_nodes;
    m["query_work_estimate"]=query_work; m["query_cross_check_enabled"]=o.cross_check;
    m["execution_structural_sharing"]=o.share_subexpressions;
    m["compiled_shared_nodes"]=program ? program->sharedNodeCount() : 0;
    m["compiled_program_bytes"]=program ? program->programBytes() : 0;
    m["execution_program_instances"]=program ? 1 : 0;
    m["query_workspace_bytes"]=workspace.capacity()*sizeof(double);
    m["unique_halfspaces"]=r.tree.leaves.size(); m["recursion_depth"]=r.recursion_depth; m["expression_depth"]=r.expression_depth;
    m["phi_full"]=r.phi_full; m["phi_folded"]=r.phi_folded;
    m["artificial_zero_faces"]=d.free_zeros[2]; m["artificial_zero_edges"]=d.free_zeros[1]; m["artificial_zero_vertices"]=d.free_zeros[0];
    m["occupied_zero_faces"]=d.occupied_zeros[2]; m["occupied_zero_edges"]=d.occupied_zeros[1]; m["occupied_zero_vertices"]=d.occupied_zeros[0];
    m["strata_count"]=d.strata_count; m["true_boundary_sign_errors"]=d.boundary_sign_errors;
    m["sign_packed_propagation"]=o.work.signs.packed_sign_propagation;
    m["sign_propagation_batches"]=o.work.signs.packed_sign_propagation ? (d.strata_count+63)/64 : d.strata_count;
    m["strict_sign_certified"]=m["sign_certificate_applied"]=m["passed"]=1; m["query_benchmark_enabled"]=o.query;
    m["query_zero_guard_fraction"]=0;
    fs::create_directories(o.output);
    if(o.query)
    {
      const auto queries=randomPoints(occupancy,1024); std::vector<HS::Vec3> boundary_points;
      for(const auto& patch:boundary.patches) for(const auto& key:octreeBoundaryCorners(patch))
        if(boundary_points.size()<1024) boundary_points.push_back({static_cast<double>(key[0])*o.h,static_cast<double>(key[1])*o.h,static_cast<double>(key[2])*o.h});
      require(queries.size()+boundary_points.size()<=o.max_query_ops/query_work/36,"Direct query benchmark work limit");
      const auto strict=[&](const auto& p) { return value(p)<0; };
      start=Clock::now(); m["strict_free_query_us"]=m["query_time"]=queryMedian(queries,strict);
      m["boundary_strict_free_query_us"]=queryMedian(boundary_points,strict); m["query_benchmark_ms"]=ms(start);
      if(smooth)
      {
        const auto extra_work=2*smooth->workPerQuery()+distance_reference->patchCount();
        const auto remaining_query_work=o.max_query_ops/36-(queries.size()+boundary_points.size())*query_work;
        require(queries.size()<=remaining_query_work/extra_work,"Smooth query benchmark work limit");
        start=Clock::now();
        m["smooth_value_query_us"]=queryMedian(queries,[&](const auto& p) { return smooth->evaluate(p,smooth_workspace); });
        m["smooth_gradient_query_us"]=queryMedian(queries,[&](const auto& p) {
          const auto s=smooth->evaluateWithGradient(p,gradient_workspace); return s.value+s.gradient.x+s.gradient.y+s.gradient.z;
        });
        m["reference_sdf_query_us"]=queryMedian(queries,[&](const auto& p) { return distance_reference->signedDistance(p,occupancy); });
        m["smooth_query_benchmark_ms"]=ms(start);
      }
      auto samples=file(o.output/"query_samples.json"); samples<<'{';
      const auto group=[&](const char* name,const auto& ps,bool comma) {
        samples<<(comma?",":"")<<'"'<<name<<"\":[";
        for(std::size_t i=0;i<ps.size();++i) samples<<(i?",":"")<<'['<<ps[i].x<<','<<ps[i].y<<','<<ps[i].z<<']';
        samples<<']';
      };
      group("random",queries,false); group("boundary",boundary_points,true); samples<<",\"artificial_seams\":[]}\n"; samples.close();
    }
    start=Clock::now(); geometryOutputs(o,occupancy,boundary,r,d,m);
    if(smooth) writeSmoothDistance(o.output,*smooth,*distance_audit);
    m["output_io_ms"]=ms(start);
    auto metrics=file(o.output/"benchmark.json");
    metrics<<"[{\"pipeline\":\"boundary-direct\",\"expression_mode\":\"direct\",\"occupancy_backend\":\"voxels\","
      <<"\"query_backend\":\""<<(o.dag ? "dag" : "recursive")<<"\",\"length_unit\":\"m\",\"time_unit\":\"ms\",\"query_time_unit\":\"us\"";
    for(const auto& item:m) metrics<<",\n\""<<item.first<<"\":"<<item.second;
    metrics<<"}]\n"; metrics.close();
    // Publish entry point LAST, only after all certification and export succeed.
    auto manifest=file(o.output/"strict_free_space.json");
    manifest<<"{\"schema\":2,\"length_unit\":\"m\",\"target\":\"R3_minus_closed_occupied_voxel_union\","
      <<"\"tree_set\":\"closure_of_nonoccupied_domain\",\"classifier\":\"root_lt_zero\",\"requires_zero_guard\":false,"
      <<"\"strict_sign_certified\":true,\"zero_level_set\":\"boundary_V\",\"fallback\":false,\"simplifier_applied\":false,"
      <<"\"tree\":\"free_direct_tree_simplified.json\",\"raw_tree\":\"free_direct_tree_raw.json\",\"is_signed_distance\":false,"
      <<"\"query_backend\":\""<<(o.dag ? "dag" : "recursive")<<"\",\"occupancy_backend\":\"voxels\",\"occupancy_reference\":\"occupancy_reference.json\",\"occupancy_reference_required_for_query\":false}\n";
    manifest.close();
    std::cout<<"Boundary direct: certified q<0 free, q=0 true voxel boundary, q>0 occupied interior; root-only query.\n";
    if(smooth) std::cout<<"Smooth distance: positive-free upper-LSE proxy with analytic gradient; shifted zero set, NOT exact Euclidean SDF.\n";
    for(const auto& item:m) std::cout<<item.first<<": "<<item.second<<'\n';
    return 0;
  }
  catch(const std::exception& e) { std::cerr<<"Boundary direct error: "<<e.what()<<'\n'; return 1; }
}
} // namespace rokae_demo
