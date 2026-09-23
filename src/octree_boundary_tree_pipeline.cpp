#include <rokae_demo/octree_boundary_tree.hpp>
#include <rokae_demo/boundary_local_tree.hpp>
#include <rokae_demo/boundary_direct_tree.hpp>

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <numeric>
#include <random>
#include <sstream>

namespace rokae_demo
{
namespace
{
namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;
using Metrics = std::map<std::string,double>;
double ms(Clock::time_point start) { return std::chrono::duration<double,std::milli>(Clock::now()-start).count(); }
void require(bool v, const char* text) { if(!v) throw std::runtime_error(text); }
struct Options
{
  fs::path input, output;
  double h = 0; unsigned depth = 20; std::size_t samples = 1000;
  bool query_benchmark = false;
  bool best_axis = false;
  std::string requested_sweep = "x";
  bool diagnostics = false;
  std::string expression = "flat";
  std::string occupancy = "octree";
  bool voxel_validation_reference = false;
  BoundaryTreeOptions work;
  BoundaryLocalOptions local;
};
std::size_t integer(const std::string& s)
{
  require(!s.empty() && s.find_first_not_of("0123456789") == std::string::npos,"Expected nonnegative integer");
  const auto n = std::stoull(s);
  require(n <= std::numeric_limits<std::size_t>::max(),"Integer overflow"); return static_cast<std::size_t>(n);
}
Options parse(int argc, char** argv)
{
  require(argc >= 3,"Expected input.xyz and output folder");
  Options o; o.input = argv[1]; o.output = argv[2]; std::set<std::string> seen;
  for(int i = 3; i < argc; i += 2)
  {
    const std::string flag = argv[i];
    require(i+1 < argc && seen.insert(flag).second,"Missing/repeated boundary option");
    const std::string value = argv[i+1];
    if(flag == "--pipeline") require(value == "octree","Expected --pipeline octree");
    else if(flag == "--tree-source") require(value == "boundary","Expected --tree-source boundary");
    else if(flag == "--boundary-sweep")
    {
      require(value == "x" || value == "y" || value == "z" || value == "best-axis","boundary-sweep must be x/y/z/best-axis");
      o.requested_sweep = value; o.best_axis = value == "best-axis";
      o.work.sweep_axis = value == "y" ? 1 : value == "z" ? 2 : 0;
    }
    else if(flag == "--octree-pruning") require(value == "rect","Boundary prototype uses the rect baseline");
    else if(flag == "--boundary-expression")
    { require(value == "flat" || value == "local","boundary-expression must be flat/local"); o.expression = value; }
    else if(flag == "--boundary-occupancy")
    { require(value == "octree" || value == "voxels","boundary-occupancy must be octree/voxels"); o.occupancy = value; }
    else if(flag == "--boundary-validation-reference")
    { require(value == "native" || value == "voxels","boundary-validation-reference must be native/voxels"); o.voxel_validation_reference = value == "voxels"; }
    else if(flag == "--boundary-diagnostics")
    { require(value == "true" || value == "false","boundary-diagnostics must be true/false"); o.diagnostics = value == "true"; }
    else if(flag == "--boundary-max-local-nodes") o.local.max_expanded_nodes = integer(value);
    else if(flag == "--boundary-max-local-attempts") o.local.max_attempts = integer(value);
    else if(flag == "--boundary-max-adjacency-checks") o.local.max_adjacency_checks = integer(value);
    else if(flag == "--boundary-max-strata") o.local.max_strata = integer(value);
    else if(flag == "--boundary-max-sign-ops") o.local.max_sign_operations = integer(value);
    else if(flag == "--sign-propagation")
    {
      require(value == "scalar" || value == "packed","sign-propagation must be scalar/packed");
      o.local.packed_sign_propagation = value == "packed";
    }
    else if(flag == "--octree-approx") require(value == "none","Boundary sweep cannot be mixed with hull approximation");
    else if(flag == "--octree-boundary") require(value == "exact","Boundary sweep requires exact exposed boundaries");
    else if(flag == "--voxel-size")
    {
      std::size_t used = 0; o.h = std::stod(value,&used);
      require(used == value.size() && std::isfinite(o.h) && o.h > 0,"voxel-size must be finite positive meters");
    }
    else if(flag == "--max-depth")
    { const auto d = integer(value); require(d <= 52,"max-depth exceeds 52"); o.depth = static_cast<unsigned>(d); }
    else if(flag == "--validation-samples") o.samples = integer(value);
    else if(flag == "--boundary-max-cells") o.work.max_compressed_cells = integer(value);
    else if(flag == "--boundary-max-events") o.work.max_event_updates = integer(value);
    else if(flag == "--boundary-max-prisms") o.work.max_free_prisms = integer(value);
    else if(flag == "--boundary-merge-passes")
    { const auto n = integer(value); require(n <= 100,"Too many merge passes"); o.work.max_merge_passes = static_cast<unsigned>(n); }
    else if(flag == "--query-benchmark")
    { require(value == "true" || value == "false","query-benchmark must be true/false"); o.query_benchmark = value == "true"; }
    else throw std::invalid_argument("Unknown boundary option: "+flag);
  }
  require(o.h > 0,"Explicit --voxel-size is required");
  require(o.work.max_compressed_cells && o.work.max_event_updates && o.work.max_free_prisms && o.work.max_merge_passes,
          "Boundary work limits must be positive");
  require(o.local.max_expanded_nodes && o.local.max_expanded_nodes < std::numeric_limits<std::size_t>::max() && o.local.max_attempts &&
          o.local.max_adjacency_checks && o.local.max_strata && o.local.max_sign_operations,"Local/sign work limits must be positive");
  if(o.expression == "local")
  {
    require(o.occupancy == "octree","Voxel occupancy comparison currently supports flat only");
    require(!o.best_axis,"Local prototype requires a fixed x/y/z sweep; best-axis ranks flat candidates");
    require(!seen.count("--boundary-diagnostics") || o.diagnostics,"Local construction requires sign diagnostics");
    o.diagnostics = true;
  }
  return o;
}
std::vector<VoxelPoint> readCloud(const fs::path& path)
{
  std::ifstream in(path); require(bool(in),"Cannot read XYZ"); std::string line; std::vector<VoxelPoint> points;
  while(std::getline(in,line))
  {
    const auto comment = line.find('#'); if(comment != std::string::npos) line.resize(comment);
    std::istringstream row(line); row >> std::ws; if(row.eof()) continue;
    double x,y,z;
    require(bool(row >> x >> y >> z) && std::isfinite(x) && std::isfinite(y) && std::isfinite(z),"Invalid XYZ point");
    row >> std::ws; require(row.eof(),"XYZ must have exactly three columns"); points.emplace_back(x,y,z);
  }
  require(!in.bad() && !points.empty(),"Empty/unreadable XYZ"); return points;
}
std::ofstream file(const fs::path& path)
{
  std::ofstream out(path); require(bool(out),"Cannot write boundary output");
  out.exceptions(std::ios::badbit | std::ios::failbit); out << std::setprecision(17) << std::boolalpha; return out;
}
template<class Range> void array(std::ostream& out,const Range& values)
{ out << '['; bool first = true; for(const auto& v : values) { if(!first) out << ','; first = false; out << v; } out << ']'; }
HS::Vec3 middle(const OctreeBox& box)
{ return {box.lower[0]+(box.upper[0]-box.lower[0])/2,box.lower[1]+(box.upper[1]-box.lower[1])/2,box.lower[2]+(box.upper[2]-box.lower[2])/2}; }

struct OccupancySource
{
  const Options& options;
  const std::vector<VoxelPoint>& raw;
  std::optional<OctreeOccupancy> octree;
  std::optional<SparseVoxelOccupancy> voxels;
  std::optional<OctreeCellPartition> partition;
  double occupancy_ms=0, partition_ms=0, lazy_fallback_ms=0;
  OccupancySource(const Options& o,const std::vector<VoxelPoint>& points): options(o),raw(points)
  {
    const auto start=Clock::now();
    if(o.occupancy=="voxels") voxels=buildSparseVoxelOccupancy(points,o.h);
    else octree=buildOctreeOccupancy(points,o.h,o.depth);
    occupancy_ms=ms(start);
    if(octree) { const auto t=Clock::now(); partition=partitionOctreeCells(*octree,pruneOctree(*octree),true); partition_ms=ms(t); }
  }
  const std::vector<VoxelIndex>& keys() const { return voxels ? voxels->occupied_voxels : octree->occupied_voxels; }
  const OctreeBox& samplingBox() const { return voxels ? voxels->sampling_box : octree->nodes[0].box; }
  bool contains(const HS::Vec3& p) const { return voxels ? voxels->contains(p) : octree->contains(p); }
  OctreeBox fineBox(std::size_t i) const
  { return voxels ? octreeGridBox(keys()[i],{1,1,1},options.h) : octree->nodes[octree->leaf_nodes[i]].box; }
  std::size_t storageBytes() const
  {
    std::size_t bytes=0;
    if(voxels)
    {
      bytes+=sizeof(SparseVoxelOccupancy)+voxels->occupied_voxels.capacity()*sizeof(VoxelIndex)+
        voxels->dense_bits.capacity()*sizeof(std::uint64_t);
      for(const auto& signs:voxels->boundary_plane_signs)
        bytes+=signs.capacity()*sizeof(unsigned char);
    }
    if(octree)
    {
      bytes+=sizeof(OctreeOccupancy)+octree->occupied_voxels.capacity()*sizeof(VoxelIndex)+
        octree->nodes.capacity()*sizeof(OctreeNode)+octree->leaf_nodes.capacity()*sizeof(int);
      for(const auto& n:octree->nodes) bytes+=n.point_indices.capacity()*sizeof(std::size_t);
    }
    return bytes; // owned occupancy payload only; not allocator overhead/RSS/peak
  }
  void ensureLegacyFallback()
  {
    if(octree) return;
    const auto start=Clock::now(); octree=buildOctreeOccupancy(raw,options.h,options.depth);
    partition=partitionOctreeCells(*octree,pruneOctree(*octree),true); lazy_fallback_ms=ms(start);
    require(octree->occupied_voxels==voxels->occupied_voxels,"Lazy fallback changed occupancy");
  }
};
std::vector<HS::Vec3> randomPoints(const OccupancySource& o,std::size_t count)
{
  std::mt19937_64 generator(20260909); std::uniform_real_distribution<double> unit(-.05,1.05);
  std::vector<HS::Vec3> points; points.reserve(count); const auto& b = o.samplingBox();
  for(std::size_t i = 0; i < count; ++i)
    points.push_back({b.lower[0]+unit(generator)*(b.upper[0]-b.lower[0]),b.lower[1]+unit(generator)*(b.upper[1]-b.lower[1]),
                      b.lower[2]+unit(generator)*(b.upper[2]-b.lower[2])});
  return points;
}
template<class F> double queryMedian(const std::vector<HS::Vec3>& points,F function)
{
  std::vector<double> times; volatile double sink = 0;
  for(const auto& p : points) sink = function(p); // warmup outside timed section
  for(int repeat = 0; repeat < 7; ++repeat)
  {
    double sum = 0; const auto start = Clock::now();
    for(int round = 0; round < 5; ++round) for(const auto& p : points) sum += function(p);
    times.push_back(ms(start)*1000/(5*points.size())); sink = sum;
  }
  (void)sink; std::sort(times.begin(),times.end()); return times[3];
}

void writeOutputs(const Options& o,const OccupancySource& occupancy,const OctreeBoundary& boundary,
                  const BoundaryTreeBatch& batch,const Metrics& m,const BoundarySignDiagnostics* flat_sign,
                  const BoundaryLocalResult* local)
{
  const auto& r = batch.candidates.at(batch.selected);
  const auto& tree = local ? local->tree : r.tree;
  const bool strict_sign = local && local->diagnostics.strict_sign_certified;
  fs::create_directories(o.output);
  auto candidates = file(o.output/"axis_candidates.json");
  candidates << "{\"requested_sweep\":\"" << o.requested_sweep << "\",\"selected_axis\":\"" << "xyz"[r.sweep_axis]
             << "\",\"selected_index\":" << batch.selected
             << ",\"selection_rule\":\"min(fallback,tree_nodes,plane_references,finite_prisms,axis); ties X,Y,Z\""
             << ",\"all_candidates_ms\":" << m.at("all_candidates_ms") << ",\"selection_ms\":" << batch.selection_ms << ",\"candidates\":[";
  for(std::size_t i = 0; i < batch.candidates.size(); ++i)
  {
    const auto& c = batch.candidates[i];
    const auto& ct = local ? local->tree : c.tree;
    candidates << (i ? "," : "") << "{\"axis\":\"" << "xyz"[c.sweep_axis] << "\",\"fallback\":" << c.fallback
      << ",\"expression_mode\":\"" << o.expression << "\",\"fallback_reason\":\"" << c.fallback_reason
      << "\",\"tree_nodes\":" << HS::serializeTree(ct.simplified).nodes.size()
      << ",\"plane_references\":" << ct.leaf_references << ",\"unique_halfspaces\":" << ct.leaves.size()
      << ",\"finite_prisms\":" << c.prisms.size() << ",\"prisms_before_merge\":" << c.prisms_before_merge
      << ",\"compressed_cells\":" << c.compressed_cells << ",\"event_updates\":" << c.event_updates
      << ",\"build_ms\":" << batch.candidate_build_ms[i]+(local ? local->core_ms : 0)
      << ",\"flat_build_ms\":" << batch.candidate_build_ms[i] << ",\"local_build_ms\":" << (local ? local->core_ms : 0)
      << ",\"registry_ms\":" << c.registry_ms
      << ",\"sweep_ms\":" << c.sweep_ms << ",\"partition_ms\":" << c.partition_ms << ",\"merge_ms\":" << c.merge_ms
      << ",\"tree_ms\":" << c.tree_ms+(local ? local->construction_ms+local->simplify_ms : 0)
      << ",\"strict_sign_certified\":" << strict_sign << ",\"certificate_passed\":true,\"probes_passed\":true}";
  }
  candidates << "]}\n"; candidates.close();
  const std::string prefix = r.fallback ? "obstacle_tree" : local ? "free_local_tree" : "free_closure_tree";
  const std::string pipeline = r.fallback ? "octree-rect-fallback" : local ? "octree-boundary-local" : "octree-boundary-closure";
  writeTreeJson(o.output/(prefix+"_raw.json"),tree,tree.raw,"raw",pipeline);
  writeTreeJson(o.output/(prefix+"_simplified.json"),tree,tree.simplified,"simplified",pipeline);
  if(local && !r.fallback)
  {
    writeTreeJson(o.output/"free_closure_tree_raw.json",r.tree,r.tree.raw,"raw","octree-boundary-closure");
    writeTreeJson(o.output/"free_closure_tree_simplified.json",r.tree,r.tree.simplified,"simplified","octree-boundary-closure");
  }
  if(flat_sign) writeBoundaryLocalReport(o.output/"boundary_sign_diagnostics.json",r,*flat_sign,local);
  auto manifest = file(o.output/"strict_free_space.json");
  manifest << "{\"schema\":1,\"length_unit\":\"m\",\"target\":\"R3_minus_closed_occupied_voxel_union\",\"tree_set\":\""
           << (r.fallback ? "closed_occupied_union" : "closure_of_nonoccupied_domain") << "\",\"classifier\":\""
           << (r.fallback ? "root_gt_zero" : strict_sign ? "root_lt_zero" : "root_lt_zero_or_root_eq_zero_and_not_occupied")
           << "\",\"tree\":\"" << prefix << "_simplified.json\",\"raw_tree\":\"" << prefix
           << "_raw.json\",\"occupancy_reference\":\"occupancy_reference.json\",\"requires_zero_guard\":" << (!r.fallback && !strict_sign)
           << ",\"expression_mode\":\"" << o.expression << "\",\"occupancy_backend\":\"" << o.occupancy
           << "\",\"occupancy_schema\":" << (occupancy.voxels ? 2 : 1) << ",\"strict_sign_certified\":" << strict_sign
           << ",\"fallback\":" << r.fallback << ",\"fallback_reason\":\"" << r.fallback_reason
           << "\",\"requested_sweep\":\"" << o.requested_sweep << "\",\"selected_axis\":\"" << "xyz"[r.sweep_axis]
           << "\",\"known_free_space_claim\":false,\"signed_distance_claim\":false}\n";
  manifest.close();
  auto reference = file(o.output/"occupancy_reference.json");
  reference << "{\"length_unit\":\"m\",\"closed_boxes\":true,\"voxel_size\":" << o.h;
  if(occupancy.voxels) reference << ",\"schema\":2,\"representation\":\"sorted_voxel_indices\",\"root_id\":null";
  else reference << ",\"root_id\":0";
  reference << ",\"occupied_voxel_indices\":[";
  for(std::size_t i = 0; i < occupancy.keys().size(); ++i) { if(i) reference << ','; array(reference,occupancy.keys()[i]); }
  reference << "],\"nodes\":[";
  for(std::size_t i = 0; !occupancy.voxels && i < occupancy.octree->nodes.size(); ++i)
  {
    const auto& n = occupancy.octree->nodes[i]; reference << (i ? "," : "") << "{\"id\":" << i << ",\"lower\":";
    array(reference,n.box.lower); reference << ",\"upper\":"; array(reference,n.box.upper);
    reference << ",\"leaf\":" << (n.width == 1) << ",\"children\":"; array(reference,n.children); reference << '}';
  }
  reference << "]}\n"; reference.close();
  auto supports = file(o.output/"boundary_supports.json");
  supports << "{\"length_unit\":\"m\",\"planes\":[";
  for(std::size_t i = 0; i < r.supports.size(); ++i)
  {
    const auto& s = r.supports[i]; supports << (i ? "," : "") << "{\"id\":" << i << ",\"axis\":" << s.axis
      << ",\"grid_coordinate\":" << s.coordinate << ",\"coordinate_m\":" << static_cast<double>(s.coordinate)*o.h
      << ",\"outward_sign_bits\":" << s.outward_signs << ",\"source_patch_ids\":";
    array(supports,s.patches); supports << '}';
  }
  supports << "],\"leaf_sources\":[";
  for(std::size_t i = 0; i < r.leaf_sources.size(); ++i)
    supports << (i ? "," : "") << "{\"leaf_id\":" << i << ",\"support_id\":" << r.leaf_sources[i].first
             << ",\"coefficient_sign\":" << r.leaf_sources[i].second << '}';
  supports << "],\"fallback_leaves_are_original_volume_planes\":" << r.fallback << "}\n"; supports.close();
  auto prisms = file(o.output/"free_prisms.json");
  prisms << "{\"length_unit\":\"m\",\"voxel_size\":" << o.h << ",\"geometry\":\"bounded_closed_free_prisms_plus_six_exterior_halfspaces\",\"fallback\":" << r.fallback << ",\"prisms\":[";
  for(std::size_t i = 0; i < r.prisms.size(); ++i)
  { prisms << (i ? "," : "") << "{\"lower_index\":"; array(prisms,r.prisms[i].lower); prisms << ",\"upper_index\":"; array(prisms,r.prisms[i].upper); prisms << '}'; }
  prisms << "]}\n"; prisms.close();
  if(occupancy.voxels) writeVoxelBoundary(o.output,o.h,boundary);
  else writeOctreeBoundary(o.output,*occupancy.octree,boundary); // OBSTACLE boundary, not free-space volume.
  auto validation = file(o.output/"validation_report.json");
  validation << "{\"passed\":true,\"integer_coverage_certificate\":true,\"strict_classifier_checked\":true,\"fallback\":" << r.fallback
             << ",\"strict_sign_certified\":" << strict_sign << ",\"local_target_complete\":" << strict_sign
             << ",\"sweep_certificate_applied\":" << !r.fallback << ",\"sampled_false_free\":" << m.at("false_free")
             << ",\"sampled_missed_free\":" << m.at("missed_free") << ",\"scalar_mismatches\":" << m.at("scalar_mismatches")
             << ",\"raw_points_classified_free\":" << m.at("raw_points_classified_free") << ",\"tested_points\":" << m.at("tested_points")
             << ",\"validated_candidates\":" << batch.candidates.size() << ",\"candidate_point_checks\":" << m.at("candidate_point_checks")
             << ",\"zero_guard_calls\":" << m.at("zero_guard_calls") << ",\"standalone_tree_is_strict_complement\":false,\"failures\":[]}\n";
  validation.close();
}
} // namespace

int runOctreeBoundaryTreePipeline(int argc,char** argv)
{
  for(int i=3;i+1<argc;i+=2)
    if(std::string(argv[i])=="--boundary-expression" && std::string(argv[i+1])=="direct")
      return runBoundaryDirectPipeline(argc,argv);
  try
  {
    const auto o = parse(argc,argv); Metrics m;
    auto start = Clock::now(); const auto raw = readCloud(o.input); m["input_io_ms"] = ms(start);
    const auto core = Clock::now(); start = Clock::now();
    OccupancySource occupancy(o,raw);
    m["occupancy_ms"]=occupancy.occupancy_ms; m["baseline_partition_ms"]=occupancy.partition_ms;
    const auto boundary=occupancy.voxels ? extractVoxelBoundary(occupancy.keys()) :
      extractOctreeBoundary(*occupancy.octree,*occupancy.partition);
    const auto batch=occupancy.voxels ? buildVoxelBoundaryTreeBatch(o.h,boundary,o.work,o.best_axis,[&] {
      occupancy.ensureLegacyFallback();
      std::vector<HS::Vec3> vertices; auto clusters=octreeAabbClusters(*occupancy.partition,vertices);
      auto t=buildClusterTree(clusters,0,"octree_aabb"); t.simplified=HS::simplifyToFixedPoint(t.raw); return t;
    }) : buildBoundaryTreeBatch(*occupancy.octree,*occupancy.partition,boundary,o.work,o.best_axis);
    const auto& r = batch.candidates.at(batch.selected);
    std::optional<BoundaryLocalResult> local;
    if(o.expression == "local") local = buildBoundaryLocalTree(r,o.local);
    const auto& tree = local ? local->tree : r.tree;
    const auto strict_query = [&](const HS::Vec3& p,bool* zero = nullptr) {
      return local ? boundaryLocalStrictFree(*local,*occupancy.octree,p,zero) : occupancy.voxels ?
        boundaryStrictFree(r,*occupancy.voxels,p,zero) : boundaryStrictFree(r,*occupancy.octree,p,zero);
    };
    m["core_runtime_ms"] = ms(core);
    m["lazy_fallback_ms"]=occupancy.lazy_fallback_ms;
    m["octree_nodes"]=occupancy.octree ? occupancy.octree->nodes.size() : 0;
    m["octree_built"]=occupancy.octree.has_value();
    m["occupancy_storage_bytes"]=occupancy.storageBytes();
    m["boundary_extraction_ms"] = boundary.extraction_ms; m["boundary_merge_ms"] = boundary.merge_ms;
    // Phase times SUM all attempted candidates, including work discarded by selection/fallback.
    for(const auto& c : batch.candidates)
    {
      m["registry_compression_ms"] += c.registry_ms; m["sweep_ms"] += c.sweep_ms; m["partition_2d_ms"] += c.partition_ms;
      m["prism_merge_ms"] += c.merge_ms; m["tree_build_simplify_ms"] += c.tree_ms;
    }
    m["all_candidates_ms"] = batch.all_candidates_ms+(local ? local->core_ms : 0); m["selection_ms"] = batch.selection_ms;
    m["candidate_count"] = batch.candidates.size();
    m["selected_candidate_build_ms"] = batch.candidate_build_ms[batch.selected]+(local ? local->core_ms : 0);
    if(local) m["tree_build_simplify_ms"] += local->construction_ms+local->simplify_ms;
    start = Clock::now();
    if(occupancy.voxels) validateSparseVoxelOccupancy(*occupancy.voxels,raw);
    if(occupancy.octree) validateOctreeOccupancy(*occupancy.octree,raw);
    std::optional<OctreeBoundary> fallback_boundary;
    for(const auto& c : batch.candidates)
    {
      if(occupancy.voxels && !c.fallback) validateVoxelBoundaryTree(*occupancy.voxels,boundary,c);
      else
      {
        // Fallback uses rect provenance; voxel boundary owner IDs intentionally differ.
        if(occupancy.voxels && !fallback_boundary)
          fallback_boundary=extractOctreeBoundary(*occupancy.octree,*occupancy.partition);
        validateBoundaryTree(*occupancy.octree,*occupancy.partition,
                             fallback_boundary ? *fallback_boundary : boundary,c);
      }
    }
    m["certificate_ms"] = ms(start);
    BoundarySignDiagnostics flat_sign;
    m["sign_diagnostics_ms"] = 0;
    if(o.diagnostics)
    {
      flat_sign = occupancy.voxels ? diagnoseBoundarySigns(*occupancy.voxels,r,r.tree,o.local) :
        diagnoseBoundarySigns(*occupancy.octree,r,r.tree,o.local); m["sign_diagnostics_ms"] += flat_sign.time_ms;
      if(local)
      {
        local->diagnostics = diagnoseBoundarySigns(*occupancy.octree,r,local->tree,o.local); m["sign_diagnostics_ms"] += local->diagnostics.time_ms;
        if(!r.fallback)
          require(local->diagnostics.applied && local->diagnostics.closure_set_certified && !local->diagnostics.unsafe_free_count &&
                  !local->diagnostics.boundary_sign_errors,"Local candidate lacks a valid geometry/sign-safety certificate (possibly work cap)");
      }
      const auto& d = local ? local->diagnostics : flat_sign;
      m["sign_certificate_applied"] = d.applied; m["strict_sign_certified"] = d.strict_sign_certified;
      m["artificial_zero_faces"] = d.free_zeros[2]; m["artificial_zero_edges"] = d.free_zeros[1]; m["artificial_zero_vertices"] = d.free_zeros[0];
      m["unsafe_free_count"] = d.unsafe_free_count; m["missed_free_count"] = d.missed_free_count;
      m["true_boundary_sign_errors"] = d.boundary_sign_errors;
      m["occupied_interior_zero_count"] = std::accumulate(d.occupied_zeros.begin(),d.occupied_zeros.end(),std::size_t(0));
      m["flat_artificial_zero_faces"] = flat_sign.free_zeros[2]; m["flat_artificial_zero_edges"] = flat_sign.free_zeros[1];
      m["flat_artificial_zero_vertices"] = flat_sign.free_zeros[0];
    }
    const auto sign_strata=flat_sign.strata_count+(local ? local->diagnostics.strata_count : 0);
    const auto sign_batches=[&](std::size_t count) { return o.local.packed_sign_propagation ? (count+63)/64 : count; };
    m["sign_packed_propagation"]=o.local.packed_sign_propagation;
    m["sign_strata_count"]=sign_strata;
    m["sign_propagation_batches"]=sign_batches(flat_sign.strata_count)+
      (local ? sign_batches(local->diagnostics.strata_count) : 0);
    start = Clock::now(); std::vector<HS::Vec3> ignored;
    const bool fine_reference=o.voxel_validation_reference || occupancy.voxels.has_value();
    std::vector<ConvexCluster> base_clusters;
    if(fine_reference)
    {
      std::vector<OctreeBox> boxes; boxes.reserve(occupancy.keys().size());
      for(std::size_t i=0;i<occupancy.keys().size();++i) boxes.push_back(occupancy.fineBox(i));
      std::vector<std::size_t> ids(boxes.size()); std::iota(ids.begin(),ids.end(),0);
      base_clusters=aabbClusters(boxes,ids,ignored);
    }
    else base_clusters=octreeAabbClusters(*occupancy.partition,ignored);
    auto base_tree = buildClusterTree(base_clusters,0,"octree_aabb"); base_tree.simplified = HS::simplifyToFixedPoint(base_tree.raw);
    m["validation_baseline_tree_ms"] = ms(start);
    start = Clock::now();
    m["false_free"] = m["missed_free"] = m["raw_points_classified_free"] = m["scalar_mismatches"] = m["tested_points"] = m["zero_guard_calls"] = 0;
    m["candidate_point_checks"] = m["all_candidate_zero_guard_calls"] = 0;
    const auto check = [&](const HS::Vec3& p, bool is_raw) {
      const bool expected = !occupancy.contains(p); m["tested_points"]++;
      for(std::size_t i = 0; i < batch.candidates.size(); ++i)
      {
        const auto& c = batch.candidates[i];
        const auto& ct = local ? local->tree : c.tree;
        bool zero = false; const bool free = local ? strict_query(p,&zero) : occupancy.voxels ? boundaryStrictFree(c,*occupancy.voxels,p,&zero) : boundaryStrictFree(c,*occupancy.octree,p,&zero);
        m["candidate_point_checks"]++; m["all_candidate_zero_guard_calls"] += zero;
        if(i == batch.selected) m["zero_guard_calls"] += zero;
        m["false_free"] += free && !expected; m["missed_free"] += !free && expected;
        m["raw_points_classified_free"] += is_raw && free;
        m["scalar_mismatches"] += HS::evaluate(ct.raw,ct.leaves,p) != HS::evaluate(ct.simplified,ct.leaves,p);
      }
      require((HS::evaluate(base_tree.simplified,base_tree.leaves,p) > 0) == expected,"Baseline/oracle disagreement");
    };
    for(const auto& p : raw) check({p.x(),p.y(),p.z()},true);
    for(std::size_t id=0;id<occupancy.keys().size();++id)
    {
      const auto box = occupancy.fineBox(id); check(middle(box),false);
      for(unsigned mask = 0; mask < 8; ++mask)
        check({(mask&1) ? box.upper[0] : box.lower[0],(mask&2) ? box.upper[1] : box.lower[1],(mask&4) ? box.upper[2] : box.lower[2]},false);
    }
    for(const auto& patch : boundary.patches)
    {
      const auto corners = octreeBoundaryCorners(patch); double p[3]{};
      for(unsigned k = 0; k < 3; ++k)
      { const double a = static_cast<double>(corners[0][k])*o.h,b = static_cast<double>(corners[2][k])*o.h; p[k] = a+(b-a)/2; }
      check({p[0],p[1],p[2]},false);
      p[patch.axis] -= o.h*.25; check({p[0],p[1],p[2]},false);
      p[patch.axis] += o.h*.5; check({p[0],p[1],p[2]},false);
    }
    if(!r.fallback)
      for(std::size_t x = 0; x+1 < r.coordinates[0].size(); ++x)
        for(std::size_t y = 0; y+1 < r.coordinates[1].size(); ++y)
          for(std::size_t z = 0; z+1 < r.coordinates[2].size(); ++z)
          {
            const std::size_t ix[3]{x,y,z}; double p[3]{};
            for(unsigned k = 0; k < 3; ++k)
            { const double a = static_cast<double>(r.coordinates[k][ix[k]])*o.h,b = static_cast<double>(r.coordinates[k][ix[k]+1])*o.h; p[k] = a+(b-a)/2; }
            check({p[0],p[1],p[2]},false);
          }
    for(const auto& p : randomPoints(occupancy,o.samples)) check(p,false);
    m["probe_validation_ms"] = ms(start);
    require(m["false_free"] == 0 && m["missed_free"] == 0 && m["raw_points_classified_free"] == 0 && m["scalar_mismatches"] == 0,
            "Strict free/obstacle complement validation failed");
    m["validation_ms"] = m["certificate_ms"]+m["sign_diagnostics_ms"]+m["validation_baseline_tree_ms"]+m["probe_validation_ms"];
    m["total_ms"] = ms(core); // construction + validation, excludes input/output I/O and optional query benchmark
    m["raw_points"] = raw.size(); m["occupied_voxels"] = occupancy.keys().size(); m["voxel_size"] = o.h;
    m["raw_exposed_faces"] = boundary.faces.size(); m["merged_boundary_patches"] = boundary.patches.size();
    m["unique_boundary_support_planes"] = r.supports.size(); m["unique_oriented_boundary_halfspaces"] = 0;
    for(const auto& s : r.supports) m["unique_oriented_boundary_halfspaces"] += bool(s.outward_signs&1)+bool(s.outward_signs&2);
    m["x_slabs"] = r.coordinates[0].size()-1; m["compressed_y_cells"] = r.coordinates[1].size()-1; m["compressed_z_cells"] = r.coordinates[2].size()-1;
    m["compressed_x_cells"] = r.coordinates[0].size()-1; // x_slabs retained as world-X legacy alias
    m["sweep_slabs"] = r.coordinates[r.sweep_axis].size()-1;
    m["compressed_u_cells"] = r.coordinates[(r.sweep_axis+1)%3].size()-1;
    m["compressed_v_cells"] = r.coordinates[(r.sweep_axis+2)%3].size()-1;
    m["compressed_cells"] = r.compressed_cells; m["event_updates"] = r.event_updates;
    m["free_prisms_before_merge"] = r.prisms_before_merge; m["free_prisms_after_merge"] = r.prisms.size();
    m["merge_passes"] = r.merge_passes; m["merge_fixed_point"] = r.merge_fixed_point;
    m["max_compressed_cells"] = o.work.max_compressed_cells; m["max_event_updates"] = o.work.max_event_updates;
    m["max_free_prisms"] = o.work.max_free_prisms; m["max_merge_passes"] = o.work.max_merge_passes;
    m["tree_nodes"] = HS::serializeTree(tree.simplified).nodes.size(); m["raw_tree_nodes"] = HS::serializeTree(tree.raw).nodes.size();
    m["plane_references"] = tree.leaf_references; m["unique_halfspaces"] = tree.leaves.size();
    if(local)
    {
      m["adjacency_face_count"] = local->adjacency.size(); m["adjacency_checks"] = local->adjacency_checks;
      m["cancelled_interface_pairs"] = local->cancelled_interface_pairs; m["logical_merge_count"] = local->logical_merges;
      m["factoring_failed_count"] = local->factoring_failed; m["local_attempts"] = local->attempts;
      m["local_build_ms"] = local->core_ms; m["adjacency_ms"] = local->adjacency_ms;
      m["local_construction_ms"] = local->construction_ms; m["local_simplify_ms"] = local->simplify_ms;
    }
    m["constructed_tree_nodes"] = m["raw_tree_nodes"]; m["simplified_tree_nodes"] = m["tree_nodes"];
    m["expanded_node_references"] = HS::treeNodeReferenceCount(tree.simplified); m["core_build_time"] = m["core_runtime_ms"];
    m["baseline_tree_nodes"] = HS::serializeTree(base_tree.simplified).nodes.size(); m["baseline_regions"] = base_clusters.size();
    m["fallback"] = r.fallback; m["boundary_leaf_provenance_applicable"] = !r.fallback;
    m["artificial_boundary_leaf_count"] = 0; m["passed"] = 1;
    m["query_benchmark_enabled"] = o.query_benchmark;
    if(o.query_benchmark)
    {
      start = Clock::now(); const auto queries = randomPoints(occupancy,1024);
      m["raw_tree_query_us"] = queryMedian(queries,[&](const auto& p) { return HS::evaluate(tree.raw,tree.leaves,p); });
      m["simplified_tree_query_us"] = queryMedian(queries,[&](const auto& p) { return HS::evaluate(tree.simplified,tree.leaves,p); });
      m["strict_free_query_us"] = queryMedian(queries,[&](const auto& p) { return strict_query(p); });
      m["query_time"] = m["strict_free_query_us"];
      m["baseline_strict_free_query_us"] = queryMedian(queries,[&](const auto& p) { return HS::evaluate(base_tree.simplified,base_tree.leaves,p)>0; });
      std::vector<HS::Vec3> seams;
      for(const auto& patch : boundary.patches)
      {
        const auto corners = octreeBoundaryCorners(patch);
        for(const auto& key : corners) seams.push_back({static_cast<double>(key[0])*o.h,static_cast<double>(key[1])*o.h,static_cast<double>(key[2])*o.h});
      }
      if(seams.size() > 1024) seams.resize(1024);
      m["boundary_strict_free_query_us"] = queryMedian(seams,[&](const auto& p) { return strict_query(p); });
      m["occupancy_query_us"] = queryMedian(queries,[&](const auto& p) { return occupancy.contains(p); });
      m["boundary_occupancy_query_us"] = queryMedian(seams,[&](const auto& p) { return occupancy.contains(p); });
      std::vector<HS::Vec3> artificial;
      for(const auto& prism:r.prisms)
      {
        const auto box=octreeGridBox(prism.lower,{prism.upper[0]-prism.lower[0],prism.upper[1]-prism.lower[1],
                                                prism.upper[2]-prism.lower[2]},o.h);
        for(unsigned axis=0;axis<3;++axis) for(unsigned high=0;high<2;++high)
        {
          auto p=middle(box); const auto c=high ? box.upper[axis] : box.lower[axis];
          if(axis==0) p.x=c; else if(axis==1) p.y=c; else p.z=c;
          if(artificial.size()<1024 && HS::evaluate(tree.simplified,tree.leaves,p)==0 && !occupancy.contains(p)) artificial.push_back(p);
        }
      }
      m["artificial_seam_query_count"]=artificial.size();
      m["artificial_seam_strict_query_us"]=artificial.empty() ? 0 : queryMedian(artificial,[&](const auto& p) { return strict_query(p); });
      std::size_t zeros = 0;
      for(const auto& p : queries) { bool used; strict_query(p,&used); zeros += used; }
      m["query_zero_guard_fraction"] = static_cast<double>(zeros)/queries.size();
      m["query_benchmark_ms"] = ms(start);
      // Export exactly the timed query corpora so before/after comparison can
      // check them byte-for-byte, including guard-heavy artificial seams.
      const auto io=Clock::now(); fs::create_directories(o.output);
      auto samples=file(o.output/"query_samples.json"); samples << '{';
      const auto group=[&](const char* name,const std::vector<HS::Vec3>& ps,bool comma) {
        samples << (comma ? "," : "") << '"' << name << "\":[";
        for(std::size_t i=0;i<ps.size();++i) samples << (i ? "," : "") << '[' << ps[i].x << ',' << ps[i].y << ',' << ps[i].z << ']';
        samples << ']';
      };
      group("random",queries,false); group("boundary",seams,true); group("artificial_seams",artificial,true);
      samples << "}\n"; samples.close(); m["query_samples_io_ms"]=ms(io);
    }
    start = Clock::now(); writeOutputs(o,occupancy,boundary,batch,m,o.diagnostics ? &flat_sign : nullptr,local ? &*local : nullptr); m["output_io_ms"] = ms(start);
    auto metrics = file(o.output/"benchmark.json");
    metrics << "[{\"pipeline\":\"octree-boundary\",\"sweep_axis\":\"" << o.requested_sweep << "\",\"selected_axis\":\"" << "xyz"[r.sweep_axis]
            << "\",\"expression_mode\":\"" << o.expression << "\",\"occupancy_backend\":\"" << o.occupancy
            << "\",\"validation_reference\":\"" << (fine_reference ? "fine_voxels" : "rect")
            << "\",\"length_unit\":\"m\",\"time_unit\":\"ms\",\"query_time_unit\":\"us\",\"fallback_reason\":\"" << r.fallback_reason << '"';
    for(const auto& item : m) metrics << ",\n\"" << item.first << "\":" << item.second;
    metrics << "}]\n"; metrics.close();
    std::cout << "Boundary sweep " << o.requested_sweep << "; selected " << "xyz"[r.sweep_axis] << "; strict query uses "
              << (r.fallback ? "original obstacle complement" : local && local->diagnostics.strict_sign_certified ? "certified root<0" : "closure tree + zero occupancy guard") << '\n';
    if(local && !r.fallback && !local->diagnostics.strict_sign_certified)
      std::cout << "Local target INCOMPLETE: residual artificial zeros; occupancy guard retained. See boundary_sign_diagnostics.json.\n";
    for(const auto& item : m) std::cout << item.first << ": " << item.second << '\n';
    return 0;
  }
  catch(const std::exception& e)
  { std::cerr << "Boundary tree error: " << e.what() << '\n'; return 1; }
}
} // namespace rokae_demo
