#include <rokae_demo/octree_convex_approximation.hpp>

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
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

struct Options
{
  fs::path input,output;
  double voxel_size = 0;
  unsigned max_depth = 20;
  std::size_t samples = 1000;
  OctreeConvexOptions convex;
};

std::size_t integer(const std::string& value)
{
  if(value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
    throw std::invalid_argument("Expected a nonnegative integer");
  const auto number = std::stoull(value);
  if(number > std::numeric_limits<std::size_t>::max()) throw std::out_of_range("Integer too large");
  return static_cast<std::size_t>(number);
}

double nonnegative(const std::string& value)
{
  std::size_t consumed = 0;
  const double number = std::stod(value,&consumed);
  if(consumed != value.size() || !std::isfinite(number) || number < 0) throw std::invalid_argument("Expected a finite nonnegative parameter");
  return number;
}

Options parse(int argc, char** argv)
{
  if(argc < 3) throw std::invalid_argument("Expected input.xyz and output directory");
  Options result;
  result.input = argv[1]; result.output = argv[2];
  std::set<std::string> seen;
  for(int i = 3; i < argc; i += 2)
  {
    const std::string option = argv[i];
    if(i+1 >= argc || !seen.insert(option).second) throw std::invalid_argument("Missing/repeated option: " + option);
    const std::string value = argv[i+1];
    if(option == "--pipeline") { if(value != "octree") throw std::invalid_argument("Expected --pipeline octree"); }
    else if(option == "--octree-approx") { if(value != "hull") throw std::invalid_argument("Expected --octree-approx hull"); }
    else if(option == "--octree-pruning") { if(value != "rect") throw std::invalid_argument("Hull mode requires the exact rect baseline"); }
    else if(option == "--octree-boundary")
    {
      if(value != "none") throw std::invalid_argument("Occupied-voxel boundary is not the approximated geometry; use convex_cells.off in hull mode");
    }
    else if(option == "--voxel-size") result.voxel_size = nonnegative(value);
    else if(option == "--max-volume-inflation") result.convex.max_volume_inflation = nonnegative(value);
    else if(option == "--max-fill-distance") result.convex.max_fill_distance = nonnegative(value);
    else if(option == "--convex-max-voxels") result.convex.max_candidate_voxels = integer(value);
    else if(option == "--convex-test-max-boxes") result.convex.max_certificate_boxes = integer(value);
    else if(option == "--validation-samples") result.samples = integer(value);
    else if(option == "--max-depth")
    {
      const auto depth = integer(value);
      if(depth > 52) throw std::invalid_argument("max-depth must be in [0,52]");
      result.max_depth = static_cast<unsigned>(depth);
    }
    else throw std::invalid_argument("Unsupported hull option: " + option);
  }
  if(!seen.count("--max-volume-inflation") || !seen.count("--max-fill-distance"))
    throw std::invalid_argument("Hull approximation requires explicit --max-volume-inflation and --max-fill-distance budgets (zero is allowed)");
  if(result.voxel_size <= 0 || !result.convex.max_candidate_voxels || !result.convex.max_certificate_boxes)
    throw std::invalid_argument("voxel-size and convex work limits must be positive");
  return result;
}

std::vector<VoxelPoint> readCloud(const fs::path& path)
{
  std::ifstream in(path);
  if(!in) throw std::runtime_error("Cannot read XYZ input");
  std::vector<VoxelPoint> result;
  std::string line;
  while(std::getline(in,line))
  {
    const auto comment = line.find('#');
    if(comment != std::string::npos) line.resize(comment);
    std::istringstream row(line); row >> std::ws;
    if(row.eof()) continue;
    double x,y,z;
    if(!(row >> x >> y >> z) || !std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) throw std::runtime_error("Invalid XYZ coordinates");
    row >> std::ws;
    if(!row.eof()) throw std::runtime_error("Expected exactly three XYZ columns");
    result.emplace_back(x,y,z);
  }
  if(in.bad() || result.empty()) throw std::runtime_error("Empty or unreadable XYZ input");
  return result;
}

std::ofstream file(const fs::path& path)
{
  std::ofstream out(path);
  if(!out) throw std::runtime_error("Cannot write " + path.string());
  out.exceptions(std::ios::badbit | std::ios::failbit);
  out << std::setprecision(17);
  return out;
}

template<class Range> void array(std::ostream& out, const Range& items)
{
  out << '['; bool first = true;
  for(const auto& item : items) { if(!first) out << ','; first = false; out << item; }
  out << ']';
}

struct Validation : ClusterValidation
{
  bool certificate = false;
  std::size_t tested = 0, raw_outside = 0, source_outside = 0, added = 0;
  std::size_t scalar_mismatch = 0, decision_mismatch = 0, cluster_mismatch = 0;
};

Validation validate(const Options& options, const OctreeOccupancy& octree, const OctreeCellPartition& baseline,
                    OctreeConvexResult& convex, const TreeRepresentation& tree, const std::vector<VoxelPoint>& raw)
{
  Validation result;
  validateOctreeOccupancy(octree,raw);
  validateOctreeConvex(octree,baseline,options.convex,convex);
  result.certificate = true;
  validateClusters(convex.clusters,convex.source_vertices,0,result);
  // A near-normal global leaf reuse would invalidate the per-cluster proof.
  // Verify exact exported coefficients in addition to requesting strict reuse.
  for(const auto& c : convex.clusters)
  {
    if(c.leaf_ids.size() != c.reduced_planes.size()) throw std::runtime_error("Missing tree plane references");
    for(std::size_t i = 0; i < c.leaf_ids.size(); ++i)
    {
      const auto& p = c.reduced_planes[i]; const auto& leaf = tree.leaves.at(c.leaf_ids[i]);
      if(p.normal.x != leaf.normal.x || p.normal.y != leaf.normal.y || p.normal.z != leaf.normal.z || p.offset != leaf.offset)
        throw std::runtime_error("Tree leaf reuse changed a certified plane");
    }
  }
  const auto check = [&](const HS::Vec3& p, bool raw_point, bool known_source) {
    const auto a = HS::evaluate(tree.raw,tree.leaves,p), b = HS::evaluate(tree.simplified,tree.leaves,p);
    if(!std::isfinite(a) || !std::isfinite(b)) throw std::runtime_error("Nonfinite tree evaluation");
    const bool old_inside = octree.contains(p), inside = b <= 0;
    ++result.tested;
    result.raw_outside += raw_point && !inside;
    result.source_outside += (old_inside || known_source) && !inside;
    result.added += !old_inside && inside;
    result.scalar_mismatch += a != b;
    result.decision_mismatch += (a <= 0) != inside;
    result.cluster_mismatch += insideClusterUnion(convex.clusters,p,0) != inside;
  };
  for(const auto& p : raw) check({p.x(),p.y(),p.z()},true,true);
  for(const auto& p : convex.source_vertices) check(p,false,true);
  for(const auto id : octree.leaf_nodes)
  {
    const auto& box = octree.nodes[id].box;
    check({box.lower[0]+(box.upper[0]-box.lower[0])/2,box.lower[1]+(box.upper[1]-box.lower[1])/2,
           box.lower[2]+(box.upper[2]-box.lower[2])/2},false,true);
  }
  std::mt19937_64 generator(20260908);
  std::uniform_real_distribution<double> fraction(-.05,1.05);
  const auto& domain = octree.nodes.front().box;
  for(std::size_t i = 0; i < options.samples; ++i)
  {
    HS::Vec3 p{domain.lower[0]+fraction(generator)*(domain.upper[0]-domain.lower[0]),
               domain.lower[1]+fraction(generator)*(domain.upper[1]-domain.lower[1]),
               domain.lower[2]+fraction(generator)*(domain.upper[2]-domain.lower[2])};
    check(p,false,false);
  }
  if(result.raw_outside || result.source_outside) result.failures.push_back("Raw/occupied source containment failed");
  if(result.scalar_mismatch || result.decision_mismatch || result.cluster_mismatch) result.failures.push_back("Cluster/tree/simplification inconsistency");
  return result;
}

void writeOutputs(const Options& options, const OctreeOccupancy& octree, const OctreeConvexResult& convex,
                  const TreeRepresentation& tree, const Validation& validation)
{
  fs::create_directories(options.output);
  auto reference = file(options.output / "octree.json");
  reference << "{\"pipeline\":\"octree-hull\",\"representation\":\"fine_occupancy_reference\",\"length_unit\":\"m\",\"voxel_size\":"
            << octree.voxel_size << ",\"occupied_voxel_indices\":[";
  for(std::size_t i = 0; i < octree.occupied_voxels.size(); ++i) { if(i) reference << ','; array(reference,octree.occupied_voxels[i]); }
  reference << "],\"voxel_cluster_ids\":"; array(reference,convex.voxel_cluster_ids);
  reference << ",\"root_id\":0,\"nodes\":[";
  for(std::size_t i = 0; i < octree.nodes.size(); ++i)
  {
    const auto& n = octree.nodes[i];
    reference << (i ? "," : "") << "{\"id\":" << i << ",\"lower_index\":"; array(reference,n.lower);
    reference << ",\"width_voxels\":" << n.width << ",\"depth\":" << n.depth << ",\"children\":"; array(reference,n.children);
    reference << ",\"occupied_count\":" << n.occupied_count << ",\"point_indices\":"; array(reference,n.point_indices);
    reference << '}';
  }
  reference << "]}\n"; reference.close();
  auto clusters = file(options.output / "convex_clusters.json");
  clusters << "{\"pipeline\":\"octree-hull\",\"length_unit\":\"m\",\"geometry_definition\":\"intersection_of_exported_halfspaces\",\"clusters\":[\n";
  std::size_t vertex_count = 0, triangle_count = 0;
  for(std::size_t i = 0; i < convex.clusters.size(); ++i)
  {
    const auto& c = convex.clusters[i]; const auto& cell = convex.cells[i];
    vertex_count += cell.mesh_vertices.size(); triangle_count += cell.mesh_triangles.size();
    clusters << (i ? ",\n" : "") << "{\"id\":" << c.id << ",\"type\":\"" << (cell.non_aabb ? "convex_polyhedron" : "aabb")
             << "\",\"source_node\":" << cell.source_node << ",\"baseline_cell\":" << cell.baseline_cell
             << ",\"source_voxel_ids\":"; array(clusters,c.cells);
    clusters << ",\"volume_m3\":" << cell.volume << ",\"fill_distance_bound_m\":" << cell.fill_distance_bound
             << ",\"leaf_ids\":"; array(clusters,c.leaf_ids);
    clusters << ",\"planes\":[";
    for(std::size_t j = 0; j < c.reduced_planes.size(); ++j)
    {
      const auto& p = c.reduced_planes[j];
      clusters << (j ? "," : "") << "{\"normal\":[" << p.normal.x << ',' << p.normal.y << ',' << p.normal.z << "],\"offset\":" << p.offset << '}';
    }
    clusters << "]}";
  }
  clusters << "\n]}\n"; clusters.close();
  auto off = file(options.output / "convex_cells.off");
  off << "OFF\n" << vertex_count << ' ' << triangle_count << " 0\n";
  for(const auto& cell : convex.cells) for(const auto& p : cell.mesh_vertices) off << p.x << ' ' << p.y << ' ' << p.z << '\n';
  std::size_t shift = 0;
  for(const auto& cell : convex.cells)
  {
    for(const auto& f : cell.mesh_triangles) off << "3 " << shift+f[0] << ' ' << shift+f[1] << ' ' << shift+f[2] << '\n';
    shift += cell.mesh_vertices.size();
  }
  off.close();
  writeTreeJson(options.output / "logic_tree_raw.json",tree,tree.raw,"raw","octree-hull");
  writeTreeJson(options.output / "logic_tree_simplified.json",tree,tree.simplified,"simplified","octree-hull");
  for(const auto& item : std::vector<std::pair<std::string,HS::ExpressionPtr>>{{"raw",tree.raw},{"simplified",tree.simplified}})
  {
    auto expression = file(options.output / ("logic_expression_"+item.first+".txt"));
    expression << "inside(x) =\n" << HS::booleanExpression(item.second,0) << "\n\nphi(x) =\n" << HS::scalarExpression(item.second,0) << '\n';
    expression.close();
  }
  auto report = file(options.output / "validation_report.json");
  report << std::boolalpha << "{\"pipeline\":\"octree-hull\",\"passed\":" << validation.failures.empty()
         << ",\"source_coverage_and_error_certificate\":" << validation.certificate
         << ",\"reference\":\"occupied_voxel_union\",\"classification_epsilon\":0,\"raw_points_outside_tree\":" << validation.raw_outside
         << ",\"source_inside_tree_outside\":" << validation.source_outside << ",\"tested_points\":" << validation.tested
         << ",\"added_occupied_samples\":" << validation.added << ",\"raw_simplified_mismatch\":" << validation.decision_mismatch
         << ",\"scalar_mismatch\":" << validation.scalar_mismatch << ",\"cluster_mismatch\":" << validation.cluster_mismatch
         << ",\"old_new_scalar_equality_required\":false,\"known_free_space_checked\":false,\"failures\":[";
  for(std::size_t i = 0; i < validation.failures.size(); ++i) report << (i ? "," : "") << '"' << validation.failures[i] << '"';
  report << "]}\n"; report.close();
}
}  // namespace

int runOctreeConvexPipeline(int argc, char** argv)
{
  try
  {
    const auto options = parse(argc,argv);
    Metrics metrics;
    auto start = Clock::now();
    const auto raw = readCloud(options.input);
    metrics["input_io_ms"] = ms(start);
    const auto total_start = Clock::now();
    start = Clock::now();
    const auto octree = buildOctreeOccupancy(raw,options.voxel_size,options.max_depth);
    metrics["split_time"] = ms(start);
    start = Clock::now();
    const auto baseline = partitionOctreeCells(octree,pruneOctree(octree),true);
    metrics["pruning_time"] = ms(start);
    start = Clock::now();
    auto convex = approximateOctreeConvex(octree,baseline,options.convex);
    metrics["convex_approximation_ms"] = ms(start);
    start = Clock::now();
    auto tree = buildClusterTree(convex.clusters,0,"octree_hull",true);
    tree.simplified = HS::simplifyToFixedPoint(tree.raw);
    metrics["tree_nodes"] = HS::serializeTree(tree.simplified).nodes.size();
    metrics["tree_time"] = ms(start);
    metrics["core_runtime_ms"] = ms(total_start);
    start = Clock::now();
    const auto validation = validate(options,octree,baseline,convex,tree,raw);
    metrics["validation_ms"] = ms(start);
    metrics["total_ms"] = ms(total_start);
    metrics["raw_points"] = raw.size();
    metrics["occupied_voxels"] = octree.occupied_voxels.size();
    metrics["octree_nodes"] = octree.nodes.size();
    metrics["voxel_resolution"] = options.voxel_size;
    metrics["baseline_convex_cells"] = baseline.cells.size();
    metrics["baseline_plane_references"] = 6*baseline.cells.size();
    metrics["baseline_tree_nodes"] = 7*baseline.cells.size()+(baseline.cells.size()>1);
    metrics["convex_cells"] = convex.clusters.size();
    metrics["non_aabb_cells"] = std::count_if(convex.cells.begin(),convex.cells.end(),[](const auto& c) { return c.non_aabb; });
    metrics["planes_before_reduce"] = metrics["planes_after_reduce"] = 0;
    for(const auto& c : convex.clusters)
    { metrics["planes_before_reduce"] += c.exact_planes.size(); metrics["planes_after_reduce"] += c.reduced_planes.size(); }
    metrics["unique_planes"] = tree.leaves.size();
    metrics["tree_leaf_references"] = tree.leaf_references;
    metrics["candidate_nodes"] = convex.candidate_nodes;
    metrics["candidate_hulls"] = convex.candidate_hulls;
    metrics["accepted_merges"] = convex.accepted_merges;
    metrics["rejected_volume"] = convex.rejected_volume;
    metrics["rejected_distance"] = convex.rejected_distance;
    metrics["rejected_plane_gain"] = convex.rejected_gain;
    metrics["rejected_work_limit"] = convex.rejected_work_limit;
    metrics["certificate_boxes"] = convex.certificate_boxes;
    metrics["max_volume_inflation"] = options.convex.max_volume_inflation;
    metrics["max_fill_distance_m"] = options.convex.max_fill_distance;
    metrics["max_candidate_voxels"] = options.convex.max_candidate_voxels;
    metrics["max_certificate_boxes"] = options.convex.max_certificate_boxes;
    metrics["occupied_volume_m3"] = convex.occupied_volume;
    metrics["output_volume_m3"] = convex.output_volume;
    metrics["volume_inflation"] = convex.volume_inflation;
    metrics["maximum_fill_distance_bound_m"] = convex.maximum_fill_distance_bound;
    metrics["raw_points_outside_tree"] = validation.raw_outside;
    metrics["false_negative"] = validation.source_outside;
    metrics["added_occupied_samples"] = validation.added;
    metrics["passed"] = validation.failures.empty();
    start = Clock::now();
    writeOutputs(options,octree,convex,tree,validation);
    metrics["output_geometry_io_ms"] = ms(start);
    auto json = file(options.output / "benchmark.json");
    json << "[{\"pipeline\":\"octree-hull\",\"approximation_mode\":\"hull\",\"length_unit\":\"m\",\"time_unit\":\"ms\",\"volume_reference\":\"occupied_voxel_union\"";
    auto csv = file(options.output / "benchmark.csv"); csv << "pipeline";
    std::cout << "Pipeline: octree-hull (explicit bounded approximation)\n";
    for(const auto& item : metrics)
    {
      json << ",\n\"" << item.first << "\":" << item.second;
      csv << ',' << item.first;
      std::cout << item.first << ": " << item.second << '\n';
    }
    json << "}]\n"; json.close();
    csv << "\noctree-hull"; for(const auto& item : metrics) csv << ',' << item.second;
    csv << '\n'; csv.close();
    return validation.failures.empty() ? EXIT_SUCCESS : EXIT_FAILURE;
  }
  catch(const std::exception& error)
  {
    std::cerr << "Octree hull error: " << error.what() << '\n'
              << "Usage: build_halfspace_tree input.xyz output --pipeline octree --voxel-size <m> --octree-approx hull "
                 "--max-volume-inflation <ratio> --max-fill-distance <m> "
                 "[--convex-max-voxels 128] [--convex-test-max-boxes 4096] [--validation-samples 1000]\n";
    return EXIT_FAILURE;
  }
}
}  // namespace rokae_demo
