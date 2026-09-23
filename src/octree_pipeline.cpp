#include <rokae_demo/octree_pipeline.hpp>
#include <rokae_demo/octree_occupancy.hpp>
#include <rokae_demo/octree_rectangular_pruning.hpp>
#include <rokae_demo/octree_boundary.hpp>
#include <rokae_demo/octree_convex_approximation.hpp>
#include <rokae_demo/octree_boundary_tree.hpp>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>

namespace rokae_demo
{
namespace
{
namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;
using Metrics = std::map<std::string, double>;

double milliseconds(Clock::time_point start)
{
  return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

struct Options
{
  fs::path input, output;
  double voxel_size = 0;
  unsigned max_depth = 20;
  std::size_t samples = 1000;
  std::string pruning_mode = "rect";
  bool boundary = false;
};

const char* pruningMode(const OctreeCellPartition& partition)
{
  return partition.rectangles ? "rect" : (partition.base.exact ? "exact" : "none");
}

std::size_t integer(const std::string& value)
{
  if(value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
    throw std::invalid_argument("Expected a nonnegative integer, got: " + value);
  const auto parsed = std::stoull(value);
  if(parsed > std::numeric_limits<std::size_t>::max())
    throw std::out_of_range("Integer parameter is too large");
  return static_cast<std::size_t>(parsed);
}

Options parse(int argc, char** argv)
{
  if(argc < 3)
    throw std::invalid_argument("Expected input.xyz and output directory");
  Options options;
  options.input = argv[1];
  options.output = argv[2];
  std::set<std::string> seen;
  for(int i = 3; i < argc; i += 2)
  {
    const std::string name = argv[i];
    if(i + 1 >= argc || !seen.insert(name).second)
      throw std::invalid_argument("Missing value or repeated option: " + name);
    const std::string value = argv[i + 1];
    if(name == "--pipeline")
    {
      if(value != "octree") throw std::invalid_argument("Expected --pipeline octree");
    }
    else if(name == "--voxel-size")
    {
      std::size_t consumed = 0;
      options.voxel_size = std::stod(value, &consumed);
      if(consumed != value.size() || !std::isfinite(options.voxel_size) || options.voxel_size <= 0)
        throw std::invalid_argument("voxel-size must be finite and positive (m)");
    }
    else if(name == "--max-depth")
    {
      const auto depth = integer(value);
      if(depth > 52) throw std::invalid_argument("max-depth must be in [0, 52]");
      options.max_depth = static_cast<unsigned>(depth);
    }
    else if(name == "--validation-samples") options.samples = integer(value);
    else if(name == "--octree-pruning")
    {
      if(value != "rect" && value != "exact" && value != "none")
        throw std::invalid_argument("octree-pruning must be rect, exact or none; approximate filling is not enabled");
      options.pruning_mode = value;
    }
    else if(name == "--octree-boundary")
    {
      if(value != "exact" && value != "none")
        throw std::invalid_argument("octree-boundary must be exact or none");
      options.boundary = value == "exact";
    }
    else if(name == "--octree-approx")
    {
      if(value != "none") throw std::invalid_argument("octree-approx must be none or hull");
    }
    else if(name == "--tree-source")
    {
      if(value != "volume") throw std::invalid_argument("tree-source must be volume or boundary");
    }
    else throw std::invalid_argument("Unsupported octree option: " + name);
  }
  if(options.voxel_size <= 0)
    throw std::invalid_argument("Octree requires explicit --voxel-size in meters");
  return options;
}

std::vector<VoxelPoint> readCloud(const fs::path& path)
{
  std::ifstream in(path);
  if(!in) throw std::runtime_error("Cannot read XYZ input: " + path.string());
  std::vector<VoxelPoint> points;
  std::string line;
  std::size_t line_number = 0;
  while(std::getline(in, line))
  {
    ++line_number;
    const auto comment = line.find('#');
    if(comment != std::string::npos) line.resize(comment);
    std::istringstream row(line);
    row >> std::ws;
    if(row.eof()) continue;
    double x, y, z;
    if(!(row >> x >> y >> z) || !std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
      throw std::runtime_error("Invalid/nonfinite XYZ at line " + std::to_string(line_number));
    row >> std::ws;
    if(!row.eof())
      throw std::runtime_error("Expected exactly 3 XYZ columns at line " + std::to_string(line_number));
    points.emplace_back(x, y, z);
  }
  if(in.bad()) throw std::runtime_error("XYZ input read failed");
  if(points.empty()) throw std::invalid_argument("Empty XYZ input");
  return points;
}

std::ofstream outputFile(const fs::path& path)
{
  std::ofstream out(path);
  if(!out) throw std::runtime_error("Cannot write output: " + path.string());
  out.exceptions(std::ios::badbit | std::ios::failbit);
  out << std::setprecision(17);
  return out;
}

template<class Range> void arrayJson(std::ostream& out, const Range& values)
{
  out << '[';
  bool first = true;
  for(const auto value : values)
  {
    if(!first) out << ',';
    first = false;
    out << value;
  }
  out << ']';
}

struct Validation : ClusterValidation
{
  std::size_t raw_outside_voxels = 0, raw_outside_tree = 0;
  std::size_t tested_points = 0, false_negative = 0, false_positive = 0;
  std::size_t cluster_mismatch = 0, tree_mismatch = 0, scalar_mismatch = 0;
  std::size_t pruning_tree_samples = 0, pruning_tree_mismatch = 0;
  double maximum_pruning_scalar_change = 0;
  double baseline_tree_build_ms = 0;
  bool occupied_union_certificate = false;
  bool pruning_certificate = false;
  bool rectangular_certificate = false;
  bool boundary_enabled = false, boundary_certificate = false;
};

Validation validate(const OctreeOccupancy& octree, const std::vector<VoxelPoint>& raw,
                    const OctreeCellPartition& partition,
                    std::vector<ConvexCluster>& clusters, const std::vector<HS::Vec3>& vertices,
                    const TreeRepresentation& tree, std::size_t samples)
{
  Validation result;
  validateOctreeOccupancy(octree, raw);
  validateOctreeCellPartition(octree, partition);
  result.pruning_certificate = true;
  result.rectangular_certificate = true;
  validateClusters(clusters, vertices, 0, result);
  // Exhaustive, constructive coverage certificate: exactly the six directed
  // planes of each full selected box, with every fine voxel accounted for once.
  result.occupied_union_certificate = clusters.size() == partition.cells.size();
  std::vector<unsigned char> covered(octree.leaf_nodes.size(), 0);
  for(std::size_t i = 0; i < clusters.size(); ++i)
  {
    const auto& c = clusters[i];
    const auto& node = partition.cells.at(i);
    const auto& box = node.box;
    bool valid = c.cells.size() == node.occupied_count && c.reduced_planes.size() == 6 &&
                 c.exact_planes.size() == 6 && c.vertices.size() == 8 && c.leaf_ids.size() == 6;
    for(const auto voxel : c.cells)
    {
      if(voxel >= covered.size() || covered[voxel] || partition.voxel_cell_ids[voxel] != i)
        valid = false;
      else covered[voxel] = 1;
    }
    for(unsigned j = 0; valid && j < 6; ++j)
    {
      const unsigned axis = j / 2;
      const double sign = j % 2 ? 1 : -1;
      const double offset = j % 2 ? box.upper[axis] : -box.lower[axis];
      const std::array<double, 3> n{axis == 0 ? sign : 0, axis == 1 ? sign : 0, axis == 2 ? sign : 0};
      const auto matches = [&](const auto& p) {
        return p.normal.x == n[0] && p.normal.y == n[1] && p.normal.z == n[2] && p.offset == offset;
      };
      valid = matches(c.exact_planes[j]) && matches(c.reduced_planes[j]) &&
              matches(tree.leaves.at(c.leaf_ids[j]));
    }
    result.occupied_union_certificate &= valid;
  }
  for(const auto seen : covered) result.occupied_union_certificate &= seen == 1;

  // Old/new *decision* consistency. Exact union-preserving pruning can change
  // the scalar field at removed internal faces, so scalar equality is NOT a
  // pruning invariant. Build the Phase 1 tree only as a validation reference.
  TreeRepresentation baseline;
  const bool compare_baseline = partition.cells.size() < octree.leaf_nodes.size() && samples > 0;
  if(compare_baseline)
  {
    const auto start = Clock::now();
    std::vector<HS::Vec3> fine_vertices;
    auto fine_clusters = octreeAabbClusters(octree, fine_vertices);
    baseline = buildClusterTree(fine_clusters, 0, "octree_aabb");
    result.baseline_tree_build_ms = milliseconds(start);
  }
  const auto check = [&](const HS::Vec3& p, bool is_raw, bool compare_pruning = false) {
    const bool occupied = octree.contains(p);
    const bool cluster = insideClusterUnion(clusters, p, 0);
    const double before = HS::evaluate(tree.raw, tree.leaves, p);
    const double after = HS::evaluate(tree.simplified, tree.leaves, p);
    const bool inside = after <= 0;
    ++result.tested_points;
    result.cluster_mismatch += occupied != cluster;
    result.false_negative += occupied && !inside;
    result.false_positive += !occupied && inside;
    result.tree_mismatch += (before <= 0) != inside;
    result.scalar_mismatch += before != after;
    if(compare_pruning)
    {
      const double old_value = compare_baseline ? HS::evaluate(baseline.raw, baseline.leaves, p) : before;
      ++result.pruning_tree_samples;
      result.pruning_tree_mismatch += (old_value <= 0) != inside;
      result.maximum_pruning_scalar_change =
          std::max(result.maximum_pruning_scalar_change, std::abs(old_value - after));
    }
    if(is_raw)
    {
      result.raw_outside_voxels += !occupied;
      result.raw_outside_tree += !(before <= 0) || !inside;
    }
  };
  for(const auto& p : raw) check({p.x(), p.y(), p.z()}, true);
  // Keep the same mandatory fine-grid probes before/after pruning, including
  // eliminated internal faces. Pruning must not reduce validation coverage.
  for(const int id : octree.leaf_nodes)
  {
    const auto& box = octree.nodes[id].box;
    for(unsigned mask = 0; mask < 8; ++mask)
      check({(mask & 1) ? box.upper[0] : box.lower[0],
             (mask & 2) ? box.upper[1] : box.lower[1],
             (mask & 4) ? box.upper[2] : box.lower[2]}, false);
    check({box.lower[0] + (box.upper[0] - box.lower[0]) / 2,
           box.lower[1] + (box.upper[1] - box.lower[1]) / 2,
           box.lower[2] + (box.upper[2] - box.lower[2]) / 2}, false);
  }
  std::mt19937_64 generator(20260907);
  const auto& bounds = octree.nodes.front().box;
  std::uniform_real_distribution<double> fraction(-0.05, 1.05);
  for(std::size_t i = 0; i < samples; ++i)
  {
    std::array<double, 3> p{};
    for(unsigned axis = 0; axis < 3; ++axis)
      p[axis] = bounds.lower[axis] + fraction(generator) * (bounds.upper[axis] - bounds.lower[axis]);
    check({p[0], p[1], p[2]}, false, true);
  }
  // Probe recovered-parent centers too: uniform random sampling may miss all
  // of a small merged region and conceal the scalar change on internal faces.
  if(compare_baseline)
    for(std::size_t i = 0; i < clusters.size(); ++i)
      if(partition.cells[i].occupied_count > 1)
        check(clusters[i].centroid, false, true);
  if(!result.occupied_union_certificate) result.failures.push_back("Occupied voxel coverage certificate failed");
  if(result.raw_outside_voxels || result.raw_outside_tree) result.failures.push_back("Raw point containment failed");
  if(result.cluster_mismatch || result.false_negative || result.false_positive)
    result.failures.push_back("Occupied geometry and halfspace tree disagree");
  if(result.tree_mismatch || result.scalar_mismatch) result.failures.push_back("Tree simplification changed evaluation");
  if(result.pruning_tree_mismatch) result.failures.push_back("Pruning changed tree inside/outside decision");
  return result;
}

void writeGeometry(const fs::path& output, const OctreeOccupancy& octree,
                   const OctreeCellPartition& partition,
                   const std::vector<ConvexCluster>& clusters, const TreeRepresentation& tree, bool boundary)
{
  const auto& pruning = partition.base;
  auto nodes = outputFile(output / "octree.json");
  nodes << "{\n\"pipeline\":\"octree\",\"phase\":" << (boundary ? 3 : (pruning.exact ? 2 : 1))
        << ",\"boundary_mode\":\"" << (boundary ? "exact" : "none") << '"'
        << ",\"pruning_mode\":\"" << pruningMode(partition)
        << "\",\"representation\":\"cubic_hierarchy_with_convex_cell_mapping"
        << "\",\"length_unit\":\"m\",\"voxel_size\":"
        << octree.voxel_size << ",\"required_depth\":" << octree.required_depth
        << ",\"source_octree_nodes\":" << octree.nodes.size()
        << ",\"occupied_voxel_indices\":[";
  for(std::size_t i = 0; i < octree.occupied_voxels.size(); ++i)
  {
    if(i) nodes << ',';
    arrayJson(nodes, octree.occupied_voxels[i]);
  }
  nodes << "],\"root_id\":0,\"nodes\":[\n";
  std::vector<int> cell_id(octree.nodes.size(), -1);
  for(std::size_t i = 0; i < pruning.cell_nodes.size(); ++i)
    cell_id[pruning.cell_nodes[i]] = static_cast<int>(i);
  std::vector<std::vector<std::size_t>> base_voxels(pruning.cell_nodes.size());
  for(std::size_t i = 0; i < pruning.voxel_cell_ids.size(); ++i)
    base_voxels[pruning.voxel_cell_ids[i]].push_back(i);
  bool first_node = true;
  for(std::size_t i = 0; i < octree.nodes.size(); ++i)
  {
    if(!pruning.active_nodes[i]) continue;
    if(!first_node) nodes << ",\n";
    first_node = false;
    const auto& n = octree.nodes[i];
    nodes << "{\"id\":" << i << ",\"lower_index\":";
    arrayJson(nodes, n.lower);
    nodes << ",\"width_voxels\":" << n.width << ",\"depth\":" << n.depth << ",\"lower\":";
    arrayJson(nodes, n.box.lower);
    nodes << ",\"upper\":"; arrayJson(nodes, n.box.upper);
    nodes << ",\"convex_cell_id\":"
          << (cell_id[i] < 0 ? -1 : static_cast<int>(partition.base_cell_ids[cell_id[i]])) << ",\"children\":";
    if(cell_id[i] >= 0) arrayJson(nodes, std::array<int, 8>{{-1,-1,-1,-1,-1,-1,-1,-1}});
    else arrayJson(nodes, n.children);
    nodes << ",\"occupied_count\":" << n.occupied_count << ",\"point_count\":" << n.point_count
          << ",\"point_indices\":";
    nodes << '[';
    bool first_point = true;
    if(cell_id[i] >= 0)
      for(const auto voxel : base_voxels[cell_id[i]])
        for(const auto point : octree.nodes[octree.leaf_nodes[voxel]].point_indices)
        {
          if(!first_point) nodes << ',';
          first_point = false;
          nodes << point;
        }
    nodes << "]}";
  }
  nodes << "\n]}\n";
  nodes.close();
  auto cells = outputFile(output / "convex_clusters.json");
  cells << "{\"length_unit\":\"m\",\"source\":\"occupied_voxel_aabb\",\"clusters\":[\n";
  for(std::size_t i = 0; i < clusters.size(); ++i)
  {
    const auto& node = partition.cells[i];
    cells << "{\"id\":" << clusters[i].id << ",\"octree_node\":" << node.source_node
          << ",\"source_parent_node\":" << node.parent_node << ",\"sibling_mask\":" << node.sibling_mask
          << ",\"voxel_index\":";
    arrayJson(cells, node.lower);
    cells << ",\"extent_voxels\":"; arrayJson(cells, node.extent);
    cells << ",\"width_voxels\":";
    if(node.extent[0] == node.extent[1] && node.extent[1] == node.extent[2]) cells << node.extent[0];
    else cells << "null";
    cells << ",\"occupied_count\":" << node.occupied_count;
    cells << ",\"source_voxel_ids\":"; arrayJson(cells, clusters[i].cells);
    cells << ",\"lower\":"; arrayJson(cells, node.box.lower);
    cells << ",\"upper\":"; arrayJson(cells, node.box.upper);
    cells << ",\"leaf_ids\":"; arrayJson(cells, clusters[i].leaf_ids);
    cells << '}' << (i + 1 == clusters.size() ? "\n" : ",\n");
  }
  cells << "]}\n";
  cells.close();
  writeTreeJson(output / "logic_tree_raw.json", tree, tree.raw, "raw", "octree");
  writeTreeJson(output / "logic_tree_simplified.json", tree, tree.simplified, "simplified", "octree");
  for(const auto& item : std::vector<std::pair<std::string, HS::ExpressionPtr>>{
          {"raw", tree.raw}, {"simplified", tree.simplified}})
  {
    auto out = outputFile(output / ("logic_expression_" + item.first + ".txt"));
    out << "inside(x) =\n" << HS::booleanExpression(item.second, 0)
        << "\n\nphi(x) =\n" << HS::scalarExpression(item.second, 0) << '\n';
    out.close();
  }
}

void writeValidation(const fs::path& output, const Validation& v)
{
  auto out = outputFile(output / "validation_report.json");
  out << std::boolalpha << "{\n\"pipeline\":\"octree\",\"passed\":" << v.failures.empty()
      << ",\"occupied_union_certificate\":" << v.occupied_union_certificate
      << ",\"pruning_certificate\":" << v.pruning_certificate
      << ",\"rectangular_partition_certificate\":" << v.rectangular_certificate
      << ",\"boundary_enabled\":" << v.boundary_enabled
      << ",\"boundary_certificate\":" << (v.boundary_enabled ? (v.boundary_certificate ? "true" : "false") : "null")
      << ",\"pruning_tree_comparison_samples\":" << v.pruning_tree_samples
      << ",\"pruning_tree_decision_mismatch\":" << v.pruning_tree_mismatch
      << ",\"maximum_pruning_scalar_change\":" << v.maximum_pruning_scalar_change
      << ",\"pruning_scalar_equality_required\":false"
      << ",\"baseline_tree_build_ms\":" << v.baseline_tree_build_ms
      << ",\"all_clusters_convex\":" << v.all_clusters_convex
      << ",\"all_cluster_vertices_satisfy\":" << v.all_cluster_vertices_satisfy
      << ",\"all_cluster_centroids_strictly_inside\":" << v.all_cluster_centroids_strictly_inside
      << ",\"classification_epsilon\":0,\"raw_points_outside_voxels\":" << v.raw_outside_voxels
      << ",\"raw_points_outside_tree\":" << v.raw_outside_tree
      << ",\"tested_points\":" << v.tested_points << ",\"false_negative\":" << v.false_negative
      << ",\"false_positive\":" << v.false_positive << ",\"cluster_mismatch\":" << v.cluster_mismatch
      << ",\"raw_simplified_mismatch\":" << v.tree_mismatch << ",\"scalar_mismatch\":" << v.scalar_mismatch
      << ",\"failures\":[";
  for(std::size_t i = 0; i < v.failures.size(); ++i)
    out << (i ? "," : "") << '"' << v.failures[i] << '"';
  out << "]}\n";
  out.close();
}

void writeMetrics(const fs::path& output, const Metrics& metrics, const OctreeCellPartition& partition, bool boundary)
{
  auto json = outputFile(output / "benchmark.json");
  json << "[{\"pipeline\":\"octree\",\"pruning_mode\":\"" << pruningMode(partition)
       << "\",\"boundary_mode\":\"" << (boundary ? "exact" : "none")
       << "\",\"length_unit\":\"m\",\"time_unit\":\"ms\","
          "\"volume_unit\":\"m^3\",\"volume_reference\":\"occupied_voxel_union\"";
  auto csv = outputFile(output / "benchmark.csv");
  csv << "pipeline";
  for(const auto& item : metrics)
  {
    json << ",\n\"" << item.first << "\":" << item.second;
    csv << ',' << item.first;
  }
  json << "}]\n";
  csv << "\noctree";
  for(const auto& item : metrics) csv << ',' << item.second;
  csv << '\n';
  json.close(); csv.close();
}
}  // namespace

int runOctreePipeline(int argc, char** argv)
{
  for(int i = 3; i+1 < argc; i += 2)
    if(std::string(argv[i]) == "--tree-source" && std::string(argv[i+1]) == "boundary")
      return runOctreeBoundaryTreePipeline(argc,argv);
  for(int i = 3; i+1 < argc; i += 2)
    if(std::string(argv[i]) == "--octree-approx" && std::string(argv[i+1]) == "hull")
      return runOctreeConvexPipeline(argc,argv);
  try
  {
    const auto options = parse(argc, argv);
    Metrics metrics;
    auto start = Clock::now();
    const auto raw = readCloud(options.input);
    metrics["input_io_ms"] = milliseconds(start);
    const auto total_start = Clock::now();
    start = Clock::now();
    auto octree = buildOctreeOccupancy(raw, options.voxel_size, options.max_depth);
    metrics["split_time"] = milliseconds(start);  // includes indexing/sorting occupancy
    start = Clock::now();
    auto base = pruneOctree(octree, options.pruning_mode != "none");
    const auto rectangular_start = Clock::now();
    const auto partition = partitionOctreeCells(octree, std::move(base), options.pruning_mode == "rect");
    metrics["rectangular_partition_ms"] = milliseconds(rectangular_start);
    metrics["pruning_time"] = milliseconds(start);  // includes constructing the cut in none mode
    const auto& pruning = partition.base;
    metrics["phase"] = options.boundary ? 3 : (pruning.exact ? 2 : 1);
    OctreeBoundary boundary;
    metrics["boundary_enabled"] = options.boundary;
    metrics["boundary_time"] = metrics["boundary_validation_ms"] = metrics["boundary_io_ms"] = 0;
    if(options.boundary)
    {
      start = Clock::now();
      boundary = extractOctreeBoundary(octree, partition);
      metrics["boundary_time"] = milliseconds(start);
      metrics["boundary_extraction_ms"] = boundary.extraction_ms;
      metrics["boundary_merge_ms"] = boundary.merge_ms;
      metrics["boundary_faces_before"] = boundary.voxel_faces_before;
      metrics["internal_faces_removed"] = boundary.internal_faces_removed;
      metrics["boundary_faces_exposed"] = boundary.faces.size();
      metrics["boundary_patches"] = boundary.patches.size();
      metrics["boundary_faces_merged"] = boundary.faces.size() - boundary.patches.size();
      metrics["boundary_triangles"] = 2 * boundary.faces.size();
    }
    start = Clock::now();
    std::vector<HS::Vec3> vertices;
    auto clusters = octreeAabbClusters(partition, vertices);
    metrics["plane_time"] = milliseconds(start);
    start = Clock::now();
    auto tree = buildClusterTree(clusters, 0, "octree_aabb");
    tree.simplified = HS::simplifyToFixedPoint(tree.raw);
    metrics["tree_nodes"] = HS::serializeTree(tree.simplified).nodes.size();
    metrics["raw_tree_nodes"] = HS::serializeTree(tree.raw).nodes.size();
    metrics["tree_time"] = milliseconds(start);
    metrics["core_runtime_ms"] = milliseconds(total_start);
    start = Clock::now();
    auto validation = validate(octree, raw, partition, clusters, vertices, tree, options.samples);
    validation.boundary_enabled = options.boundary;
    if(options.boundary)
    {
      const auto boundary_start = Clock::now();
      try
      {
        validateOctreeBoundary(octree, partition, boundary);
        validation.boundary_certificate = true;
      }
      catch(const std::exception&)
      {
        validation.failures.push_back("Boundary exact coverage certificate failed");
      }
      metrics["boundary_validation_ms"] = milliseconds(boundary_start);
    }
    metrics["validation_time"] = metrics["validation_ms"] = milliseconds(start);
    metrics["total_time"] = metrics["total_ms"] = milliseconds(total_start);
    metrics["raw_points"] = raw.size();
    metrics["occupied_voxels"] = metrics["leaf_nodes"] = octree.leaf_nodes.size();
    metrics["octree_nodes"] = octree.nodes.size();
    metrics["octree_nodes_after_prune"] =
        std::count(pruning.active_nodes.begin(), pruning.active_nodes.end(), static_cast<unsigned char>(1));
    metrics["convex_cells_before_prune"] = octree.leaf_nodes.size();
    metrics["output_leaf_nodes"] = clusters.size();
    metrics["pruning_merges"] = pruning.merge_operations;
    metrics["octree_leaf_nodes_after_prune"] = pruning.cell_nodes.size();
    metrics["convex_cells_before_rectangles"] = pruning.cell_nodes.size();
    metrics["rectangular_groups"] = partition.rectangular_groups;
    metrics["rectangular_merge_count"] = partition.rectangular_merge_count;
    metrics["rectangular_cells_removed"] = pruning.cell_nodes.size() - partition.cells.size();
    metrics["convex_cells"] = clusters.size();
    metrics["planes_before_prune"] = 6 * octree.leaf_nodes.size();
    metrics["planes_after_prune"] = 6 * clusters.size();
    // Exact reference count for the unpruned MIN-of-six-plane-MAX tree.
    metrics["tree_nodes_before_prune"] = 7 * octree.leaf_nodes.size() + (octree.leaf_nodes.size() > 1);
    metrics["planes_before_reduce"] = metrics["planes_after_reduce"] = 6 * clusters.size();
    metrics["unique_planes"] = tree.leaves.size();
    metrics["tree_leaf_references"] = tree.leaf_references;
    metrics["tree_depth"] = HS::treeDepth(tree.simplified);
    metrics["voxel_resolution"] = options.voxel_size;
    metrics["required_depth"] = octree.required_depth;
    metrics["max_depth"] = options.max_depth;
    metrics["validation_samples"] = options.samples;
    metrics["raw_points_outside_voxels"] = validation.raw_outside_voxels;
    metrics["raw_points_outside_tree"] = validation.raw_outside_tree;
    metrics["false_negative"] = validation.false_negative;
    metrics["false_positive"] = validation.false_positive;
    metrics["pruning_tree_decision_mismatch"] = validation.pruning_tree_mismatch;
    metrics["pruning_tree_comparison_samples"] = validation.pruning_tree_samples;
    metrics["maximum_pruning_scalar_change"] = validation.maximum_pruning_scalar_change;
    metrics["baseline_tree_build_ms"] = validation.baseline_tree_build_ms;
    double volume = 0;
    for(const auto id : octree.leaf_nodes) volume += octree.nodes[id].box.volume();
    if(!std::isfinite(volume) || volume <= 0) throw std::runtime_error("Total volume is not representable");
    double output_volume = 0;
    for(const auto& cell : partition.cells) output_volume += cell.box.volume();
    if(!std::isfinite(output_volume) || output_volume <= 0) throw std::runtime_error("Output volume is not representable");
    metrics["occupied_volume_m3"] = volume;
    metrics["output_volume_m3"] = output_volume;
    metrics["volume_roundoff_m3"] = output_volume - volume;
    // C = V by the exhaustive certificate above, not a point-cloud volume estimate.
    metrics["volume_inflation"] = 0;
    metrics["passed"] = validation.failures.empty();
    start = Clock::now();
    fs::create_directories(options.output);
    if(!options.boundary && (fs::exists(options.output / "boundary.off") ||
                             fs::exists(options.output / "boundary_patches.json")))
      std::cerr << "Boundary disabled: existing boundary files are from a previous run and are not updated. "
                   "Use a fresh output directory to avoid stale artifacts.\n";
    writeGeometry(options.output, octree, partition, clusters, tree, options.boundary);
    if(options.boundary && validation.boundary_certificate)
    {
      const auto boundary_start = Clock::now();
      writeOctreeBoundary(options.output, octree, boundary);
      metrics["boundary_io_ms"] = milliseconds(boundary_start);
    }
    writeValidation(options.output, validation);
    metrics["output_geometry_io_ms"] = milliseconds(start);
    writeMetrics(options.output, metrics, partition, options.boundary);
    std::cout << "Pipeline: octree (" << pruningMode(partition)
              << ")\n";
    for(const auto& item : metrics) std::cout << item.first << ": " << item.second << '\n';
    std::cout << "Output: " << options.output << '\n';
    return validation.failures.empty() ? EXIT_SUCCESS : EXIT_FAILURE;
  }
  catch(const std::exception& error)
  {
    std::cerr << "Octree error: " << error.what() << '\n'
              << "Usage: build_halfspace_tree input.xyz output --pipeline octree --voxel-size <m> "
                 "[--octree-pruning rect|exact|none] [--octree-boundary exact|none] "
                 "[--max-depth 20] [--validation-samples 1000]\n";
    return EXIT_FAILURE;
  }
}
}  // namespace rokae_demo
