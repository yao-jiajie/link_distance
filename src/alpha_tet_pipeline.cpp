#include <rokae_demo/alpha_tet_pipeline.hpp>
#include <rokae_demo/convex_cluster.hpp>

#include <rokae_demo/alpha_wrap_volume_adapter.hpp>
#include <rokae_demo/halfspace_logic_tree.hpp>
#include <rokae_demo/point_cloud_preprocessing.hpp>
#include <rokae_demo/units.hpp>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/IO/write_points.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/number_utils.h>
#include <CGAL/version.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace rokae_demo
{
namespace
{

namespace fs = std::filesystem;
namespace HS = rokae_demo::halfspace;
namespace PMP = CGAL::Polygon_mesh_processing;

using Adapter = AlphaWrapVolumeAdapter;
using Kernel = Adapter::Kernel;
using Point = Adapter::Point;
using Mesh = Adapter::SurfaceMesh;
using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using AabbTraits = CGAL::AABB_traits<Kernel, Primitive>;
using AabbTree = CGAL::AABB_tree<AabbTraits>;

#define ROKAE_STRINGIFY_IMPL(value) #value
#define ROKAE_STRINGIFY(value) ROKAE_STRINGIFY_IMPL(value)
constexpr const char* kCgalVersion = ROKAE_STRINGIFY(CGAL_VERSION);

enum class PlaneReductionMode
{
  EXACT,
  CONSERVATIVE
};

struct Options
{
  std::string input_path;
  std::string output_directory;
  double alpha = 0.0;
  double offset = 0.0;
  double alpha_ratio = 20.0;
  double offset_ratio = 100.0;
  bool absolute_parameters = false;
  bool voxel_parameters = false;
  double alpha_voxel_factor = 3.0;
  double offset_voxel_factor = 1.0;
  PlaneReductionMode plane_reduction = PlaneReductionMode::EXACT;
  double plane_angle_degrees = 0.1;
  double plane_distance_tolerance = 1e-8;
  double support_epsilon = 0.0;
  double classification_epsilon = 0.0;
  std::size_t validation_samples = 100000;
  std::size_t max_cluster_cells = 0;
  std::size_t max_cluster_planes = 0;
  std::size_t target_cluster_count = 0;
  std::vector<double> benchmark_alpha_ratios;
  PointCloudPreprocessOptions preprocessing;
  bool auto_offset = true;
  double offset_displacement_factor = 1.75;
  double offset_growth_factor = 1.5;
  std::size_t max_offset_retries = 4;
};


struct Candidate
{
  int a = -1;
  int b = -1;
  std::size_t version_a = 0;
  std::size_t version_b = 0;
  double shared_area = 0.0;
  std::size_t boundary_faces = 0;
  std::size_t support_planes = 0;
};

struct CandidatePriority
{
  bool operator()(const Candidate& lhs, const Candidate& rhs) const
  {
    if(lhs.shared_area != rhs.shared_area)
      return lhs.shared_area < rhs.shared_area;
    if(lhs.boundary_faces != rhs.boundary_faces)
      return lhs.boundary_faces > rhs.boundary_faces;
    if(lhs.support_planes != rhs.support_planes)
      return lhs.support_planes > rhs.support_planes;
    if(lhs.a != rhs.a)
      return lhs.a > rhs.a;
    return lhs.b > rhs.b;
  }
};

struct Metrics
{
  std::string pipeline = "alpha-tet";
  std::string parameter_mode;
  double alpha_voxel_factor = 0.0;
  double offset_voxel_factor = 0.0;
  std::size_t point_count = 0;
  double alpha = 0.0;
  double offset = 0.0;
  double requested_offset = 0.0;
  std::size_t offset_retries = 0;
  std::size_t wrap_faces = 0;
  std::size_t finite_cells = 0;
  std::size_t outside_cells = 0;
  std::size_t tetra_count = 0;
  std::size_t initial_clusters = 0;
  std::size_t final_clusters = 0;
  std::size_t raw_cluster_planes = 0;
  std::size_t reduced_cluster_planes = 0;
  std::size_t unique_planes = 0;
  std::size_t tetra_tree_leaf_references = 0;
  std::size_t tree_leaf_references = 0;
  std::size_t raw_tree_nodes = 0;
  std::size_t tree_nodes = 0;
  std::size_t tree_depth = 0;
  double alpha_wrap_ms = 0.0;
  double tet_extract_ms = 0.0;
  double agglomeration_ms = 0.0;
  double plane_reduction_ms = 0.0;
  double tree_build_ms = 0.0;
  double simplification_ms = 0.0;
  double validation_ms = 0.0;
  double total_ms = 0.0;
  double core_runtime_ms = 0.0;
  double containment_recovery_ms = 0.0;
  std::size_t raw_points_outside_wrap = 0;
  std::size_t raw_points_outside_tree = 0;
  std::size_t false_positive = 0;
  std::size_t false_negative = 0;
  double maximum_observed_expansion = 0.0;
  double maximum_support_expansion = 0.0;
  PointCloudPreprocessStats preprocessing;
};

struct Validation : ClusterValidation
{
  bool triangle_mesh = false;
  bool closed = false;
  bool self_intersects = true;
  bool bounds_volume = false;
  bool outward = false;
  // This count always refers to the unmodified, raw point cloud.
  std::size_t input_points_outside = 0;
  std::size_t wrap_input_points_outside = 0;
  std::size_t raw_points_outside_tree = 0;
  std::size_t tested_points = 0;
  std::size_t boundary_ignored = 0;
  std::size_t mesh_inside_tet_outside = 0;
  std::size_t mesh_outside_tet_inside = 0;
  std::size_t tet_inside_cluster_outside = 0;
  std::size_t tet_outside_cluster_inside = 0;
  std::size_t tetra_tree_mismatch = 0;
  std::size_t mesh_inside_tree_outside = 0;
  std::size_t mesh_outside_tree_inside = 0;
  std::size_t raw_simplified_mismatch = 0;
};


struct RunResult
{
  Metrics metrics;
  Validation validation;
  bool passed = false;
};

double elapsedMilliseconds(const std::chrono::steady_clock::time_point& start)
{
  return std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - start)
      .count();
}

HS::Vec3 toVec(const Point& point)
{
  return {CGAL::to_double(point.x()), CGAL::to_double(point.y()),
          CGAL::to_double(point.z())};
}

HS::Vec3 add(const HS::Vec3& lhs, const HS::Vec3& rhs)
{
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

HS::Vec3 subtract(const HS::Vec3& lhs, const HS::Vec3& rhs)
{
  return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

HS::Vec3 scale(const HS::Vec3& value, double factor)
{
  return {value.x * factor, value.y * factor, value.z * factor};
}

double dot(const HS::Vec3& lhs, const HS::Vec3& rhs)
{
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

HS::Vec3 cross(const HS::Vec3& lhs, const HS::Vec3& rhs)
{
  return {lhs.y * rhs.z - lhs.z * rhs.y,
          lhs.z * rhs.x - lhs.x * rhs.z,
          lhs.x * rhs.y - lhs.y * rhs.x};
}

double length(const HS::Vec3& value)
{
  return std::sqrt(dot(value, value));
}

bool parsePositiveDouble(const std::string& text, double& value)
{
  try
  {
    std::size_t parsed = 0;
    value = std::stod(text, &parsed);
    return parsed == text.size() && std::isfinite(value) && value > 0.0;
  }
  catch(const std::exception&)
  {
    return false;
  }
}

bool parseNonnegativeDouble(const std::string& text, double& value)
{
  try
  {
    std::size_t parsed = 0;
    value = std::stod(text, &parsed);
    return parsed == text.size() && std::isfinite(value) && value >= 0.0;
  }
  catch(const std::exception&)
  {
    return false;
  }
}

bool parseSize(const std::string& text, std::size_t& value)
{
  try
  {
    std::size_t parsed = 0;
    const unsigned long long parsed_value = std::stoull(text, &parsed);
    if(parsed != text.size())
      return false;
    value = static_cast<std::size_t>(parsed_value);
    return true;
  }
  catch(const std::exception&)
  {
    return false;
  }
}

bool parseRatioList(const std::string& text, std::vector<double>& values)
{
  std::stringstream stream(text);
  std::string token;
  while(std::getline(stream, token, ','))
  {
    double value = 0.0;
    if(!parsePositiveDouble(token, value))
      return false;
    values.push_back(value);
  }
  return !values.empty();
}

void printUsage(const char* program)
{
  std::cerr
      << "Usage: " << program
      << " points.{xyz,ply,off} output_dir --pipeline alpha-tet [options]\n"
      << "  Input coordinates and all length parameters use meters (m).\n"
      << "  --alpha VALUE_METERS --offset VALUE_METERS\n"
      << "  --alpha-ratio VALUE --offset-ratio VALUE\n"
      << "  --voxel-size METERS [--alpha-voxel-factor 3 --offset-voxel-factor 1]\n"
      << "  Absolute, voxel-factor, and bbox-ratio parameters are mutually exclusive.\n"
      << "  --plane-reduction exact|conservative\n"
      << "  --plane-angle-deg VALUE --plane-distance-tol VALUE\n"
      << "  --support-eps VALUE --classification-eps VALUE\n"
      << "  --max-cluster-cells COUNT --max-cluster-planes COUNT\n"
      << "  --target-cluster-count COUNT --validation-samples COUNT\n"
      << "  --preprocess none|jet|voxel --smooth-neighbors auto|COUNT\n"
      << "  --smooth-iterations COUNT\n"
      << "  --auto-offset true|false --offset-displacement-factor VALUE\n"
      << "  --offset-growth-factor VALUE --max-offset-retries COUNT\n"
      << "  --benchmark-alpha-ratios 20,15,10,8,5\n";
}

bool parseOptions(int argc, char** argv, Options& options)
{
  if(argc < 3)
    return false;
  options.input_path = argv[1];
  options.output_directory = argv[2];
  bool has_alpha = false;
  bool has_offset = false;
  bool has_alpha_ratio = false;
  bool has_offset_ratio = false;
  bool has_factors = false;
  bool has_preprocess = false;

  for(int i = 3; i < argc; ++i)
  {
    const std::string option = argv[i];
    if(i + 1 >= argc)
      return false;
    const std::string value = argv[++i];
    if(option == "--pipeline")
    {
      if(value != "alpha-tet" && value != "alpha_tet")
        return false;
    }
    else if(option == "--alpha")
    {
      if(!parsePositiveDouble(value, options.alpha))
        return false;
      has_alpha = true;
    }
    else if(option == "--offset")
    {
      if(!parsePositiveDouble(value, options.offset))
        return false;
      has_offset = true;
    }
    else if(option == "--alpha-ratio")
    {
      if(!parsePositiveDouble(value, options.alpha_ratio))
        return false;
      has_alpha_ratio = true;
    }
    else if(option == "--offset-ratio")
    {
      if(!parsePositiveDouble(value, options.offset_ratio))
        return false;
      has_offset_ratio = true;
    }
    else if(option == "--voxel-size")
    {
      if(!parsePositiveDouble(value, options.preprocessing.voxel_size))
        return false;
    }
    else if(option == "--alpha-voxel-factor" || option == "--offset-voxel-factor")
    {
      double& factor = option == "--alpha-voxel-factor"
                           ? options.alpha_voxel_factor : options.offset_voxel_factor;
      if(!parsePositiveDouble(value, factor))
        return false;
      has_factors = true;
    }
    else if(option == "--plane-reduction")
    {
      if(value == "exact")
        options.plane_reduction = PlaneReductionMode::EXACT;
      else if(value == "conservative")
        options.plane_reduction = PlaneReductionMode::CONSERVATIVE;
      else
        return false;
    }
    else if(option == "--plane-angle-deg")
    {
      if(!parseNonnegativeDouble(value, options.plane_angle_degrees) ||
         options.plane_angle_degrees >= 90.0)
        return false;
    }
    else if(option == "--plane-distance-tol")
    {
      if(!parseNonnegativeDouble(value, options.plane_distance_tolerance))
        return false;
    }
    else if(option == "--support-eps")
    {
      if(!parseNonnegativeDouble(value, options.support_epsilon))
        return false;
    }
    else if(option == "--classification-eps")
    {
      if(!parseNonnegativeDouble(value, options.classification_epsilon))
        return false;
    }
    else if(option == "--validation-samples")
    {
      if(!parseSize(value, options.validation_samples))
        return false;
    }
    else if(option == "--max-cluster-cells")
    {
      if(!parseSize(value, options.max_cluster_cells))
        return false;
    }
    else if(option == "--max-cluster-planes")
    {
      if(!parseSize(value, options.max_cluster_planes))
        return false;
    }
    else if(option == "--target-cluster-count")
    {
      if(!parseSize(value, options.target_cluster_count))
        return false;
    }
    else if(option == "--preprocess")
    {
      has_preprocess = true;
      if(value == "none")
        options.preprocessing.mode = PointCloudPreprocessMode::NONE;
      else if(value == "jet")
        options.preprocessing.mode = PointCloudPreprocessMode::JET;
      else if(value == "voxel")
        options.preprocessing.mode = PointCloudPreprocessMode::VOXEL;
      else
        return false;
    }
    else if(option == "--smooth-neighbors")
    {
      if(value == "auto")
        options.preprocessing.neighbors = 0;
      else if(!parseSize(value, options.preprocessing.neighbors) ||
              options.preprocessing.neighbors < 6)
        return false;
    }
    else if(option == "--smooth-iterations")
    {
      if(!parseSize(value, options.preprocessing.iterations) ||
         options.preprocessing.iterations == 0)
        return false;
    }
    else if(option == "--auto-offset")
    {
      if(value == "true")
        options.auto_offset = true;
      else if(value == "false")
        options.auto_offset = false;
      else
        return false;
    }
    else if(option == "--offset-displacement-factor")
    {
      if(!parsePositiveDouble(value, options.offset_displacement_factor))
        return false;
    }
    else if(option == "--offset-growth-factor")
    {
      if(!parsePositiveDouble(value, options.offset_growth_factor) ||
         options.offset_growth_factor <= 1.0)
        return false;
    }
    else if(option == "--max-offset-retries")
    {
      if(!parseSize(value, options.max_offset_retries))
        return false;
    }
    else if(option == "--benchmark-alpha-ratios")
    {
      if(!parseRatioList(value, options.benchmark_alpha_ratios))
        return false;
    }
    else
    {
      std::cerr << "Unknown alpha-tet option: " << option << "\n";
      return false;
    }
  }

  const bool any_absolute = has_alpha || has_offset;
  const bool sweep = !options.benchmark_alpha_ratios.empty();
  const bool any_relative = has_alpha_ratio || has_offset_ratio || sweep;
  if((any_absolute && any_relative) || (has_factors && (any_absolute || any_relative)))
  {
    std::cerr << "Absolute, voxel-factor, and bbox-ratio modes cannot be mixed.\n";
    return false;
  }
  if(any_absolute && !(has_alpha && has_offset))
  {
    std::cerr << "Absolute mode requires --alpha and --offset.\n";
    return false;
  }
  if(any_relative && !((has_alpha_ratio || sweep) && has_offset_ratio))
  {
    std::cerr << "BBox mode requires --alpha-ratio (or a sweep) and --offset-ratio.\n";
    return false;
  }
  const bool has_voxel = options.preprocessing.voxel_size > 0.0;
  if(has_voxel)
  {
    if(has_preprocess && options.preprocessing.mode != PointCloudPreprocessMode::VOXEL)
    {
      std::cerr << "--voxel-size cannot be combined with --preprocess none/jet.\n";
      return false;
    }
    options.preprocessing.mode = PointCloudPreprocessMode::VOXEL;
  }
  if(!has_voxel && (has_factors || options.preprocessing.mode == PointCloudPreprocessMode::VOXEL))
  {
    std::cerr << "Voxel mode requires --voxel-size in meters.\n";
    return false;
  }
  options.absolute_parameters = any_absolute;
  options.voxel_parameters = has_voxel && !any_absolute && !any_relative;
  return true;
}

FaceKey faceKey(const Adapter::Tetrahedron& tetrahedron, int opposite)
{
  FaceKey result;
  int cursor = 0;
  for(int i = 0; i < 4; ++i)
  {
    if(i != opposite)
      result.vertices[static_cast<std::size_t>(cursor++)] =
          tetrahedron.vertices[static_cast<std::size_t>(i)];
  }
  std::sort(result.vertices.begin(), result.vertices.end());
  return result;
}

EdgeKey edgeKey(std::size_t a, std::size_t b)
{
  if(b < a)
    std::swap(a, b);
  return {{{a, b}}};
}

Plane planeFromFace(const BoundaryFace& face,
                    const std::vector<HS::Vec3>& vertices)
{
  const HS::Vec3& a = vertices.at(face.key.vertices[0]);
  const HS::Vec3& b = vertices.at(face.key.vertices[1]);
  const HS::Vec3& c = vertices.at(face.key.vertices[2]);
  HS::Vec3 normal = cross(subtract(b, a), subtract(c, a));
  const double twice_area = length(normal);
  if(!std::isfinite(twice_area) ||
     twice_area <= std::numeric_limits<double>::epsilon())
    throw std::runtime_error("Degenerate tetrahedral boundary face");
  normal = scale(normal, 1.0 / twice_area);
  double offset = dot(normal, a);
  if(dot(normal, vertices.at(face.opposite_vertex)) - offset > 0.0)
  {
    normal = scale(normal, -1.0);
    offset = -offset;
  }
  Plane plane;
  plane.normal = normal;
  plane.point = a;
  plane.offset = offset;
  plane.area = 0.5 * twice_area;
  if(face.source_face_id >= 0)
    plane.source_face_ids.push_back(face.source_face_id);
  return plane;
}

bool samePlane(const Plane& lhs,
               const Plane& rhs,
               double cosine_tolerance,
               double distance_tolerance)
{
  return dot(lhs.normal, rhs.normal) >= cosine_tolerance &&
         std::abs(lhs.offset - rhs.offset) <= distance_tolerance;
}

std::size_t approximatePlaneCount(
    const std::map<FaceKey, BoundaryFace>& faces,
    const std::vector<HS::Vec3>& vertices,
    double distance_tolerance)
{
  std::vector<Plane> unique;
  for(const auto& entry : faces)
  {
    const Plane candidate = planeFromFace(entry.second, vertices);
    bool duplicate = false;
    for(const Plane& existing : unique)
    {
      if(samePlane(candidate, existing, 1.0 - 1e-12,
                   distance_tolerance))
      {
        duplicate = true;
        break;
      }
    }
    if(!duplicate)
      unique.push_back(candidate);
  }
  return unique.size();
}

class TetraClusterGraph
{
public:
  static std::vector<ConvexCluster> build(
      const std::vector<Adapter::Tetrahedron>& tetrahedra)
  {
    std::vector<ConvexCluster> clusters(tetrahedra.size());
    int source_face_id = 0;
    for(const Adapter::Tetrahedron& tetrahedron : tetrahedra)
    {
      ConvexCluster& cluster = clusters.at(tetrahedron.id);
      cluster.id = static_cast<int>(tetrahedron.id);
      cluster.cells.push_back(tetrahedron.id);
      cluster.vertices.insert(tetrahedron.vertices.begin(),
                              tetrahedron.vertices.end());
      for(int i = 0; i < 4; ++i)
      {
        BoundaryFace face;
        face.key = faceKey(tetrahedron, i);
        face.opposite_vertex =
            tetrahedron.vertices[static_cast<std::size_t>(i)];
        face.source_face_id = source_face_id++;
        cluster.boundary_faces.emplace(face.key, face);
        const int neighbor = tetrahedron.neighbors[static_cast<std::size_t>(i)];
        if(neighbor >= 0)
          cluster.neighbors.insert(neighbor);
      }
    }
    return clusters;
  }
};

class ConvexAgglomerator
{
  const std::vector<HS::Vec3>& vertices_;
  const Options& options_;
  double convexity_epsilon_ = 0.0;
  std::vector<ConvexCluster> clusters_;
  std::size_t active_count_ = 0;

  struct MergeData
  {
    std::vector<std::size_t> cells;
    std::map<FaceKey, BoundaryFace> boundary_faces;
    std::set<std::size_t> vertices;
    std::set<int> neighbors;
  };

public:
  ConvexAgglomerator(const std::vector<HS::Vec3>& vertices,
                     const std::vector<Adapter::Tetrahedron>& tetrahedra,
                     const Options& options,
                     double convexity_epsilon)
      : vertices_(vertices), options_(options),
        convexity_epsilon_(convexity_epsilon),
        clusters_(TetraClusterGraph::build(tetrahedra)),
        active_count_(tetrahedra.size())
  {
  }

  std::vector<ConvexCluster> run()
  {
    std::priority_queue<Candidate, std::vector<Candidate>,
                        CandidatePriority>
        queue;
    for(const ConvexCluster& cluster : clusters_)
    {
      for(const int neighbor : cluster.neighbors)
      {
        if(cluster.id < neighbor)
          queue.push(makeCandidate(cluster.id, neighbor));
      }
    }

    while(!queue.empty())
    {
      if(options_.target_cluster_count > 0 &&
         active_count_ <= options_.target_cluster_count)
        break;
      const Candidate candidate = queue.top();
      queue.pop();
      if(!candidateStillValid(candidate))
        continue;

      MergeData merged;
      if(!buildMerge(candidate.a, candidate.b, merged) ||
         !isConvex(merged.boundary_faces, merged.vertices))
        continue;

      commitMerge(candidate.a, candidate.b, std::move(merged));
      for(const int neighbor : clusters_[candidate.a].neighbors)
        queue.push(makeCandidate(candidate.a, neighbor));
    }

    std::vector<ConvexCluster> result;
    result.reserve(active_count_);
    for(ConvexCluster& cluster : clusters_)
    {
      if(cluster.active)
        result.push_back(std::move(cluster));
    }
    std::sort(result.begin(), result.end(),
              [](const ConvexCluster& lhs, const ConvexCluster& rhs) {
                return lhs.id < rhs.id;
              });
    std::map<int, int> remapped_ids;
    for(std::size_t i = 0; i < result.size(); ++i)
      remapped_ids.emplace(result[i].id, static_cast<int>(i));
    for(std::size_t i = 0; i < result.size(); ++i)
    {
      std::set<int> remapped_neighbors;
      for(const int neighbor : result[i].neighbors)
      {
        const auto remapped = remapped_ids.find(neighbor);
        if(remapped != remapped_ids.end())
          remapped_neighbors.insert(remapped->second);
      }
      result[i].neighbors = std::move(remapped_neighbors);
      result[i].id = static_cast<int>(i);
    }
    return result;
  }

private:
  Candidate makeCandidate(int first, int second) const
  {
    if(second < first)
      std::swap(first, second);
    const ConvexCluster& a = clusters_.at(static_cast<std::size_t>(first));
    const ConvexCluster& b = clusters_.at(static_cast<std::size_t>(second));
    Candidate candidate;
    candidate.a = first;
    candidate.b = second;
    candidate.version_a = a.version;
    candidate.version_b = b.version;
    for(const auto& entry : a.boundary_faces)
    {
      if(b.boundary_faces.count(entry.first) != 0)
        candidate.shared_area += planeFromFace(entry.second, vertices_).area;
    }
    candidate.boundary_faces = a.boundary_faces.size() +
                               b.boundary_faces.size();
    std::size_t shared_count = 0;
    for(const auto& entry : a.boundary_faces)
      shared_count += b.boundary_faces.count(entry.first);
    candidate.boundary_faces -= 2 * shared_count;
    std::map<FaceKey, BoundaryFace> merged_faces = a.boundary_faces;
    for(const auto& entry : b.boundary_faces)
    {
      const auto existing = merged_faces.find(entry.first);
      if(existing == merged_faces.end())
        merged_faces.emplace(entry);
      else
        merged_faces.erase(existing);
    }
    candidate.support_planes = approximatePlaneCount(
        merged_faces, vertices_, convexity_epsilon_);
    return candidate;
  }

  bool candidateStillValid(const Candidate& candidate) const
  {
    const ConvexCluster& a = clusters_.at(candidate.a);
    const ConvexCluster& b = clusters_.at(candidate.b);
    return a.active && b.active && a.version == candidate.version_a &&
           b.version == candidate.version_b &&
           a.neighbors.count(candidate.b) != 0;
  }

  bool buildMerge(int first, int second, MergeData& merged) const
  {
    const ConvexCluster& a = clusters_.at(first);
    const ConvexCluster& b = clusters_.at(second);
    if(options_.max_cluster_cells > 0 &&
       a.cells.size() + b.cells.size() > options_.max_cluster_cells)
      return false;

    merged.cells = a.cells;
    merged.cells.insert(merged.cells.end(), b.cells.begin(), b.cells.end());
    merged.vertices = a.vertices;
    merged.vertices.insert(b.vertices.begin(), b.vertices.end());
    merged.boundary_faces = a.boundary_faces;
    for(const auto& entry : b.boundary_faces)
    {
      const auto existing = merged.boundary_faces.find(entry.first);
      if(existing == merged.boundary_faces.end())
        merged.boundary_faces.emplace(entry);
      else
        merged.boundary_faces.erase(existing);
    }
    merged.neighbors = a.neighbors;
    merged.neighbors.insert(b.neighbors.begin(), b.neighbors.end());
    merged.neighbors.erase(first);
    merged.neighbors.erase(second);

    if(options_.max_cluster_planes > 0 &&
       approximatePlaneCount(merged.boundary_faces, vertices_,
                             convexity_epsilon_) >
           options_.max_cluster_planes)
      return false;
    return true;
  }

  bool isConvex(const std::map<FaceKey, BoundaryFace>& boundary_faces,
                const std::set<std::size_t>& vertices) const
  {
    for(const auto& entry : boundary_faces)
    {
      const BoundaryFace& face = entry.second;
      const HS::Vec3& a = vertices_.at(face.key.vertices[0]);
      const HS::Vec3& b = vertices_.at(face.key.vertices[1]);
      const HS::Vec3& c = vertices_.at(face.key.vertices[2]);
      const HS::Vec3& opposite = vertices_.at(face.opposite_vertex);
      const Point p0(a.x, a.y, a.z);
      const Point p1(b.x, b.y, b.z);
      const Point p2(c.x, c.y, c.z);
      const CGAL::Orientation interior_side = CGAL::orientation(
          p0, p1, p2, Point(opposite.x, opposite.y, opposite.z));
      if(interior_side == CGAL::COPLANAR)
        return false;
      for(const std::size_t vertex : vertices)
      {
        const HS::Vec3& value = vertices_.at(vertex);
        const CGAL::Orientation side = CGAL::orientation(
            p0, p1, p2, Point(value.x, value.y, value.z));
        if(side != CGAL::COPLANAR && side != interior_side)
          return false;
      }

      // Retain a scaled numerical guard because the exported halfspaces are
      // doubles even though the convexity predicate above is exact.
      const Plane plane = planeFromFace(entry.second, vertices_);
      for(const std::size_t vertex : vertices)
      {
        if(plane.evaluate(vertices_.at(vertex)) > convexity_epsilon_)
          return false;
      }
    }
    return true;
  }

  void commitMerge(int first, int second, MergeData merged)
  {
    ConvexCluster& a = clusters_.at(first);
    ConvexCluster& b = clusters_.at(second);
    for(const int neighbor : merged.neighbors)
    {
      ConvexCluster& adjacent = clusters_.at(neighbor);
      adjacent.neighbors.erase(first);
      adjacent.neighbors.erase(second);
      adjacent.neighbors.insert(first);
    }
    a.cells = std::move(merged.cells);
    a.boundary_faces = std::move(merged.boundary_faces);
    a.vertices = std::move(merged.vertices);
    a.neighbors = std::move(merged.neighbors);
    ++a.version;
    b.active = false;
    b.neighbors.clear();
    ++b.version;
    --active_count_;
  }
};

class DisjointSet
{
  std::vector<std::size_t> parent_;

public:
  explicit DisjointSet(std::size_t size) : parent_(size)
  {
    std::iota(parent_.begin(), parent_.end(), 0);
  }

  std::size_t find(std::size_t value)
  {
    if(parent_[value] != value)
      parent_[value] = find(parent_[value]);
    return parent_[value];
  }

  void unite(std::size_t a, std::size_t b)
  {
    a = find(a);
    b = find(b);
    if(a != b)
      parent_[b] = a;
  }
};

class SupportPlaneReducer
{
  const std::vector<HS::Vec3>& vertices_;
  const Options& options_;
  double support_epsilon_ = 0.0;

public:
  SupportPlaneReducer(const std::vector<HS::Vec3>& vertices,
                      const Options& options,
                      double support_epsilon)
      : vertices_(vertices), options_(options),
        support_epsilon_(support_epsilon)
  {
  }

  std::vector<Plane> extractExactPlanes(const ConvexCluster& cluster) const
  {
    std::vector<Plane> planes;
    planes.reserve(cluster.boundary_faces.size());
    for(const auto& entry : cluster.boundary_faces)
      planes.push_back(planeFromFace(entry.second, vertices_));
    return planes;
  }

  std::vector<Plane> reduce(const ConvexCluster& cluster,
                            double& maximum_expansion) const
  {
    std::vector<BoundaryFace> faces;
    std::vector<Plane> planes;
    for(const auto& entry : cluster.boundary_faces)
    {
      faces.push_back(entry.second);
      planes.push_back(planeFromFace(entry.second, vertices_));
    }
    DisjointSet sets(faces.size());
    std::map<EdgeKey, std::vector<std::size_t>> edge_faces;
    for(std::size_t i = 0; i < faces.size(); ++i)
    {
      const auto& v = faces[i].key.vertices;
      edge_faces[edgeKey(v[0], v[1])].push_back(i);
      edge_faces[edgeKey(v[1], v[2])].push_back(i);
      edge_faces[edgeKey(v[0], v[2])].push_back(i);
    }

    const double pi = std::acos(-1.0);
    const double cosine_tolerance =
        std::cos(options_.plane_angle_degrees * pi / 180.0);
    for(const auto& entry : edge_faces)
    {
      const std::vector<std::size_t>& adjacent = entry.second;
      for(std::size_t i = 0; i < adjacent.size(); ++i)
      {
        for(std::size_t j = i + 1; j < adjacent.size(); ++j)
        {
          const std::size_t a = adjacent[i];
          const std::size_t b = adjacent[j];
          bool merge = false;
          if(options_.plane_reduction == PlaneReductionMode::EXACT)
          {
            const auto& fa = faces[a].key.vertices;
            const auto& fb = faces[b].key.vertices;
            const Point p0(vertices_[fa[0]].x, vertices_[fa[0]].y,
                           vertices_[fa[0]].z);
            const Point p1(vertices_[fa[1]].x, vertices_[fa[1]].y,
                           vertices_[fa[1]].z);
            const Point p2(vertices_[fa[2]].x, vertices_[fa[2]].y,
                           vertices_[fa[2]].z);
            merge = dot(planes[a].normal, planes[b].normal) > 0.0;
            for(const std::size_t vertex : fb)
            {
              const HS::Vec3& q = vertices_[vertex];
              merge &= CGAL::coplanar(p0, p1, p2,
                                      Point(q.x, q.y, q.z));
            }
          }
          else
          {
            merge = samePlane(planes[a], planes[b], cosine_tolerance,
                              options_.plane_distance_tolerance);
          }
          if(merge)
            sets.unite(a, b);
        }
      }
    }

    std::map<std::size_t, std::vector<std::size_t>> groups;
    for(std::size_t i = 0; i < faces.size(); ++i)
      groups[sets.find(i)].push_back(i);

    std::vector<Plane> result;
    result.reserve(groups.size());
    for(const auto& entry : groups)
    {
      const std::vector<std::size_t>& group = entry.second;
      if(options_.plane_reduction == PlaneReductionMode::EXACT)
      {
        Plane representative = planes[group.front()];
        representative.source_face_ids.clear();
        for(const std::size_t index : group)
          representative.source_face_ids.insert(
              representative.source_face_ids.end(),
              planes[index].source_face_ids.begin(),
              planes[index].source_face_ids.end());
        result.push_back(std::move(representative));
        continue;
      }

      HS::Vec3 average;
      std::set<std::size_t> patch_vertices;
      Plane representative;
      for(const std::size_t index : group)
      {
        average = add(average,
                      scale(planes[index].normal, planes[index].area));
        patch_vertices.insert(faces[index].key.vertices.begin(),
                              faces[index].key.vertices.end());
        representative.source_face_ids.insert(
            representative.source_face_ids.end(),
            planes[index].source_face_ids.begin(),
            planes[index].source_face_ids.end());
      }
      const double average_length = length(average);
      if(average_length <= std::numeric_limits<double>::epsilon())
        average = planes[group.front()].normal;
      else
        average = scale(average, 1.0 / average_length);
      representative.normal = average;
      representative.area = 0.0;
      double patch_support = -std::numeric_limits<double>::infinity();
      for(const std::size_t vertex : patch_vertices)
        patch_support = std::max(patch_support,
                                 dot(average, vertices_[vertex]));
      double cluster_support = -std::numeric_limits<double>::infinity();
      for(const std::size_t vertex : cluster.vertices)
        cluster_support = std::max(cluster_support,
                                   dot(average, vertices_[vertex]));
      representative.offset = cluster_support + support_epsilon_;
      representative.point = scale(average, representative.offset);
      maximum_expansion =
          std::max(maximum_expansion,
                   representative.offset - patch_support);
      result.push_back(std::move(representative));
    }
    return result;
  }
};



TreeRepresentation buildTetraTree(
    const std::vector<Adapter::Tetrahedron>& tetrahedra,
    const std::vector<HS::Vec3>& vertices,
    double plane_tolerance)
{
  TreeRepresentation tree;
  std::vector<HS::ExpressionPtr> tetra_expressions;
  int source_id = 0;
  for(const Adapter::Tetrahedron& tetrahedron : tetrahedra)
  {
    std::vector<HS::ExpressionPtr> halfspaces;
    for(int i = 0; i < 4; ++i)
    {
      BoundaryFace face;
      face.key = faceKey(tetrahedron, i);
      face.opposite_vertex = tetrahedron.vertices[i];
      face.source_face_id = source_id++;
      const Plane plane = planeFromFace(face, vertices);
      const int leaf_id = findOrAddLeaf(
          plane, static_cast<int>(tetrahedron.id), plane_tolerance,
          tree.leaves);
      halfspaces.push_back(HS::makeLeaf(leaf_id));
      ++tree.leaf_references;
    }
    tetra_expressions.push_back(
        HS::makeNode(HS::TreeOperator::MAX, std::move(halfspaces)));
  }
  if(tetra_expressions.empty())
    throw std::runtime_error("No retained tetrahedra for raw tetra tree");
  tree.raw = tetra_expressions.size() == 1
                 ? tetra_expressions.front()
                 : HS::makeNode(HS::TreeOperator::MIN,
                                std::move(tetra_expressions));
  tree.simplified = tree.raw;
  return tree;
}


std::array<double, 6> pointBoundingBox(const std::vector<Point>& points)
{
  std::array<double, 6> bbox{
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity()};
  for(const Point& point : points)
  {
    bbox[0] = std::min(bbox[0], CGAL::to_double(point.x()));
    bbox[1] = std::min(bbox[1], CGAL::to_double(point.y()));
    bbox[2] = std::min(bbox[2], CGAL::to_double(point.z()));
    bbox[3] = std::max(bbox[3], CGAL::to_double(point.x()));
    bbox[4] = std::max(bbox[4], CGAL::to_double(point.y()));
    bbox[5] = std::max(bbox[5], CGAL::to_double(point.z()));
  }
  return bbox;
}

double bboxDiagonal(const std::array<double, 6>& bbox)
{
  const double dx = bbox[3] - bbox[0];
  const double dy = bbox[4] - bbox[1];
  const double dz = bbox[5] - bbox[2];
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}


std::size_t countOutsidePoints(const Mesh& mesh,
                               const std::vector<Point>& points)
{
  CGAL::Side_of_triangle_mesh<Mesh, Kernel> side(mesh);
  return static_cast<std::size_t>(std::count_if(
      points.begin(), points.end(), [&](const Point& point) {
        return side(point) == CGAL::ON_UNBOUNDED_SIDE;
      }));
}

void validateGeometry(const std::vector<Point>& input_points,
                      const std::vector<Point>& wrap_input_points,
                      const Adapter& adapter,
                      const std::vector<ConvexCluster>& clusters,
                      const TreeRepresentation& tetra_tree,
                      const TreeRepresentation& cluster_tree,
                      const Options& options,
                      const std::array<double, 6>& bbox,
                      double diagonal,
                      Validation& validation,
                      Metrics& metrics)
{
  const Mesh& mesh = adapter.surfaceMesh();
  validation.triangle_mesh = CGAL::is_triangle_mesh(mesh);
  validation.closed = CGAL::is_closed(mesh);
  validation.self_intersects = PMP::does_self_intersect(mesh);
  validation.bounds_volume = PMP::does_bound_a_volume(mesh);
  validation.outward = PMP::is_outward_oriented(mesh);
  if(!validation.triangle_mesh || !validation.closed ||
     validation.self_intersects || !validation.bounds_volume ||
     !validation.outward)
    validation.failures.push_back("Alpha Wrap surface validity failed");

  CGAL::Side_of_triangle_mesh<Mesh, Kernel> side(mesh);
  for(const Point& point : input_points)
  {
    if(side(point) == CGAL::ON_UNBOUNDED_SIDE)
      ++validation.input_points_outside;
    if(!HS::evaluateInside(cluster_tree.simplified, cluster_tree.leaves,
                           toVec(point), options.classification_epsilon))
      ++validation.raw_points_outside_tree;
  }
  metrics.raw_points_outside_wrap = validation.input_points_outside;
  metrics.raw_points_outside_tree = validation.raw_points_outside_tree;
  if(validation.raw_points_outside_tree != 0)
    validation.failures.push_back("Raw input points lie outside final expression tree");
  if(validation.input_points_outside != 0)
    validation.failures.push_back("Raw input points lie outside Alpha Wrap");
  for(const Point& point : wrap_input_points)
  {
    if(side(point) == CGAL::ON_UNBOUNDED_SIDE)
      ++validation.wrap_input_points_outside;
  }
  if(validation.wrap_input_points_outside != 0)
    validation.failures.push_back(
        "Preprocessed Alpha Wrap input points lie outside Alpha Wrap");

  AabbTree distance_tree(faces(mesh).first, faces(mesh).second, mesh);
  distance_tree.accelerate_distance_queries();
  const double expansion = 0.1 * diagonal;
  std::mt19937_64 generator(20260904);
  std::uniform_real_distribution<double> random_x(bbox[0] - expansion,
                                                   bbox[3] + expansion);
  std::uniform_real_distribution<double> random_y(bbox[1] - expansion,
                                                   bbox[4] + expansion);
  std::uniform_real_distribution<double> random_z(bbox[2] - expansion,
                                                   bbox[5] + expansion);
  const double epsilon_squared = options.classification_epsilon *
                                 options.classification_epsilon;
  const std::size_t tetra_tree_sample_limit =
      std::min<std::size_t>(options.validation_samples, 10000);
  for(std::size_t i = 0; i < options.validation_samples; ++i)
  {
    const Point point(random_x(generator), random_y(generator),
                      random_z(generator));
    const HS::Vec3 value = toVec(point);
    const CGAL::Bounded_side mesh_side = side(point);
    if(mesh_side == CGAL::ON_BOUNDARY ||
       distance_tree.squared_distance(point) <= epsilon_squared)
    {
      ++validation.boundary_ignored;
      continue;
    }

    const bool mesh_inside = mesh_side == CGAL::ON_BOUNDED_SIDE;
    const bool tetra_inside = adapter.containsInRetainedCells(point);
    const bool cluster_inside = insideClusterUnion(
        clusters, value, options.classification_epsilon);
    const bool tree_inside = HS::evaluateInside(
        cluster_tree.raw, cluster_tree.leaves, value,
        options.classification_epsilon);
    const bool simplified_inside = HS::evaluateInside(
        cluster_tree.simplified, cluster_tree.leaves, value,
        options.classification_epsilon);

    if(mesh_inside && !tetra_inside)
      ++validation.mesh_inside_tet_outside;
    if(!mesh_inside && tetra_inside)
      ++validation.mesh_outside_tet_inside;
    if(tetra_inside && !cluster_inside)
      ++validation.tet_inside_cluster_outside;
    if(!tetra_inside && cluster_inside)
      ++validation.tet_outside_cluster_inside;
    if(mesh_inside && !tree_inside)
      ++validation.mesh_inside_tree_outside;
    if(!mesh_inside && tree_inside)
    {
      ++validation.mesh_outside_tree_inside;
      metrics.maximum_observed_expansion =
          std::max(metrics.maximum_observed_expansion,
                   std::sqrt(CGAL::to_double(
                       distance_tree.squared_distance(point))));
    }
    if(tree_inside != simplified_inside)
      ++validation.raw_simplified_mismatch;
    if(i < tetra_tree_sample_limit)
    {
      const bool tetra_tree_inside = HS::evaluateInside(
          tetra_tree.raw, tetra_tree.leaves, value,
          options.classification_epsilon);
      if(tetra_tree_inside != tetra_inside)
        ++validation.tetra_tree_mismatch;
    }
    ++validation.tested_points;
  }

  metrics.false_positive = validation.mesh_outside_tree_inside;
  metrics.false_negative = validation.mesh_inside_tree_outside;
  if(validation.mesh_inside_tet_outside != 0 ||
     validation.mesh_outside_tet_inside != 0)
    validation.failures.push_back("Mesh/tetra union random validation failed");
  if(validation.tet_inside_cluster_outside != 0 ||
     validation.tet_outside_cluster_inside != 0)
    validation.failures.push_back("Tetra/cluster union validation failed");
  if(validation.tetra_tree_mismatch != 0)
    validation.failures.push_back("Raw tetra DNF validation failed");
  if(validation.raw_simplified_mismatch != 0)
    validation.failures.push_back("Raw/simplified tree validation failed");
  if(validation.mesh_inside_tree_outside != 0)
    validation.failures.push_back("Final tree has false negatives");
  if(options.plane_reduction == PlaneReductionMode::EXACT &&
     validation.mesh_outside_tree_inside != 0)
    validation.failures.push_back("Exact final tree has false positives");
}

Mesh clusterMesh(const ConvexCluster& cluster,
                 const std::vector<HS::Vec3>& vertices)
{
  Mesh mesh;
  std::map<std::size_t, Mesh::Vertex_index> vertex_map;
  for(const std::size_t vertex : cluster.vertices)
  {
    const HS::Vec3& point = vertices.at(vertex);
    vertex_map.emplace(vertex,
                       mesh.add_vertex(Point(point.x, point.y, point.z)));
  }
  for(const auto& entry : cluster.boundary_faces)
  {
    std::array<std::size_t, 3> ids = entry.first.vertices;
    const Plane plane = planeFromFace(entry.second, vertices);
    const HS::Vec3 ab = subtract(vertices[ids[1]], vertices[ids[0]]);
    const HS::Vec3 ac = subtract(vertices[ids[2]], vertices[ids[0]]);
    if(dot(cross(ab, ac), plane.normal) < 0.0)
      std::swap(ids[1], ids[2]);
    const auto face = mesh.add_face(vertex_map.at(ids[0]),
                                    vertex_map.at(ids[1]),
                                    vertex_map.at(ids[2]));
    if(face == Mesh::null_face())
      throw std::runtime_error("Failed to construct convex cluster mesh");
  }
  return mesh;
}



void writeValidationJson(const fs::path& path,
                         const Options& options,
                         const Metrics& metrics,
                         const Validation& validation)
{
  std::ofstream out(path);
  if(!out)
    throw std::runtime_error("Failed to write validation JSON");
  const auto boolean = [](bool value) { return value ? "true" : "false"; };
  out << std::setprecision(17)
      << "{\n  \"cgal_version\": \"" << kCgalVersion << "\",\n"
      << "  \"alpha_wrap_api\": \"CGAL internal API isolated by AlphaWrapVolumeAdapter\",\n"
      << "  \"pipeline\": \"alpha-tet\",\n"
      << "  \"units\": {\"length\": \"" << kLengthUnit
      << "\", \"area\": \"" << kAreaUnit << "\", \"volume\": \""
      << kVolumeUnit << "\"},\n"
      << "  \"preprocessing\": {\"mode\": \""
      << pointCloudPreprocessModeName(metrics.preprocessing.mode)
      << "\", \"input_points\": " << metrics.preprocessing.input_points
      << ", \"output_points\": " << metrics.preprocessing.output_points
      << ", \"requested_neighbors\": "
      << metrics.preprocessing.requested_neighbors
      << ", \"effective_neighbors\": "
      << metrics.preprocessing.effective_neighbors
      << ", \"iterations\": " << metrics.preprocessing.iterations
      << ", \"average_spacing\": "
      << metrics.preprocessing.average_spacing
      << ", \"mean_displacement\": "
      << metrics.preprocessing.mean_displacement
      << ", \"rms_displacement\": "
      << metrics.preprocessing.rms_displacement
      << ", \"percentile95_displacement\": "
      << metrics.preprocessing.percentile95_displacement
      << ", \"maximum_displacement\": "
      << metrics.preprocessing.maximum_displacement
      << ", \"scale_estimation_ms\": "
      << metrics.preprocessing.scale_estimation_ms
      << ", \"smoothing_ms\": " << metrics.preprocessing.smoothing_ms
      << ", \"total_ms\": " << metrics.preprocessing.total_ms << "},\n"
      << "  \"voxelization\": {\"voxel_size\": " << metrics.preprocessing.voxelization.voxel_size
      << ", \"occupied_voxels\": " << metrics.preprocessing.voxelization.occupied_voxels
      << ", \"voxel_points\": " << metrics.preprocessing.voxelization.voxel_points
      << ", \"raw_points_outside_voxels\": "
      << metrics.preprocessing.voxelization.raw_points_outside_voxels << "},\n"
      << "  \"raw_points_outside_wrap\": " << metrics.raw_points_outside_wrap << ",\n"
      << "  \"raw_points_outside_tree\": " << metrics.raw_points_outside_tree << ",\n"
      << "  \"offset_policy\": {\"automatic\": "
      << boolean(options.auto_offset) << ", \"requested_offset\": "
      << metrics.requested_offset << ", \"final_offset\": "
      << metrics.offset << ", \"displacement_factor\": "
      << options.offset_displacement_factor << ", \"growth_factor\": "
      << options.offset_growth_factor << ", \"retries\": "
      << metrics.offset_retries << "},\n"
      << "  \"plane_reduction\": \""
      << (options.plane_reduction == PlaneReductionMode::EXACT
              ? "exact"
              : "conservative")
      << "\",\n"
      << "  \"surface\": {\"triangle_mesh\": "
      << boolean(validation.triangle_mesh) << ", \"closed\": "
      << boolean(validation.closed) << ", \"self_intersects\": "
      << boolean(validation.self_intersects) << ", \"bounds_volume\": "
      << boolean(validation.bounds_volume) << ", \"outward\": "
      << boolean(validation.outward) << ", \"input_points_outside\": "
      << validation.input_points_outside
      << ", \"wrap_input_points_outside\": "
      << validation.wrap_input_points_outside << "},\n"
      << "  \"clusters\": {\"all_convex\": "
      << boolean(validation.all_clusters_convex)
      << ", \"all_vertices_satisfy\": "
      << boolean(validation.all_cluster_vertices_satisfy)
      << ", \"centroids_strictly_inside\": "
      << boolean(validation.all_cluster_centroids_strictly_inside)
      << "},\n"
      << "  \"classification\": {\"tested_points\": "
      << validation.tested_points << ", \"boundary_ignored\": "
      << validation.boundary_ignored
      << ", \"mesh_inside_tet_outside\": "
      << validation.mesh_inside_tet_outside
      << ", \"mesh_outside_tet_inside\": "
      << validation.mesh_outside_tet_inside
      << ", \"tet_inside_cluster_outside\": "
      << validation.tet_inside_cluster_outside
      << ", \"tet_outside_cluster_inside\": "
      << validation.tet_outside_cluster_inside
      << ", \"tetra_tree_mismatch\": "
      << validation.tetra_tree_mismatch
      << ", \"mesh_inside_tree_outside\": "
      << validation.mesh_inside_tree_outside
      << ", \"mesh_outside_tree_inside\": "
      << validation.mesh_outside_tree_inside
      << ", \"raw_simplified_mismatch\": "
      << validation.raw_simplified_mismatch << "},\n"
      << "  \"conservative_expansion\": {\"maximum_support_error\": "
      << metrics.maximum_support_expansion
      << ", \"maximum_observed_false_positive_distance\": "
      << metrics.maximum_observed_expansion << "},\n"
      << "  \"failures\": [";
  for(std::size_t i = 0; i < validation.failures.size(); ++i)
  {
    if(i != 0)
      out << ", ";
    out << '"' << validation.failures[i] << '"';
  }
  out << "]\n}\n";
}

void writeBenchmarkCsvHeader(std::ostream& out)
{
  out << "pipeline,length_unit,point_count,alpha,offset,wrap_faces,tetra_count,"
         "convex_parts,unique_planes,tree_nodes,tree_depth,alpha_wrap_ms,"
         "tet_extract_ms,agglomeration_ms,plane_reduction_ms,tree_build_ms,"
         "validation_ms,total_ms,false_positive,false_negative,"
         "preprocess_mode,smooth_neighbors,smooth_iterations,"
         "preprocessing_ms,smoothing_ms,max_smoothing_displacement,"
         "requested_offset,offset_retries,raw_points,occupied_voxels,voxel_points,"
         "voxel_point_ratio,voxel_size,parameter_mode,alpha_voxel_factor,offset_voxel_factor,"
         "convex_clusters,planes_before_reduce,planes_after_reduce,raw_tree_nodes,"
         "voxelization_ms,core_runtime_ms,total_runtime_ms,raw_points_outside_voxels,"
         "raw_points_outside_wrap,raw_points_outside_tree,containment_recovery_ms\n";
}

void writeBenchmarkCsvRow(std::ostream& out, const Metrics& metrics)
{
  out << std::setprecision(17) << metrics.pipeline << ',' << kLengthUnit
      << ',' << metrics.point_count << ',' << metrics.alpha << ',' << metrics.offset
      << ',' << metrics.wrap_faces << ',' << metrics.tetra_count << ','
      << metrics.final_clusters << ',' << metrics.unique_planes << ','
      << metrics.tree_nodes << ',' << metrics.tree_depth << ','
      << metrics.alpha_wrap_ms << ',' << metrics.tet_extract_ms << ','
      << metrics.agglomeration_ms << ',' << metrics.plane_reduction_ms << ','
      << metrics.tree_build_ms << ',' << metrics.validation_ms << ','
      << metrics.total_ms << ',' << metrics.false_positive << ','
      << metrics.false_negative << ','
      << pointCloudPreprocessModeName(metrics.preprocessing.mode) << ','
      << metrics.preprocessing.effective_neighbors << ','
      << metrics.preprocessing.iterations << ','
      << metrics.preprocessing.total_ms << ','
      << metrics.preprocessing.smoothing_ms << ','
      << metrics.preprocessing.maximum_displacement << ','
      << metrics.requested_offset << ',' << metrics.offset_retries << ','
      << metrics.point_count << ',' << metrics.preprocessing.voxelization.occupied_voxels << ','
      << metrics.preprocessing.voxelization.voxel_points << ','
      << metrics.preprocessing.voxelization.voxel_point_ratio << ','
      << metrics.preprocessing.voxelization.voxel_size << ',' << metrics.parameter_mode << ','
      << metrics.alpha_voxel_factor << ',' << metrics.offset_voxel_factor << ','
      << metrics.final_clusters << ',' << metrics.raw_cluster_planes << ','
      << metrics.reduced_cluster_planes << ',' << metrics.raw_tree_nodes << ','
      << metrics.preprocessing.voxelization.construction_ms << ','
      << metrics.core_runtime_ms << ',' << metrics.total_ms << ','
      << metrics.preprocessing.voxelization.raw_points_outside_voxels << ','
      << metrics.raw_points_outside_wrap << ',' << metrics.raw_points_outside_tree << ','
      << metrics.containment_recovery_ms << '\n';
}

void writeBenchmarkJsonRow(std::ostream& out, const Metrics& metrics)
{
  out << std::setprecision(17)
      << "{\"pipeline\": \"" << metrics.pipeline
      << "\", \"length_unit\": \"" << kLengthUnit
      << "\", \"point_count\": " << metrics.point_count
      << ", \"alpha\": " << metrics.alpha
      << ", \"requested_offset\": " << metrics.requested_offset
      << ", \"offset\": " << metrics.offset
      << ", \"offset_retries\": " << metrics.offset_retries
      << ", \"wrap_faces\": " << metrics.wrap_faces
      << ", \"finite_cells\": " << metrics.finite_cells
      << ", \"outside_cells\": " << metrics.outside_cells
      << ", \"tetra_count\": " << metrics.tetra_count
      << ", \"initial_tetra_clusters\": " << metrics.initial_clusters
      << ", \"convex_parts\": " << metrics.final_clusters
      << ", \"raw_cluster_planes\": " << metrics.raw_cluster_planes
      << ", \"reduced_cluster_planes\": "
      << metrics.reduced_cluster_planes << ", \"unique_planes\": "
      << metrics.unique_planes << ", \"tetra_tree_leaf_references\": "
      << metrics.tetra_tree_leaf_references
      << ", \"tree_leaf_references\": " << metrics.tree_leaf_references
      << ", \"raw_tree_nodes\": " << metrics.raw_tree_nodes
      << ", \"tree_nodes\": " << metrics.tree_nodes
      << ", \"tree_depth\": " << metrics.tree_depth
      << ", \"preprocess_mode\": \""
      << pointCloudPreprocessModeName(metrics.preprocessing.mode)
      << "\", \"smooth_neighbors\": "
      << metrics.preprocessing.effective_neighbors
      << ", \"smooth_iterations\": " << metrics.preprocessing.iterations
      << ", \"average_spacing\": "
      << metrics.preprocessing.average_spacing
      << ", \"maximum_smoothing_displacement\": "
      << metrics.preprocessing.maximum_displacement
      << ", \"scale_estimation_ms\": "
      << metrics.preprocessing.scale_estimation_ms
      << ", \"smoothing_ms\": " << metrics.preprocessing.smoothing_ms
      << ", \"preprocessing_ms\": " << metrics.preprocessing.total_ms
      << ", \"alpha_wrap_ms\": " << metrics.alpha_wrap_ms
      << ", \"tet_extract_ms\": " << metrics.tet_extract_ms
      << ", \"agglomeration_ms\": " << metrics.agglomeration_ms
      << ", \"plane_reduction_ms\": " << metrics.plane_reduction_ms
      << ", \"tree_build_ms\": " << metrics.tree_build_ms
      << ", \"simplification_ms\": " << metrics.simplification_ms
      << ", \"validation_ms\": " << metrics.validation_ms
      << ", \"total_ms\": " << metrics.total_ms
      << ", \"false_positive\": " << metrics.false_positive
      << ", \"false_negative\": " << metrics.false_negative
      << ", \"maximum_support_expansion\": "
      << metrics.maximum_support_expansion
      << ", \"maximum_observed_expansion\": "
      << metrics.maximum_observed_expansion
      << ", \"raw_points\": " << metrics.point_count
      << ", \"occupied_voxels\": " << metrics.preprocessing.voxelization.occupied_voxels
      << ", \"voxel_points\": " << metrics.preprocessing.voxelization.voxel_points
      << ", \"voxel_point_ratio\": " << metrics.preprocessing.voxelization.voxel_point_ratio
      << ", \"voxel_size\": " << metrics.preprocessing.voxelization.voxel_size
      << ", \"parameter_mode\": \"" << metrics.parameter_mode << '"'
      << ", \"alpha_voxel_factor\": " << metrics.alpha_voxel_factor
      << ", \"offset_voxel_factor\": " << metrics.offset_voxel_factor
      << ", \"convex_clusters\": " << metrics.final_clusters
      << ", \"planes_before_reduce\": " << metrics.raw_cluster_planes
      << ", \"planes_after_reduce\": " << metrics.reduced_cluster_planes
      << ", \"voxelization_ms\": " << metrics.preprocessing.voxelization.construction_ms
      << ", \"core_runtime_ms\": " << metrics.core_runtime_ms
      << ", \"total_runtime_ms\": " << metrics.total_ms
      << ", \"raw_points_outside_voxels\": "
      << metrics.preprocessing.voxelization.raw_points_outside_voxels
      << ", \"raw_points_outside_wrap\": " << metrics.raw_points_outside_wrap
      << ", \"raw_points_outside_tree\": " << metrics.raw_points_outside_tree
      << ", \"containment_recovery_ms\": " << metrics.containment_recovery_ms << '}';
}

void writeRunOutputs(const fs::path& output,
                     const Adapter& adapter,
                     const std::vector<HS::Vec3>& vertices,
                     const std::vector<ConvexCluster>& clusters,
                     const TreeRepresentation& tetra_tree,
                     const TreeRepresentation& cluster_tree,
                     const std::vector<Point>& wrap_input_points,
                     const Options& options,
                     const Metrics& metrics,
                     const Validation& validation)
{
  fs::create_directories(output / "convex_clusters");
  if(!CGAL::IO::write_polygon_mesh(
         (output / "wrap.off").string(), adapter.surfaceMesh(),
         CGAL::parameters::stream_precision(17)))
    throw std::runtime_error("Failed to write Alpha Wrap mesh");
  if(!CGAL::IO::write_points(
         (output / "alpha_wrap_input_points.xyz").string(),
         wrap_input_points, CGAL::parameters::stream_precision(17)))
    throw std::runtime_error("Failed to write preprocessed point cloud");
  for(const ConvexCluster& cluster : clusters)
  {
    std::ostringstream name;
    name << "cluster_" << std::setw(5) << std::setfill('0') << cluster.id
         << ".off";
    const Mesh mesh = clusterMesh(cluster, vertices);
    if(!CGAL::IO::write_polygon_mesh(
           (output / "convex_clusters" / name.str()).string(), mesh,
           CGAL::parameters::stream_precision(17)))
      throw std::runtime_error("Failed to write convex cluster mesh");
  }
  std::ofstream cluster_json(output / "convex_clusters.json");
  cluster_json << std::setprecision(17) << "{\n  \"length_unit\": \""
               << kLengthUnit << "\",\n  \"clusters\": [\n";
  for(std::size_t i = 0; i < clusters.size(); ++i)
  {
    const ConvexCluster& cluster = clusters[i];
    cluster_json << "    {\"id\": " << cluster.id
                 << ", \"cell_count\": " << cluster.cells.size()
                 << ", \"cells\": [";
    for(std::size_t j = 0; j < cluster.cells.size(); ++j)
    {
      if(j != 0)
        cluster_json << ", ";
      cluster_json << cluster.cells[j];
    }
    cluster_json << "], \"neighbors\": [";
    std::size_t neighbor_index = 0;
    for(const int neighbor : cluster.neighbors)
    {
      if(neighbor_index++ != 0)
        cluster_json << ", ";
      cluster_json << neighbor;
    }
    cluster_json << "], \"boundary_faces\": "
                 << cluster.boundary_faces.size()
                 << ", \"raw_planes\": " << cluster.exact_planes.size()
                 << ", \"reduced_planes\": "
                 << cluster.reduced_planes.size()
                 << ", \"centroid\": [" << cluster.centroid.x << ", "
                 << cluster.centroid.y << ", " << cluster.centroid.z
                 << "]}"
                 << (i + 1 == clusters.size() ? "\n" : ",\n");
  }
  cluster_json << "  ]\n}\n";
  writeTreeJson(output / "logic_tree_tetra_raw.json", tetra_tree,
                tetra_tree.raw, "tetra_raw");
  writeTreeJson(output / "logic_tree_raw.json", cluster_tree,
                cluster_tree.raw, "cluster_raw");
  writeTreeJson(output / "logic_tree_simplified.json", cluster_tree,
                cluster_tree.simplified, "simplified");
  std::ofstream raw_expression(output / "logic_expression_raw.txt");
  raw_expression << "inside(x) =\n"
                 << HS::booleanExpression(cluster_tree.raw, 0)
                 << "\n\nphi(x) =\n"
                 << HS::scalarExpression(cluster_tree.raw, 0) << '\n';
  std::ofstream simplified_expression(
      output / "logic_expression_simplified.txt");
  simplified_expression
      << "inside(x) =\n"
      << HS::booleanExpression(cluster_tree.simplified, 0)
      << "\n\nphi(x) =\n"
      << HS::scalarExpression(cluster_tree.simplified, 0) << '\n';
  writeValidationJson(output / "validation_report.json", options, metrics,
                      validation);

  std::ofstream csv(output / "benchmark.csv");
  writeBenchmarkCsvHeader(csv);
  writeBenchmarkCsvRow(csv, metrics);
  std::ofstream json(output / "benchmark.json");
  json << "[\n  ";
  writeBenchmarkJsonRow(json, metrics);
  json << "\n]\n";
}

RunResult runOnce(const Options& source_options,
                  const std::vector<Point>& input_points,
                  const PointCloudPreprocessResult& preprocessed,
                  const fs::path& output)
{
  const auto total_start = std::chrono::steady_clock::now();
  Options options = source_options;
  const std::array<double, 6> bbox = pointBoundingBox(input_points);
  const double diagonal = bboxDiagonal(bbox);
  if(!std::isfinite(diagonal) || diagonal <= 0.0)
    throw std::runtime_error("Point-set bounding box is degenerate");
  if(options.voxel_parameters)
  {
    options.alpha = options.alpha_voxel_factor * options.preprocessing.voxel_size;
    options.offset = options.offset_voxel_factor * options.preprocessing.voxel_size;
  }
  else if(!options.absolute_parameters)
  {
    options.alpha = diagonal / options.alpha_ratio;
    options.offset = diagonal / options.offset_ratio;
  }
  if(!std::isfinite(options.alpha) || !std::isfinite(options.offset) ||
     options.alpha <= 0.0 || options.offset <= 0.0)
    throw std::runtime_error("Resolved alpha/offset must be finite and positive");
  if(options.classification_epsilon == 0.0)
    options.classification_epsilon = std::max(1e-12, diagonal * 1e-10);
  if(options.support_epsilon == 0.0 &&
     options.plane_reduction == PlaneReductionMode::CONSERVATIVE)
    options.support_epsilon = diagonal * 1e-12;

  RunResult result;
  Metrics& metrics = result.metrics;
  metrics.point_count = input_points.size();
  metrics.alpha = options.alpha;
  metrics.requested_offset = options.offset;
  metrics.preprocessing = preprocessed.stats;
  metrics.parameter_mode = options.absolute_parameters ? "absolute" :
                           options.voxel_parameters ? "voxel-factor" : "bbox-ratio";
  if(options.voxel_parameters)
  {
    metrics.alpha_voxel_factor = options.alpha_voxel_factor;
    metrics.offset_voxel_factor = options.offset_voxel_factor;
  }
  if(options.auto_offset &&
     preprocessed.stats.mode == PointCloudPreprocessMode::JET)
  {
    options.offset = std::max(
        options.offset, options.offset_displacement_factor *
                            preprocessed.stats.maximum_displacement);
  }

  Adapter adapter;
  double accumulated_wrap_ms = 0.0;
  double accumulated_extract_ms = 0.0;
  while(true)
  {
    adapter.run(preprocessed.points, options.alpha, options.offset);
    Mesh& candidate_wrap = adapter.surfaceMesh();
    PMP::orient_to_bound_a_volume(
        candidate_wrap, CGAL::parameters::outward_orientation(true));
    accumulated_wrap_ms += adapter.statistics().alpha_wrap_ms;
    accumulated_extract_ms += adapter.statistics().tetra_extract_ms;

    if(!options.auto_offset ||
       preprocessed.stats.mode == PointCloudPreprocessMode::NONE)
      break;
    const auto containment_start = std::chrono::steady_clock::now();
    const std::size_t outside = countOutsidePoints(candidate_wrap, input_points);
    metrics.containment_recovery_ms += elapsedMilliseconds(containment_start);
    if(outside == 0 || metrics.offset_retries >= options.max_offset_retries)
      break;
    options.offset *= options.offset_growth_factor;
    if(!std::isfinite(options.offset))
      throw std::runtime_error("Offset recovery overflowed");
    ++metrics.offset_retries;
  }
  Mesh& wrap = adapter.surfaceMesh();
  const auto& adapter_statistics = adapter.statistics();
  metrics.offset = options.offset;
  metrics.alpha_wrap_ms = accumulated_wrap_ms;
  metrics.tet_extract_ms = accumulated_extract_ms;
  metrics.wrap_faces = wrap.number_of_faces();
  metrics.finite_cells = adapter_statistics.finite_cells;
  metrics.outside_cells = adapter_statistics.outside_cells;
  metrics.tetra_count = adapter.interiorCells().size();
  metrics.initial_clusters = metrics.tetra_count;

  std::vector<HS::Vec3> vertices;
  vertices.reserve(adapter.vertices().size());
  for(const Point& point : adapter.vertices())
    vertices.push_back(toVec(point));
  const std::array<double, 6> validation_bbox =
      pointBoundingBox(adapter.vertices());
  const double validation_diagonal =
      std::max(diagonal, bboxDiagonal(validation_bbox));

  const double plane_tolerance = std::max(1e-12, diagonal * 1e-12);
  TreeRepresentation tetra_tree = buildTetraTree(
      adapter.interiorCells(), vertices, plane_tolerance);
  metrics.tetra_tree_leaf_references = tetra_tree.leaf_references;

  auto start = std::chrono::steady_clock::now();
  ConvexAgglomerator agglomerator(vertices, adapter.interiorCells(), options,
                                  options.classification_epsilon);
  std::vector<ConvexCluster> clusters = agglomerator.run();
  metrics.agglomeration_ms = elapsedMilliseconds(start);
  metrics.final_clusters = clusters.size();

  start = std::chrono::steady_clock::now();
  SupportPlaneReducer reducer(vertices, options, options.support_epsilon);
  for(ConvexCluster& cluster : clusters)
  {
    metrics.raw_cluster_planes += cluster.boundary_faces.size();
    cluster.exact_planes = reducer.extractExactPlanes(cluster);
    cluster.reduced_planes =
        reducer.reduce(cluster, metrics.maximum_support_expansion);
    metrics.reduced_cluster_planes += cluster.reduced_planes.size();
  }
  metrics.plane_reduction_ms = elapsedMilliseconds(start);

  start = std::chrono::steady_clock::now();
  TreeRepresentation cluster_tree =
      buildClusterTree(clusters, plane_tolerance);
  metrics.tree_build_ms = elapsedMilliseconds(start);
  metrics.unique_planes = cluster_tree.leaves.size();
  metrics.tree_leaf_references = cluster_tree.leaf_references;
  metrics.raw_tree_nodes = HS::serializeTree(cluster_tree.raw).nodes.size();

  start = std::chrono::steady_clock::now();
  cluster_tree.simplified = HS::simplifyToFixedPoint(cluster_tree.raw);
  metrics.simplification_ms = elapsedMilliseconds(start);
  metrics.tree_nodes =
      HS::serializeTree(cluster_tree.simplified).nodes.size();
  metrics.tree_depth = HS::treeDepth(cluster_tree.simplified);

  metrics.core_runtime_ms = preprocessed.stats.total_ms + elapsedMilliseconds(total_start) -
      preprocessed.stats.voxelization.validation_ms - metrics.containment_recovery_ms;
  start = std::chrono::steady_clock::now();
  validateClusters(clusters, vertices, options.classification_epsilon,
                   result.validation);
  validateGeometry(input_points, preprocessed.points, adapter, clusters,
                   tetra_tree, cluster_tree, options, validation_bbox,
                   validation_diagonal, result.validation, metrics);
  metrics.validation_ms = elapsedMilliseconds(start) +
      preprocessed.stats.voxelization.validation_ms + metrics.containment_recovery_ms;
  metrics.total_ms =
      preprocessed.stats.total_ms + elapsedMilliseconds(total_start);
  result.passed = result.validation.failures.empty();

  writeRunOutputs(output, adapter, vertices, clusters, tetra_tree,
                  cluster_tree, preprocessed.points, options, metrics,
                  result.validation);

  std::cout << std::setprecision(17)
            << "CGAL version: " << kCgalVersion << "\n"
            << "Pipeline: alpha-tet\n"
            << "Length unit: " << kLengthUnit << "\n"
            << "Parameter mode: " << metrics.parameter_mode << "\n"
            << "raw_points: " << metrics.point_count << "\n"
            << "voxel_size (m): " << metrics.preprocessing.voxelization.voxel_size << "\n"
            << "occupied_voxels: " << metrics.preprocessing.voxelization.occupied_voxels << "\n"
            << "voxel_points: " << metrics.preprocessing.voxelization.voxel_points << "\n"
            << "voxel_point_ratio: " << metrics.preprocessing.voxelization.voxel_point_ratio << "\n"
            << "raw_points_outside_voxels: "
            << metrics.preprocessing.voxelization.raw_points_outside_voxels << "\n"
            << "raw_points_outside_tree: " << metrics.raw_points_outside_tree << "\n"
            << "Alpha Wrap API: internal, isolated by AlphaWrapVolumeAdapter\n"
            << "Plane reduction: "
            << (options.plane_reduction == PlaneReductionMode::EXACT
                    ? "exact"
                    : "conservative")
            << "\n"
            << "Preprocessing: "
            << pointCloudPreprocessModeName(metrics.preprocessing.mode)
            << "\n"
            << "Smooth neighbors: "
            << metrics.preprocessing.effective_neighbors << "\n"
            << "Smooth iterations: " << metrics.preprocessing.iterations
            << "\n"
            << "Average spacing (m): " << metrics.preprocessing.average_spacing
            << "\n"
            << "Maximum smoothing displacement (m): "
            << metrics.preprocessing.maximum_displacement << "\n"
            << "Raw input points outside: "
            << result.validation.input_points_outside << "\n"
            << "Wrap input points outside: "
            << result.validation.wrap_input_points_outside << "\n"
            << "alpha (m): " << metrics.alpha << "\n"
            << "requested offset (m): " << metrics.requested_offset << "\n"
            << "offset (m): " << metrics.offset << "\n"
            << "offset retries: " << metrics.offset_retries << "\n"
            << "wrap_surface_faces: " << metrics.wrap_faces << "\n"
            << "finite_cells: " << metrics.finite_cells << "\n"
            << "outside_cells: " << metrics.outside_cells << "\n"
            << "wrap_internal_tetrahedra: " << metrics.tetra_count << "\n"
            << "initial_tetra_clusters: " << metrics.initial_clusters << "\n"
            << "final_convex_clusters: " << metrics.final_clusters << "\n"
            << "raw_cluster_planes: " << metrics.raw_cluster_planes << "\n"
            << "reduced_cluster_planes: "
            << metrics.reduced_cluster_planes << "\n"
            << "reduced_unique_planes: " << metrics.unique_planes << "\n"
            << "tetra_tree_leaf_references: "
            << metrics.tetra_tree_leaf_references << "\n"
            << "tree_leaf_references: " << metrics.tree_leaf_references
            << "\n"
            << "tree_nodes: " << metrics.tree_nodes << "\n"
            << "tree_depth: " << metrics.tree_depth << "\n"
            << "mesh_inside_tree_outside: "
            << result.validation.mesh_inside_tree_outside << "\n"
            << "mesh_outside_tree_inside: "
            << result.validation.mesh_outside_tree_inside << "\n"
            << "alpha_wrap_ms: " << metrics.alpha_wrap_ms << "\n"
            << "tet_extract_ms: " << metrics.tet_extract_ms << "\n"
            << "agglomeration_ms: " << metrics.agglomeration_ms << "\n"
            << "plane_reduction_ms: " << metrics.plane_reduction_ms << "\n"
            << "tree_build_ms: " << metrics.tree_build_ms << "\n"
            << "scale_estimation_ms: "
            << metrics.preprocessing.scale_estimation_ms << "\n"
            << "smoothing_ms: " << metrics.preprocessing.smoothing_ms << "\n"
            << "validation_ms: " << metrics.validation_ms << "\n"
            << "voxelization_ms: " << metrics.preprocessing.voxelization.construction_ms << "\n"
            << "core_runtime_ms: " << metrics.core_runtime_ms << "\n"
            << "total_ms: " << metrics.total_ms << "\n"
            << "Output: " << output.string() << "\n";
  for(const std::string& failure : result.validation.failures)
    std::cerr << "Validation failure: " << failure << "\n";
  return result;
}

std::string ratioDirectoryName(double ratio)
{
  std::ostringstream stream;
  stream << std::setprecision(8) << ratio;
  std::string result = "alpha_ratio_" + stream.str();
  std::replace(result.begin(), result.end(), '.', '_');
  return result;
}

void writeAggregateBenchmarks(const fs::path& output,
                              const std::vector<RunResult>& runs)
{
  fs::create_directories(output);
  std::ofstream csv(output / "benchmark.csv");
  writeBenchmarkCsvHeader(csv);
  for(const RunResult& run : runs)
    writeBenchmarkCsvRow(csv, run.metrics);

  std::ofstream json(output / "benchmark.json");
  json << "[\n";
  for(std::size_t i = 0; i < runs.size(); ++i)
  {
    json << "  ";
    writeBenchmarkJsonRow(json, runs[i].metrics);
    json << (i + 1 == runs.size() ? "\n" : ",\n");
  }
  json << "]\n";
}

}  // namespace

int runAlphaTetPipeline(int argc, char** argv)
{
  Options options;
  if(!parseOptions(argc, argv, options))
  {
    printUsage(argv[0]);
    return EXIT_FAILURE;
  }

  try
  {
    std::vector<Point> points;
    if(!CGAL::IO::read_points(options.input_path,
                              std::back_inserter(points)) ||
       points.empty())
      throw std::runtime_error("Failed to read a non-empty point set");

    const fs::path output(options.output_directory);
    const PointCloudPreprocessResult preprocessed =
        preprocessPointCloud(points, options.preprocessing);
    if(preprocessed.stats.mode == PointCloudPreprocessMode::VOXEL &&
       preprocessed.points.size() > points.size())
      std::cerr << "Voxel corner input expanded by "
                << preprocessed.stats.voxelization.voxel_point_ratio << "x\n";
    if(options.benchmark_alpha_ratios.empty())
    {
      const RunResult result = runOnce(options, points, preprocessed, output);
      return result.passed ? EXIT_SUCCESS : EXIT_FAILURE;
    }

    std::vector<RunResult> runs;
    bool passed = true;
    for(const double ratio : options.benchmark_alpha_ratios)
    {
      Options run_options = options;
      run_options.alpha_ratio = ratio;
      run_options.benchmark_alpha_ratios.clear();
      RunResult result = runOnce(run_options, points, preprocessed,
                                 output / ratioDirectoryName(ratio));
      passed &= result.passed;
      runs.push_back(std::move(result));
    }
    writeAggregateBenchmarks(output, runs);
    return passed ? EXIT_SUCCESS : EXIT_FAILURE;
  }
  catch(const std::exception& error)
  {
    std::cerr << "alpha-tet pipeline failed: " << error.what() << "\n";
    return EXIT_FAILURE;
  }
}

}  // namespace rokae_demo
