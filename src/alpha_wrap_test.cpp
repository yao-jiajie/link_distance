#include <rokae_demo/alpha_wrap_volume_adapter.hpp>
#include <rokae_demo/point_cloud_preprocessing.hpp>
#include <rokae_demo/units.hpp>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/IO/write_points.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/alpha_wrap_3.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/version.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace PMP = CGAL::Polygon_mesh_processing;

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point3 = Kernel::Point_3;
using Vector3 = Kernel::Vector_3;
using SurfaceMesh = CGAL::Surface_mesh<Point3>;
using AlphaWrapVolumeAdapter = rokae_demo::AlphaWrapVolumeAdapter;
using Primitive = CGAL::AABB_face_graph_triangle_primitive<SurfaceMesh>;
using AabbTraits = CGAL::AABB_traits<Kernel, Primitive>;
using AabbTree = CGAL::AABB_tree<AabbTraits>;

#define ROKAE_STRINGIFY_IMPL(value) #value
#define ROKAE_STRINGIFY(value) ROKAE_STRINGIFY_IMPL(value)
static constexpr const char* kCgalVersion = ROKAE_STRINGIFY(CGAL_VERSION);

struct Options
{
  struct Probe
  {
    Point3 point;
    bool expected_inside = false;
  };

  std::string input_path;
  std::string output_off_path;
  double alpha = 0.0;
  double offset = 0.0;
  double alpha_ratio = 20.0;
  double offset_ratio = 100.0;
  bool use_absolute = false;
  bool voxel_parameters = false;
  double alpha_voxel_factor = 3.0;
  double offset_voxel_factor = 1.0;
  std::size_t volume_validation_samples = 10000;
  std::vector<Probe> probes;
  rokae_demo::PointCloudPreprocessOptions preprocessing;
  bool auto_offset = true;
  double offset_displacement_factor = 1.75;
  double offset_growth_factor = 1.5;
  std::size_t max_offset_retries = 4;
};

struct Validation
{
  bool triangle_mesh = false;
  bool closed = false;
  bool self_intersects = true;
  bool outward_oriented = false;
  bool bounds_volume = false;
  std::size_t bounded = 0;
  std::size_t boundary = 0;
  std::size_t unbounded = 0;
  std::size_t wrap_input_bounded = 0;
  std::size_t wrap_input_boundary = 0;
  std::size_t wrap_input_unbounded = 0;
};

struct VolumeValidation
{
  std::size_t tested = 0;
  std::size_t boundary_ignored = 0;
  std::size_t mesh_inside_tet_outside = 0;
  std::size_t mesh_outside_tet_inside = 0;
};

static void printUsage(const char* program)
{
  std::cerr
      << "Usage: " << program << " INPUT.{xyz,ply,off} OUTPUT.off [parameters]\n"
      << "Input coordinates and all length parameters use meters (m).\n"
      << "Absolute mode: --alpha VALUE_METERS --offset VALUE_METERS\n"
      << "Relative mode: --alpha-ratio VALUE --offset-ratio VALUE\n"
      << "Voxel mode: --voxel-size METERS [--alpha-voxel-factor 3 --offset-voxel-factor 1]\n"
      << "Preprocessing: --preprocess none|jet|voxel "
         "--smooth-neighbors auto|COUNT --smooth-iterations COUNT\n"
      << "Offset recovery: --auto-offset true|false "
         "--offset-displacement-factor VALUE\n"
      << "  --offset-growth-factor VALUE --max-offset-retries COUNT\n"
      << "Validation: --volume-validation-samples COUNT\n";
}

static bool parseNonnegativeInteger(const std::string& text,
                                    std::size_t& value)
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

static bool parsePositiveDouble(const std::string& text, double& value)
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

static bool parseFiniteDouble(const std::string& text, double& value)
{
  try
  {
    std::size_t parsed = 0;
    value = std::stod(text, &parsed);
    return parsed == text.size() && std::isfinite(value);
  }
  catch(const std::exception&)
  {
    return false;
  }
}

static bool parseOptions(int argc, char** argv, Options& options)
{
  if(argc < 3)
  {
    printUsage(argv[0]);
    return false;
  }

  options.input_path = argv[1];
  options.output_off_path = argv[2];
  bool has_alpha = false;
  bool has_offset = false;
  bool has_alpha_ratio = false;
  bool has_offset_ratio = false;
  bool has_factors = false;
  bool has_preprocess = false;

  for(int i = 3; i < argc; ++i)
  {
    const std::string option = argv[i];
    if(option == "--probe")
    {
      if(i + 4 >= argc)
      {
        std::cerr << "--probe requires X Y Z inside|outside\n";
        return false;
      }
      double x = 0.0;
      double y = 0.0;
      double z = 0.0;
      if(!parseFiniteDouble(argv[++i], x) ||
         !parseFiniteDouble(argv[++i], y) ||
         !parseFiniteDouble(argv[++i], z))
      {
        std::cerr << "--probe coordinates must be finite numbers\n";
        return false;
      }
      const std::string expected = argv[++i];
      if(expected != "inside" && expected != "outside")
      {
        std::cerr << "--probe expectation must be inside or outside\n";
        return false;
      }
      options.probes.push_back(
          {Point3(x, y, z), expected == "inside"});
      continue;
    }
    if(option == "--volume-validation-samples")
    {
      if(i + 1 >= argc ||
         !parseNonnegativeInteger(argv[++i],
                                  options.volume_validation_samples))
      {
        std::cerr << "--volume-validation-samples requires an integer\n";
        return false;
      }
      continue;
    }
    if(option == "--preprocess")
    {
      has_preprocess = true;
      if(i + 1 >= argc)
        return false;
      const std::string mode = argv[++i];
      if(mode == "none")
        options.preprocessing.mode =
            rokae_demo::PointCloudPreprocessMode::NONE;
      else if(mode == "jet")
        options.preprocessing.mode =
            rokae_demo::PointCloudPreprocessMode::JET;
      else if(mode == "voxel")
        options.preprocessing.mode =
            rokae_demo::PointCloudPreprocessMode::VOXEL;
      else
      {
        std::cerr << "--preprocess requires none, jet, or voxel\n";
        return false;
      }
      continue;
    }
    if(option == "--smooth-neighbors")
    {
      if(i + 1 >= argc)
        return false;
      const std::string neighbors = argv[++i];
      if(neighbors == "auto")
        options.preprocessing.neighbors = 0;
      else if(!parseNonnegativeInteger(neighbors,
                                      options.preprocessing.neighbors) ||
              options.preprocessing.neighbors < 6)
      {
        std::cerr << "--smooth-neighbors requires auto or an integer >= 6\n";
        return false;
      }
      continue;
    }
    if(option == "--smooth-iterations")
    {
      if(i + 1 >= argc ||
         !parseNonnegativeInteger(argv[++i],
                                  options.preprocessing.iterations) ||
         options.preprocessing.iterations == 0)
      {
        std::cerr << "--smooth-iterations requires a positive integer\n";
        return false;
      }
      continue;
    }
    if(option == "--auto-offset")
    {
      if(i + 1 >= argc)
        return false;
      const std::string enabled = argv[++i];
      if(enabled == "true")
        options.auto_offset = true;
      else if(enabled == "false")
        options.auto_offset = false;
      else
      {
        std::cerr << "--auto-offset requires true or false\n";
        return false;
      }
      continue;
    }
    if(option == "--max-offset-retries")
    {
      if(i + 1 >= argc ||
         !parseNonnegativeInteger(argv[++i], options.max_offset_retries))
      {
        std::cerr << "--max-offset-retries requires an integer\n";
        return false;
      }
      continue;
    }
    if(i + 1 >= argc)
    {
      std::cerr << "Missing value for " << option << "\n";
      return false;
    }

    double value = 0.0;
    if(!parsePositiveDouble(argv[++i], value))
    {
      std::cerr << "Expected a finite positive value after " << option << "\n";
      return false;
    }

    if(option == "--alpha")
    {
      options.alpha = value;
      has_alpha = true;
    }
    else if(option == "--offset")
    {
      options.offset = value;
      has_offset = true;
    }
    else if(option == "--alpha-ratio")
    {
      options.alpha_ratio = value;
      has_alpha_ratio = true;
    }
    else if(option == "--offset-ratio")
    {
      options.offset_ratio = value;
      has_offset_ratio = true;
    }
    else if(option == "--voxel-size")
    {
      options.preprocessing.voxel_size = value;
    }
    else if(option == "--alpha-voxel-factor")
    {
      options.alpha_voxel_factor = value;
      has_factors = true;
    }
    else if(option == "--offset-voxel-factor")
    {
      options.offset_voxel_factor = value;
      has_factors = true;
    }
    else if(option == "--offset-displacement-factor")
    {
      options.offset_displacement_factor = value;
    }
    else if(option == "--offset-growth-factor")
    {
      if(value <= 1.0)
      {
        std::cerr << "--offset-growth-factor must be greater than 1\n";
        return false;
      }
      options.offset_growth_factor = value;
    }
    else
    {
      std::cerr << "Unknown option: " << option << "\n";
      return false;
    }
  }

  const bool any_absolute = has_alpha || has_offset;
  const bool any_relative = has_alpha_ratio || has_offset_ratio;
  if((any_absolute && any_relative) || (has_factors && (any_absolute || any_relative)))
  {
    std::cerr << "Absolute, voxel-factor, and bbox-ratio modes cannot be mixed.\n";
    return false;
  }
  if(any_absolute && !(has_alpha && has_offset))
  {
    std::cerr << "Absolute mode requires both --alpha and --offset.\n";
    return false;
  }

  if(any_relative && !(has_alpha_ratio && has_offset_ratio))
  {
    std::cerr << "BBox mode requires both --alpha-ratio and --offset-ratio.\n";
    return false;
  }
  const bool has_voxel = options.preprocessing.voxel_size > 0.0;
  if(has_voxel)
  {
    if(has_preprocess && options.preprocessing.mode != rokae_demo::PointCloudPreprocessMode::VOXEL)
    {
      std::cerr << "--voxel-size cannot be combined with --preprocess none/jet.\n";
      return false;
    }
    options.preprocessing.mode = rokae_demo::PointCloudPreprocessMode::VOXEL;
  }
  if(!has_voxel && (has_factors || options.preprocessing.mode == rokae_demo::PointCloudPreprocessMode::VOXEL))
  {
    std::cerr << "Voxel mode requires --voxel-size in meters.\n";
    return false;
  }
  options.use_absolute = any_absolute;
  options.voxel_parameters = has_voxel && !any_absolute && !any_relative;
  return true;
}

static bool validateProbes(const SurfaceMesh& mesh,
                           const std::vector<Options::Probe>& probes)
{
  if(probes.empty())
    return true;
  CGAL::Side_of_triangle_mesh<SurfaceMesh, Kernel> side(mesh);
  bool passed = true;
  for(std::size_t index = 0; index < probes.size(); ++index)
  {
    const CGAL::Bounded_side result = side(probes[index].point);
    const bool inside =
        result == CGAL::ON_BOUNDED_SIDE || result == CGAL::ON_BOUNDARY;
    const bool matches = inside == probes[index].expected_inside;
    std::cout << "Probe[" << index << "]: "
              << (result == CGAL::ON_BOUNDED_SIDE
                      ? "inside"
                      : (result == CGAL::ON_BOUNDARY ? "boundary"
                                                     : "outside"))
              << " expected="
              << (probes[index].expected_inside ? "inside" : "outside")
              << " pass=" << (matches ? "yes" : "no") << "\n";
    passed &= matches;
  }
  return passed;
}

static std::string pathWithoutExtension(const std::string& path)
{
  const std::size_t slash = path.find_last_of("/\\");
  const std::size_t dot = path.find_last_of('.');
  if(dot == std::string::npos || (slash != std::string::npos && dot < slash))
    return path;
  return path.substr(0, dot);
}

// CGAL 5.5+ only requires alpha and offset to be strictly positive; it does
// not require offset < alpha. The latter is usually a useful operating regime,
// so main() emits a warning when it is not satisfied.
bool generateAlphaWrap(const std::vector<Point3>& points,
                       double alpha,
                       double offset,
                       AlphaWrapVolumeAdapter& adapter)
{
  if(points.empty() || !std::isfinite(alpha) || !std::isfinite(offset) ||
     alpha <= 0.0 || offset <= 0.0)
    return false;

  adapter.run(points, alpha, offset);
  SurfaceMesh& output_mesh = adapter.surfaceMesh();

  if(output_mesh.is_empty() || !CGAL::is_triangle_mesh(output_mesh) ||
     !CGAL::is_closed(output_mesh))
    return false;

  if(PMP::does_self_intersect(output_mesh))
    return false;

  // Normalize every connected component to an outward orientation before
  // normals and half-space leaves are exported.
  PMP::orient_to_bound_a_volume(
      output_mesh, CGAL::parameters::outward_orientation(true));

  return PMP::does_bound_a_volume(output_mesh) &&
         PMP::is_outward_oriented(output_mesh);
}

static VolumeValidation validateTetrahedralRepresentation(
    const AlphaWrapVolumeAdapter& adapter,
    std::size_t sample_count,
    const std::array<double, 6>& bbox,
    double diagonal,
    double boundary_epsilon)
{
  VolumeValidation result;
  const SurfaceMesh& mesh = adapter.surfaceMesh();
  CGAL::Side_of_triangle_mesh<SurfaceMesh, Kernel> side(mesh);
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
  const double epsilon_squared = boundary_epsilon * boundary_epsilon;
  for(std::size_t i = 0; i < sample_count; ++i)
  {
    const Point3 point(random_x(generator), random_y(generator),
                       random_z(generator));
    const CGAL::Bounded_side mesh_side = side(point);
    if(mesh_side == CGAL::ON_BOUNDARY ||
       distance_tree.squared_distance(point) <= epsilon_squared)
    {
      ++result.boundary_ignored;
      continue;
    }
    const bool mesh_inside = mesh_side == CGAL::ON_BOUNDED_SIDE;
    const bool tetra_inside = adapter.containsInRetainedCells(point);
    if(mesh_inside && !tetra_inside)
      ++result.mesh_inside_tet_outside;
    else if(!mesh_inside && tetra_inside)
      ++result.mesh_outside_tet_inside;
    ++result.tested;
  }
  return result;
}

static Validation validateMeshAndContainment(const SurfaceMesh& mesh,
                                             const std::vector<Point3>& points,
                                             const std::vector<Point3>&
                                                 wrap_input_points)
{
  Validation validation;
  validation.triangle_mesh = CGAL::is_triangle_mesh(mesh);
  validation.closed = CGAL::is_closed(mesh);
  if(!mesh.is_empty() && validation.triangle_mesh)
    validation.self_intersects = PMP::does_self_intersect(mesh);

  if(validation.closed && validation.triangle_mesh &&
     !validation.self_intersects)
  {
    validation.bounds_volume = PMP::does_bound_a_volume(mesh);
    validation.outward_oriented = PMP::is_outward_oriented(mesh);

    // This classifies only the supplied discrete samples. Passing this check
    // does not prove containment of an unknown surface between those samples.
    CGAL::Side_of_triangle_mesh<SurfaceMesh, Kernel> side(mesh);
    std::size_t printed_outside = 0;
    for(std::size_t i = 0; i < points.size(); ++i)
    {
      const CGAL::Bounded_side result = side(points[i]);
      if(result == CGAL::ON_BOUNDED_SIDE)
      {
        ++validation.bounded;
      }
      else if(result == CGAL::ON_BOUNDARY)
      {
        ++validation.boundary;
      }
      else
      {
        ++validation.unbounded;
        if(printed_outside < 10)
        {
          std::cerr << "Outside point[" << i << "] = ["
                    << CGAL::to_double(points[i].x()) << ", "
                    << CGAL::to_double(points[i].y()) << ", "
                    << CGAL::to_double(points[i].z()) << "]\n";
          ++printed_outside;
        }
      }
    }
    for(const Point3& point : wrap_input_points)
    {
      const CGAL::Bounded_side result = side(point);
      if(result == CGAL::ON_BOUNDED_SIDE)
        ++validation.wrap_input_bounded;
      else if(result == CGAL::ON_BOUNDARY)
        ++validation.wrap_input_boundary;
      else
        ++validation.wrap_input_unbounded;
    }
  }
  return validation;
}

static std::size_t countOutsidePoints(const SurfaceMesh& mesh,
                                      const std::vector<Point3>& points)
{
  CGAL::Side_of_triangle_mesh<SurfaceMesh, Kernel> side(mesh);
  return static_cast<std::size_t>(std::count_if(
      points.begin(), points.end(), [&](const Point3& point) {
        return side(point) == CGAL::ON_UNBOUNDED_SIDE;
      }));
}

static std::vector<std::size_t> faceVertexIds(
    const SurfaceMesh& mesh,
    SurfaceMesh::Face_index face,
    const std::map<SurfaceMesh::Vertex_index, std::size_t>& vertex_ids)
{
  std::vector<std::size_t> result;
  for(const SurfaceMesh::Vertex_index vertex :
      CGAL::vertices_around_face(mesh.halfedge(face), mesh))
    result.push_back(vertex_ids.at(vertex));
  return result;
}

static std::vector<std::size_t> adjacentFaceIds(
    const SurfaceMesh& mesh,
    SurfaceMesh::Face_index face,
    const std::map<SurfaceMesh::Face_index, std::size_t>& face_ids)
{
  std::set<std::size_t> neighbors;
  for(const SurfaceMesh::Halfedge_index halfedge :
      CGAL::halfedges_around_face(mesh.halfedge(face), mesh))
  {
    const SurfaceMesh::Face_index neighbor =
        mesh.face(mesh.opposite(halfedge));
    if(neighbor != SurfaceMesh::null_face())
      neighbors.insert(face_ids.at(neighbor));
  }
  return std::vector<std::size_t>(neighbors.begin(), neighbors.end());
}

static bool writeGeometryJson(const std::string& path,
                              const SurfaceMesh& mesh,
                              double alpha,
                              double requested_offset,
                              double offset,
                              std::size_t offset_retries,
                              const rokae_demo::PointCloudPreprocessStats&
                                  preprocessing,
                              const Validation& validation)
{
  std::ofstream out(path);
  if(!out)
    return false;
  out << std::setprecision(17);

  std::map<SurfaceMesh::Vertex_index, std::size_t> vertex_ids;
  std::map<SurfaceMesh::Face_index, std::size_t> face_ids;
  std::size_t id = 0;
  for(const SurfaceMesh::Vertex_index vertex : mesh.vertices())
    vertex_ids[vertex] = id++;
  id = 0;
  for(const SurfaceMesh::Face_index face : mesh.faces())
    face_ids[face] = id++;

  out << "{\n"
      << "  \"metadata\": {\n"
      << "    \"units\": {\"length\": \"" << rokae_demo::kLengthUnit
      << "\", \"area\": \"" << rokae_demo::kAreaUnit
      << "\", \"volume\": \"" << rokae_demo::kVolumeUnit << "\"},\n"
      << "    \"alpha\": " << alpha << ",\n"
      << "    \"requested_offset\": " << requested_offset << ",\n"
      << "    \"offset\": " << offset << ",\n"
      << "    \"offset_retries\": " << offset_retries << ",\n"
      << "    \"raw_points\": " << preprocessing.input_points << ",\n"
      << "    \"voxel_size\": " << preprocessing.voxelization.voxel_size << ",\n"
      << "    \"occupied_voxels\": " << preprocessing.voxelization.occupied_voxels << ",\n"
      << "    \"voxel_points\": " << preprocessing.voxelization.voxel_points << ",\n"
      << "    \"voxel_point_ratio\": " << preprocessing.voxelization.voxel_point_ratio << ",\n"
      << "    \"raw_points_outside_voxels\": "
      << preprocessing.voxelization.raw_points_outside_voxels << ",\n"
      << "    \"raw_points_outside_wrap\": " << validation.unbounded << ",\n"
      << "    \"preprocessing\": {\"mode\": \""
      << rokae_demo::pointCloudPreprocessModeName(preprocessing.mode)
      << "\", \"input_points\": " << preprocessing.input_points
      << ", \"output_points\": " << preprocessing.output_points
      << ", \"requested_neighbors\": "
      << preprocessing.requested_neighbors
      << ", \"effective_neighbors\": "
      << preprocessing.effective_neighbors
      << ", \"iterations\": " << preprocessing.iterations
      << ", \"average_spacing\": " << preprocessing.average_spacing
      << ", \"mean_displacement\": "
      << preprocessing.mean_displacement
      << ", \"rms_displacement\": " << preprocessing.rms_displacement
      << ", \"percentile95_displacement\": "
      << preprocessing.percentile95_displacement
      << ", \"maximum_displacement\": "
      << preprocessing.maximum_displacement
      << ", \"scale_estimation_ms\": "
      << preprocessing.scale_estimation_ms
      << ", \"smoothing_ms\": " << preprocessing.smoothing_ms
      << ", \"total_ms\": " << preprocessing.total_ms << "},\n"
      << "    \"halfspace_convention\": \"psi_i(x) = n_i^T (x - w_i) = n_i^T x + d_i\",\n"
      << "    \"validation\": {\"closed\": "
      << (validation.closed ? "true" : "false")
      << ", \"triangle_mesh\": "
      << (validation.triangle_mesh ? "true" : "false")
      << ", \"self_intersects\": "
      << (validation.self_intersects ? "true" : "false")
      << ", \"outward_oriented\": "
      << (validation.outward_oriented ? "true" : "false")
      << ", \"bounded_points\": " << validation.bounded
      << ", \"boundary_points\": " << validation.boundary
      << ", \"unbounded_points\": " << validation.unbounded << "},\n"
      << "    \"wrap_input_validation\": {\"bounded_points\": "
      << validation.wrap_input_bounded << ", \"boundary_points\": "
      << validation.wrap_input_boundary << ", \"unbounded_points\": "
      << validation.wrap_input_unbounded << "}\n"
      << "  },\n";

  out << "  \"vertices\": [\n";
  std::size_t vertex_counter = 0;
  for(const SurfaceMesh::Vertex_index vertex : mesh.vertices())
  {
    const Point3& point = mesh.point(vertex);
    out << "    {\"id\": " << vertex_ids.at(vertex)
        << ", \"position\": [" << CGAL::to_double(point.x()) << ", "
        << CGAL::to_double(point.y()) << ", "
        << CGAL::to_double(point.z()) << "]}";
    out << (++vertex_counter == mesh.number_of_vertices() ? "\n" : ",\n");
  }
  out << "  ],\n";

  out << "  \"edges\": [\n";
  std::size_t edge_counter = 0;
  for(const SurfaceMesh::Edge_index edge : mesh.edges())
  {
    const SurfaceMesh::Halfedge_index halfedge = mesh.halfedge(edge);
    const SurfaceMesh::Vertex_index source = mesh.source(halfedge);
    const SurfaceMesh::Vertex_index target = mesh.target(halfedge);
    const SurfaceMesh::Face_index face0 = mesh.face(halfedge);
    const SurfaceMesh::Face_index face1 = mesh.face(mesh.opposite(halfedge));
    const double length =
        std::sqrt(CGAL::to_double(CGAL::squared_distance(
            mesh.point(source), mesh.point(target))));

    out << "    {\"id\": " << edge_counter
        << ", \"vertices\": [" << vertex_ids.at(source) << ", "
        << vertex_ids.at(target) << "], \"adjacent_faces\": [";
    if(face0 != SurfaceMesh::null_face())
      out << face_ids.at(face0);
    else
      out << -1;
    out << ", ";
    if(face1 != SurfaceMesh::null_face())
      out << face_ids.at(face1);
    else
      out << -1;
    out << "], \"length\": " << length << "}";
    out << (++edge_counter == mesh.number_of_edges() ? "\n" : ",\n");
  }
  out << "  ],\n";

  out << "  \"faces\": [\n";
  std::size_t face_counter = 0;
  bool normals_are_unit = true;
  for(const SurfaceMesh::Face_index face : mesh.faces())
  {
    const std::vector<std::size_t> vertices =
        faceVertexIds(mesh, face, vertex_ids);
    if(vertices.size() != 3)
      return false;

    std::vector<Point3> points;
    for(const SurfaceMesh::Vertex_index vertex :
        CGAL::vertices_around_face(mesh.halfedge(face), mesh))
      points.push_back(mesh.point(vertex));

    const Vector3 cross =
        CGAL::cross_product(points[1] - points[0], points[2] - points[0]);
    const double norm = std::sqrt(CGAL::to_double(cross.squared_length()));
    if(!std::isfinite(norm) || norm <= std::numeric_limits<double>::epsilon())
      return false;
    const double nx = CGAL::to_double(cross.x()) / norm;
    const double ny = CGAL::to_double(cross.y()) / norm;
    const double nz = CGAL::to_double(cross.z()) / norm;
    const double unit_norm = std::sqrt(nx * nx + ny * ny + nz * nz);
    if(std::abs(unit_norm - 1.0) > 1e-9)
      normals_are_unit = false;

    const double wx = CGAL::to_double(points[0].x());
    const double wy = CGAL::to_double(points[0].y());
    const double wz = CGAL::to_double(points[0].z());
    const double plane_offset = -(nx * wx + ny * wy + nz * wz);
    const std::vector<std::size_t> neighbors =
        adjacentFaceIds(mesh, face, face_ids);

    out << "    {\"id\": " << face_ids.at(face)
        << ", \"vertices\": [" << vertices[0] << ", " << vertices[1]
        << ", " << vertices[2] << "], \"normal\": ["
        << nx << ", " << ny << ", " << nz << "], \"point\": ["
        << wx << ", " << wy << ", " << wz << "], \"offset\": "
        << plane_offset << ", \"adjacent_faces\": [";
    for(std::size_t i = 0; i < neighbors.size(); ++i)
    {
      if(i != 0)
        out << ", ";
      out << neighbors[i];
    }
    out << "]}";
    out << (++face_counter == mesh.number_of_faces() ? "\n" : ",\n");
  }
  out << "  ]\n}\n";
  return static_cast<bool>(out) && normals_are_unit;
}

int main(int argc, char** argv)
{
  Options options;
  if(!parseOptions(argc, argv, options))
    return EXIT_FAILURE;

  std::vector<Point3> points;
  if(!CGAL::IO::read_points(options.input_path, std::back_inserter(points)) ||
     points.empty())
  {
    std::cerr << "Failed to read a non-empty XYZ/PLY/OFF point set: "
              << options.input_path << "\n";
    return EXIT_FAILURE;
  }
  rokae_demo::PointCloudPreprocessResult preprocessed;
  try
  {
    preprocessed =
        rokae_demo::preprocessPointCloud(points, options.preprocessing);
  }
  catch(const std::exception& error)
  {
    std::cerr << "Point-cloud preprocessing failed: " << error.what()
              << "\n";
    return EXIT_FAILURE;
  }

  double xmin = std::numeric_limits<double>::infinity();
  double ymin = std::numeric_limits<double>::infinity();
  double zmin = std::numeric_limits<double>::infinity();
  double xmax = -std::numeric_limits<double>::infinity();
  double ymax = -std::numeric_limits<double>::infinity();
  double zmax = -std::numeric_limits<double>::infinity();
  for(const Point3& point : points)
  {
    const double x = CGAL::to_double(point.x());
    const double y = CGAL::to_double(point.y());
    const double z = CGAL::to_double(point.z());
    xmin = std::min(xmin, x);
    ymin = std::min(ymin, y);
    zmin = std::min(zmin, z);
    xmax = std::max(xmax, x);
    ymax = std::max(ymax, y);
    zmax = std::max(zmax, z);
  }
  const double dx = xmax - xmin;
  const double dy = ymax - ymin;
  const double dz = zmax - zmin;
  const double bbox_diagonal = std::sqrt(dx * dx + dy * dy + dz * dz);
  if(!std::isfinite(bbox_diagonal) || bbox_diagonal <= 0.0)
  {
    std::cerr << "The point-set bounding box must have a positive diagonal.\n";
    return EXIT_FAILURE;
  }

  if(options.voxel_parameters)
  {
    options.alpha = options.alpha_voxel_factor * options.preprocessing.voxel_size;
    options.offset = options.offset_voxel_factor * options.preprocessing.voxel_size;
  }
  else if(!options.use_absolute)
  {
    options.alpha = bbox_diagonal / options.alpha_ratio;
    options.offset = bbox_diagonal / options.offset_ratio;
  }
  const double requested_offset = options.offset;
  if(!std::isfinite(options.alpha) || !std::isfinite(options.offset) ||
     options.alpha <= 0.0 || options.offset <= 0.0)
  {
    std::cerr << "Resolved alpha/offset must be finite and positive\n";
    return EXIT_FAILURE;
  }
  if(options.auto_offset &&
     preprocessed.stats.mode == rokae_demo::PointCloudPreprocessMode::JET)
  {
    options.offset = std::max(
        options.offset, options.offset_displacement_factor *
                            preprocessed.stats.maximum_displacement);
  }
  if(options.offset >= options.alpha)
  {
    std::cerr << "Warning: offset >= alpha. CGAL permits this, but a small "
                 "fraction of alpha is normally recommended.\n";
  }

  std::cout << std::setprecision(17)
            << "长度单位: " << rokae_demo::kLengthUnit << "\n"
            << "点云数量: " << points.size() << "\n"
            << "包围盒范围: [" << xmin << ", " << ymin << ", " << zmin
            << "] - [" << xmax << ", " << ymax << ", " << zmax << "]\n"
            << "包围盒对角线(m): " << bbox_diagonal << "\n"
            << "alpha(m): " << options.alpha << "\n"
            << "requested offset(m): " << requested_offset << "\n"
            << "offset(m): " << options.offset << "\n";

  AlphaWrapVolumeAdapter adapter;
  const auto start = std::chrono::steady_clock::now();
  bool generated = false;
  std::size_t offset_retries = 0;
  while(true)
  {
    generated = generateAlphaWrap(
        preprocessed.points, options.alpha, options.offset, adapter);
    if(!generated || !options.auto_offset ||
       preprocessed.stats.mode == rokae_demo::PointCloudPreprocessMode::NONE ||
       countOutsidePoints(adapter.surfaceMesh(), points) == 0 ||
       offset_retries >= options.max_offset_retries)
      break;
    options.offset *= options.offset_growth_factor;
    if(!std::isfinite(options.offset))
    {
      std::cerr << "Offset recovery overflowed\n";
      return EXIT_FAILURE;
    }
    ++offset_retries;
  }
  const double runtime_ms =
      std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - start)
          .count();
  if(!generated)
  {
    std::cerr << "Alpha wrapping failed output validity checks.\n";
    return EXIT_FAILURE;
  }
  const SurfaceMesh& mesh = adapter.surfaceMesh();

  const Validation validation =
      validateMeshAndContainment(mesh, points, preprocessed.points);
  std::array<double, 6> bbox{
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity()};
  for(const Point3& point : adapter.vertices())
  {
    const double x = CGAL::to_double(point.x());
    const double y = CGAL::to_double(point.y());
    const double z = CGAL::to_double(point.z());
    bbox[0] = std::min(bbox[0], x);
    bbox[1] = std::min(bbox[1], y);
    bbox[2] = std::min(bbox[2], z);
    bbox[3] = std::max(bbox[3], x);
    bbox[4] = std::max(bbox[4], y);
    bbox[5] = std::max(bbox[5], z);
  }
  const double wrap_dx = bbox[3] - bbox[0];
  const double wrap_dy = bbox[4] - bbox[1];
  const double wrap_dz = bbox[5] - bbox[2];
  const double validation_diagonal = std::max(
      bbox_diagonal,
      std::sqrt(wrap_dx * wrap_dx + wrap_dy * wrap_dy +
                wrap_dz * wrap_dz));
  const VolumeValidation volume_validation =
      validateTetrahedralRepresentation(
          adapter, options.volume_validation_samples, bbox,
          validation_diagonal,
          std::max(1e-12, validation_diagonal * 1e-10));
  const bool probes_passed = validateProbes(mesh, options.probes);
  const std::string output_base = pathWithoutExtension(options.output_off_path);
  const std::string output_ply_path = output_base + ".ply";
  const std::string output_json_path = output_base + "_geometry.json";
  const std::string smoothed_points_path = output_base + "_smoothed.xyz";

  const bool wrote_off = CGAL::IO::write_polygon_mesh(
      options.output_off_path, mesh, CGAL::parameters::stream_precision(17));
  const bool wrote_ply = CGAL::IO::write_polygon_mesh(
      output_ply_path, mesh, CGAL::parameters::stream_precision(17));
  const bool wrote_json = writeGeometryJson(
      output_json_path, mesh, options.alpha, requested_offset,
      options.offset, offset_retries, preprocessed.stats, validation);
  const bool wrote_smoothed_points = CGAL::IO::write_points(
      smoothed_points_path, preprocessed.points,
      CGAL::parameters::stream_precision(17));

  const auto& tetra_stats = adapter.statistics();
  std::cout << "CGAL version: " << kCgalVersion << "\n"
            << "CGAL Alpha Wrap API: internal adapter (version-pinned)\n"
            << "预处理方法: "
            << rokae_demo::pointCloudPreprocessModeName(
                   preprocessed.stats.mode)
            << "\n"
            << "voxel_size (m): " << preprocessed.stats.voxelization.voxel_size << "\n"
            << "occupied_voxels: " << preprocessed.stats.voxelization.occupied_voxels << "\n"
            << "voxel_points: " << preprocessed.stats.voxelization.voxel_points << "\n"
            << "voxel_point_ratio: " << preprocessed.stats.voxelization.voxel_point_ratio << "\n"
            << "raw_points_outside_voxels: "
            << preprocessed.stats.voxelization.raw_points_outside_voxels << "\n"
            << "平滑邻居数: " << preprocessed.stats.effective_neighbors
            << "\n"
            << "平滑迭代数: " << preprocessed.stats.iterations << "\n"
            << "平均点间距(m): " << preprocessed.stats.average_spacing << "\n"
            << "平均平滑位移(m): " << preprocessed.stats.mean_displacement
            << "\n"
            << "最大平滑位移(m): " << preprocessed.stats.maximum_displacement
            << "\n"
            << "尺度估计时间(ms): "
            << preprocessed.stats.scale_estimation_ms << "\n"
            << "平滑时间(ms): " << preprocessed.stats.smoothing_ms << "\n"
            << "最终 offset(m): " << options.offset << "\n"
            << "offset 重试次数: " << offset_retries << "\n"
            << "输出顶点数: " << mesh.number_of_vertices() << "\n"
            << "输出边数: " << mesh.number_of_edges() << "\n"
            << "输出面数: " << mesh.number_of_faces() << "\n"
            << "是否闭合: " << (validation.closed ? "是" : "否") << "\n"
            << "是否三角网格: " << (validation.triangle_mesh ? "是" : "否")
            << "\n"
            << "是否自交: " << (validation.self_intersects ? "是" : "否")
            << "\n"
            << "是否方向一致且朝外: "
            << ((validation.outward_oriented && validation.bounds_volume)
                    ? "是"
                    : "否")
            << "\n"
            << "ON_BOUNDED_SIDE: " << validation.bounded << "\n"
            << "ON_BOUNDARY: " << validation.boundary << "\n"
            << "ON_UNBOUNDED_SIDE: " << validation.unbounded << "\n"
            << "Wrap input ON_UNBOUNDED_SIDE: "
            << validation.wrap_input_unbounded << "\n"
            << "Finite tetrahedral cells: " << tetra_stats.finite_cells
            << "\n"
            << "Outside tetrahedral cells: " << tetra_stats.outside_cells
            << "\n"
            << "Non-outside retained cells: "
            << tetra_stats.non_outside_cells << "\n"
            << "Tet validation tested: " << volume_validation.tested << "\n"
            << "Tet validation boundary ignored: "
            << volume_validation.boundary_ignored << "\n"
            << "Mesh inside / tet outside: "
            << volume_validation.mesh_inside_tet_outside << "\n"
            << "Mesh outside / tet inside: "
            << volume_validation.mesh_outside_tet_inside << "\n"
            << "运行时间(ms): " << runtime_ms << "\n"
            << "METRICS length_unit=" << rokae_demo::kLengthUnit
            << " alpha=" << options.alpha
            << " requested_offset=" << requested_offset
            << " offset=" << options.offset
            << " offset_retries=" << offset_retries
            << " vertices=" << mesh.number_of_vertices()
            << " edges=" << mesh.number_of_edges()
            << " faces=" << mesh.number_of_faces()
            << " finite_cells=" << tetra_stats.finite_cells
            << " outside_cells=" << tetra_stats.outside_cells
            << " non_outside_cells=" << tetra_stats.non_outside_cells
            << " mesh_inside_tet_outside="
            << volume_validation.mesh_inside_tet_outside
            << " mesh_outside_tet_inside="
            << volume_validation.mesh_outside_tet_inside
            << " runtime_ms=" << runtime_ms
            << " preprocessing_ms=" << preprocessed.stats.total_ms
            << " smooth_neighbors="
            << preprocessed.stats.effective_neighbors
            << " max_smoothing_displacement="
            << preprocessed.stats.maximum_displacement
            << " outside_points=" << validation.unbounded
            << " wrap_input_outside_points="
            << validation.wrap_input_unbounded
            << " closed=" << (validation.closed ? 1 : 0)
            << " triangle=" << (validation.triangle_mesh ? 1 : 0)
            << " self_intersect=" << (validation.self_intersects ? 1 : 0)
            << " oriented="
            << ((validation.outward_oriented && validation.bounds_volume) ? 1
                                                                          : 0)
            << "\n";

  const bool valid =
      validation.closed && validation.triangle_mesh &&
      !validation.self_intersects && validation.outward_oriented &&
      validation.bounds_volume && validation.unbounded == 0 &&
      validation.wrap_input_unbounded == 0 &&
      volume_validation.mesh_inside_tet_outside == 0 &&
      volume_validation.mesh_outside_tet_inside == 0;
  const bool all_valid = valid && probes_passed;
  if(!wrote_off || !wrote_ply || !wrote_json || !wrote_smoothed_points)
    std::cerr << "Failed to write one or more requested output files.\n";
  if(validation.unbounded != 0)
    std::cerr << "Containment acceptance failed: input points lie outside.\n";

  return (all_valid && wrote_off && wrote_ply && wrote_json &&
          wrote_smoothed_points)
             ? EXIT_SUCCESS
             : EXIT_FAILURE;
}
