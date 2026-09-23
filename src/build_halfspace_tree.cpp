#include <rokae_demo/halfspace_logic_tree.hpp>
#include <rokae_demo/units.hpp>
#include <rokae_demo/alpha_tet_pipeline.hpp>
#include <rokae_demo/octree_pipeline.hpp>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/convex_decomposition_3.h>
#include <CGAL/number_utils.h>
#include <CGAL/version.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
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
#include <tuple>
#include <utility>
#include <vector>

namespace fs = std::filesystem;
namespace PMP = CGAL::Polygon_mesh_processing;
namespace HS = rokae_demo::halfspace;

using ExactKernel = CGAL::Exact_predicates_exact_constructions_kernel;
using ExactPoint = ExactKernel::Point_3;
using ExactVector = ExactKernel::Vector_3;
using ExactPolyhedron = CGAL::Polyhedron_3<ExactKernel>;
using NefPolyhedron =
    CGAL::Nef_polyhedron_3<ExactKernel, CGAL::SNC_indexed_items>;
using VolumeIterator = NefPolyhedron::Volume_const_iterator;

using InexactKernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using InexactPoint = InexactKernel::Point_3;
using InexactMesh = CGAL::Surface_mesh<InexactPoint>;
using Primitive = CGAL::AABB_face_graph_triangle_primitive<InexactMesh>;
using AabbTraits = CGAL::AABB_traits<InexactKernel, Primitive>;
using AabbTree = CGAL::AABB_tree<AabbTraits>;

#define ROKAE_STRINGIFY_IMPL(value) #value
#define ROKAE_STRINGIFY(value) ROKAE_STRINGIFY_IMPL(value)
static constexpr const char* kCgalVersion = ROKAE_STRINGIFY(CGAL_VERSION);

struct Options
{
  std::string input_path;
  std::string output_directory;
  double plane_angle_tolerance = 1e-8;
  double plane_distance_tolerance = 1e-8;
  double classification_epsilon = 1e-8;
  double expression_tolerance = 1e-12;
  double volume_tolerance = 1e-8;
  std::size_t validation_samples = 100000;
  int expected_parts = -1;
  int expected_min_parts = -1;
  int expected_planes = -1;
  std::size_t source_point_count = 0;
  double source_alpha = 0.0;
  double source_offset = 0.0;
  struct Probe
  {
    HS::Vec3 point;
    bool expected_inside = false;
  };
  std::vector<Probe> probes;
};

struct Plane
{
  HS::Vec3 normal;
  HS::Vec3 point;
  double offset = 0.0;
  std::vector<int> source_face_ids;

  double evaluate(const HS::Vec3& value) const
  {
    return normal.x * value.x + normal.y * value.y +
           normal.z * value.z - offset;
  }
};

struct ConvexPart
{
  int id = -1;
  ExactPolyhedron mesh;
  std::vector<int> leaf_ids;
  HS::Vec3 centroid;
  double volume = 0.0;
  bool convex = false;
  bool closed = false;
  bool self_intersects = true;
  bool all_vertices_satisfy = false;
  bool centroid_strictly_inside = false;
  std::size_t face_count = 0;
};

struct Timings
{
  double mesh_conversion_ms = 0.0;
  double nef_construction_ms = 0.0;
  double convex_decomposition_ms = 0.0;
  double halfspace_extraction_ms = 0.0;
  double tree_construction_ms = 0.0;
  double simplification_ms = 0.0;
  double validation_ms = 0.0;
  double total_ms = 0.0;
};

struct Confusion
{
  std::size_t mesh_inside_tree_inside = 0;
  std::size_t mesh_inside_tree_outside = 0;
  std::size_t mesh_outside_tree_inside = 0;
  std::size_t mesh_outside_tree_outside = 0;
  std::size_t boundary_ignored = 0;
  std::size_t ambiguous = 0;
};

struct ValidationReport
{
  bool input_closed = false;
  bool input_triangle_mesh = false;
  bool input_self_intersects = true;
  bool input_bounds_volume = false;
  bool input_outward = false;
  bool all_parts_valid = false;
  bool exact_symmetric_difference_empty = false;
  bool raw_simplified_equivalent = false;
  bool special_points_passed = false;
  bool probes_passed = false;
  Confusion confusion;
  double original_volume = 0.0;
  double convex_parts_volume_sum = 0.0;
  double volume_absolute_error = 0.0;
  double maximum_expression_error = 0.0;
  std::size_t tested_points = 0;
  std::vector<std::string> failures;
};

struct Pipeline
{
  Options options;
  ExactPolyhedron input_polyhedron;
  InexactMesh input_mesh;
  NefPolyhedron original_nef;
  std::vector<Plane> original_planes;
  std::vector<ConvexPart> parts;
  std::vector<HS::HalfspaceLeaf> leaves;
  HS::ExpressionPtr raw_tree;
  HS::ExpressionPtr simplified_tree;
  std::size_t reflex_edge_count = 0;
  std::size_t inserted_cut_plane_count = 0;
  std::size_t original_boundary_plane_count = 0;
  std::size_t raw_leaf_reference_count = 0;
  Timings timings;
  ValidationReport validation;
};

static double elapsedMilliseconds(
    const std::chrono::steady_clock::time_point& start)
{
  return std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - start)
      .count();
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

static bool parseInteger(const std::string& text, int& value)
{
  try
  {
    std::size_t parsed = 0;
    value = std::stoi(text, &parsed);
    return parsed == text.size();
  }
  catch(const std::exception&)
  {
    return false;
  }
}

static void printUsage(const char* program)
{
  std::cerr
      << "Usage: " << program << " wrapped.{off,ply} tree_output [options]\n"
      << "  Input coordinates and all length parameters use meters (m).\n"
      << "  --pipeline exact-nef\n"
      << "  --decomposition exact\n"
      << "  --plane-angle-tol VALUE\n"
      << "  --plane-distance-tol VALUE\n"
      << "  --classification-eps VALUE\n"
      << "  --expression-tol VALUE\n"
      << "  --volume-tol VALUE\n"
      << "  --validation-samples COUNT\n"
      << "  --point-count COUNT --alpha VALUE --offset VALUE (benchmark metadata)\n"
      << "  --probe X Y Z inside|outside\n"
      << "Point-cloud alternatives (see README for route-specific options):\n"
      << "  " << program << " points.xyz output --pipeline alpha-tet --voxel-size <m>\n"
      << "  " << program << " points.xyz output --pipeline octree --voxel-size <m>\n"
      << "    Octree: [--octree-pruning rect|exact|none] [--octree-boundary exact|none] "
         "[--max-depth 20] [--validation-samples 1000]\n"
      << "    Optional: --octree-approx hull --max-volume-inflation <ratio> --max-fill-distance <m>\n";
}

static bool parseOptions(int argc, char** argv, Options& options)
{
  if(argc < 3)
  {
    printUsage(argv[0]);
    return false;
  }
  options.input_path = argv[1];
  options.output_directory = argv[2];

  for(int index = 3; index < argc; ++index)
  {
    const std::string option = argv[index];
    if(option == "--probe")
    {
      if(index + 4 >= argc)
      {
        std::cerr << "--probe requires X Y Z inside|outside\n";
        return false;
      }
      Options::Probe probe;
      if(!parseFiniteDouble(argv[++index], probe.point.x) ||
         !parseFiniteDouble(argv[++index], probe.point.y) ||
         !parseFiniteDouble(argv[++index], probe.point.z))
        return false;
      const std::string expected = argv[++index];
      if(expected != "inside" && expected != "outside")
        return false;
      probe.expected_inside = expected == "inside";
      options.probes.push_back(probe);
      continue;
    }

    if(index + 1 >= argc)
    {
      std::cerr << "Missing value for " << option << "\n";
      return false;
    }
    const std::string value = argv[++index];
    if(option == "--pipeline")
    {
      if(value != "exact-nef")
      {
        std::cerr << "This route requires --pipeline exact-nef.\n";
        return false;
      }
    }
    else if(option == "--decomposition")
    {
      if(value != "exact")
      {
        std::cerr << "Only --decomposition exact is supported.\n";
        return false;
      }
    }
    else if(option == "--plane-angle-tol")
    {
      if(!parsePositiveDouble(value, options.plane_angle_tolerance))
        return false;
    }
    else if(option == "--plane-distance-tol")
    {
      if(!parsePositiveDouble(value, options.plane_distance_tolerance))
        return false;
    }
    else if(option == "--classification-eps")
    {
      if(!parsePositiveDouble(value, options.classification_epsilon))
        return false;
    }
    else if(option == "--expression-tol")
    {
      if(!parsePositiveDouble(value, options.expression_tolerance))
        return false;
    }
    else if(option == "--volume-tol")
    {
      if(!parsePositiveDouble(value, options.volume_tolerance))
        return false;
    }
    else if(option == "--validation-samples")
    {
      if(!parseNonnegativeInteger(value, options.validation_samples))
        return false;
    }
    else if(option == "--point-count")
    {
      if(!parseNonnegativeInteger(value, options.source_point_count))
        return false;
    }
    else if(option == "--alpha")
    {
      if(!parsePositiveDouble(value, options.source_alpha))
        return false;
    }
    else if(option == "--offset")
    {
      if(!parsePositiveDouble(value, options.source_offset))
        return false;
    }
    else if(option == "--expect-parts")
    {
      if(!parseInteger(value, options.expected_parts))
        return false;
    }
    else if(option == "--expect-planes")
    {
      if(!parseInteger(value, options.expected_planes))
        return false;
    }
    else if(option == "--expect-min-parts")
    {
      if(!parseInteger(value, options.expected_min_parts))
        return false;
    }
    else
    {
      std::cerr << "Unknown option: " << option << "\n";
      return false;
    }
  }
  return true;
}

static HS::Vec3 toVec3(const ExactPoint& point)
{
  return {CGAL::to_double(point.x()), CGAL::to_double(point.y()),
          CGAL::to_double(point.z())};
}

static InexactPoint toInexactPoint(const ExactPoint& point)
{
  return {CGAL::to_double(point.x()), CGAL::to_double(point.y()),
          CGAL::to_double(point.z())};
}

static HS::Vec3 add(const HS::Vec3& lhs, const HS::Vec3& rhs)
{
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

static HS::Vec3 scale(const HS::Vec3& value, double factor)
{
  return {value.x * factor, value.y * factor, value.z * factor};
}

static double dot(const HS::Vec3& lhs, const HS::Vec3& rhs)
{
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

static bool sameDirectedPlane(const Plane& lhs,
                              const Plane& rhs,
                              double angle_tolerance,
                              double distance_tolerance)
{
  return dot(lhs.normal, rhs.normal) >= 1.0 - angle_tolerance &&
         std::abs(lhs.offset - rhs.offset) <= distance_tolerance;
}

static bool sameDirectedPlane(const HS::HalfspaceLeaf& lhs,
                              const Plane& rhs,
                              double angle_tolerance,
                              double distance_tolerance)
{
  return dot(lhs.normal, rhs.normal) >= 1.0 - angle_tolerance &&
         std::abs(lhs.offset - rhs.offset) <= distance_tolerance;
}

static std::vector<ExactPoint> facetPoints(
    ExactPolyhedron::Facet_const_handle facet)
{
  std::vector<ExactPoint> points;
  auto circulator = facet->facet_begin();
  const auto begin = circulator;
  do
  {
    points.push_back(circulator->vertex()->point());
    ++circulator;
  } while(circulator != begin);
  return points;
}

static HS::Vec3 vertexAverage(const ExactPolyhedron& polyhedron)
{
  if(polyhedron.size_of_vertices() == 0)
    throw std::runtime_error("Cannot average an empty vertex set");
  HS::Vec3 result;
  for(auto vertex = polyhedron.vertices_begin();
      vertex != polyhedron.vertices_end(); ++vertex)
    result = add(result, toVec3(vertex->point()));
  return scale(result, 1.0 / static_cast<double>(polyhedron.size_of_vertices()));
}

static Plane planeFromFacet(ExactPolyhedron::Facet_const_handle facet,
                            const HS::Vec3* interior_point,
                            int source_face_id)
{
  const std::vector<ExactPoint> points = facetPoints(facet);
  if(points.size() < 3)
    throw std::runtime_error("Facet has fewer than three vertices");

  ExactVector exact_normal;
  bool found = false;
  for(std::size_t i = 1; i + 1 < points.size(); ++i)
  {
    exact_normal =
        CGAL::cross_product(points[i] - points[0], points[i + 1] - points[0]);
    if(exact_normal.squared_length() != ExactKernel::FT(0))
    {
      found = true;
      break;
    }
  }
  if(!found)
    throw std::runtime_error("Degenerate facet has no valid normal");

  double nx = CGAL::to_double(exact_normal.x());
  double ny = CGAL::to_double(exact_normal.y());
  double nz = CGAL::to_double(exact_normal.z());
  const double normal_length = std::sqrt(nx * nx + ny * ny + nz * nz);
  if(!std::isfinite(normal_length) ||
     normal_length <= std::numeric_limits<double>::epsilon())
    throw std::runtime_error("Facet normal cannot be normalized");
  nx /= normal_length;
  ny /= normal_length;
  nz /= normal_length;

  Plane plane;
  plane.normal = {nx, ny, nz};
  for(const ExactPoint& point : points)
    plane.point = add(plane.point, toVec3(point));
  plane.point = scale(plane.point, 1.0 / static_cast<double>(points.size()));
  plane.offset = dot(plane.normal, plane.point);
  if(source_face_id >= 0)
    plane.source_face_ids.push_back(source_face_id);

  if(interior_point != nullptr && plane.evaluate(*interior_point) >= 0.0)
  {
    plane.normal = scale(plane.normal, -1.0);
    plane.offset = -plane.offset;
  }
  return plane;
}

static ExactPolyhedron triangulatedCopy(const ExactPolyhedron& input)
{
  ExactPolyhedron result = input;
  if(!CGAL::is_triangle_mesh(result) && !PMP::triangulate_faces(result))
    throw std::runtime_error("Failed to triangulate polyhedron faces");
  return result;
}

static InexactMesh toInexactMesh(const ExactPolyhedron& input)
{
  const ExactPolyhedron triangulated = triangulatedCopy(input);
  InexactMesh result;
  std::map<const void*, InexactMesh::Vertex_index> vertex_map;
  for(auto vertex = triangulated.vertices_begin();
      vertex != triangulated.vertices_end(); ++vertex)
  {
    vertex_map[static_cast<const void*>(&*vertex)] =
        result.add_vertex(toInexactPoint(vertex->point()));
  }

  for(auto facet = triangulated.facets_begin();
      facet != triangulated.facets_end(); ++facet)
  {
    std::vector<InexactMesh::Vertex_index> vertices;
    auto circulator = facet->facet_begin();
    const auto begin = circulator;
    do
    {
      vertices.push_back(
          vertex_map.at(static_cast<const void*>(&*circulator->vertex())));
      ++circulator;
    } while(circulator != begin);
    if(result.add_face(vertices) == InexactMesh::null_face())
      throw std::runtime_error("Failed to convert a face to Surface_mesh");
  }
  return result;
}

static void writeInexactOff(const std::string& path,
                            const ExactPolyhedron& polyhedron)
{
  const InexactMesh mesh = toInexactMesh(polyhedron);
  if(!CGAL::IO::write_polygon_mesh(
         path, mesh, CGAL::parameters::stream_precision(17)))
    throw std::runtime_error("Failed to write convex part: " + path);
}

static void validateAndOrientInput(Pipeline& pipeline)
{
  if(!CGAL::IO::read_polygon_mesh(pipeline.options.input_path,
                                  pipeline.input_polyhedron) ||
     pipeline.input_polyhedron.empty())
    throw std::runtime_error("Failed to read a non-empty OFF/PLY mesh");

  pipeline.validation.input_triangle_mesh =
      CGAL::is_triangle_mesh(pipeline.input_polyhedron);
  pipeline.validation.input_closed =
      CGAL::is_closed(pipeline.input_polyhedron);
  if(!pipeline.validation.input_triangle_mesh)
    throw std::runtime_error("Input mesh must be triangulated");
  if(!pipeline.validation.input_closed)
    throw std::runtime_error("Input mesh must be closed");
  pipeline.validation.input_self_intersects =
      PMP::does_self_intersect(pipeline.input_polyhedron);
  if(pipeline.validation.input_self_intersects)
    throw std::runtime_error("Input mesh self-intersects");

  PMP::orient_to_bound_a_volume(
      pipeline.input_polyhedron,
      CGAL::parameters::outward_orientation(true));
  pipeline.validation.input_bounds_volume =
      PMP::does_bound_a_volume(pipeline.input_polyhedron);
  pipeline.validation.input_outward =
      PMP::is_outward_oriented(pipeline.input_polyhedron);
  if(!pipeline.validation.input_bounds_volume ||
     !pipeline.validation.input_outward)
    throw std::runtime_error("Input mesh does not bound an outward volume");

  pipeline.input_mesh = toInexactMesh(pipeline.input_polyhedron);

  pipeline.original_planes.clear();
  int face_id = 0;
  for(auto facet = pipeline.input_polyhedron.facets_begin();
      facet != pipeline.input_polyhedron.facets_end(); ++facet, ++face_id)
    pipeline.original_planes.push_back(
        planeFromFacet(facet, nullptr, face_id));
}

static std::size_t countReflexEdges(const InexactMesh& mesh,
                                    double epsilon)
{
  std::vector<Plane> face_planes(mesh.number_of_faces());
  for(const auto face : mesh.faces())
  {
    std::vector<HS::Vec3> points;
    for(const auto vertex :
        CGAL::vertices_around_face(mesh.halfedge(face), mesh))
    {
      const InexactPoint& point = mesh.point(vertex);
      points.push_back({point.x(), point.y(), point.z()});
    }
    const HS::Vec3 a{points[1].x - points[0].x,
                     points[1].y - points[0].y,
                     points[1].z - points[0].z};
    const HS::Vec3 b{points[2].x - points[0].x,
                     points[2].y - points[0].y,
                     points[2].z - points[0].z};
    HS::Vec3 normal{a.y * b.z - a.z * b.y,
                    a.z * b.x - a.x * b.z,
                    a.x * b.y - a.y * b.x};
    const double length = std::sqrt(dot(normal, normal));
    normal = scale(normal, 1.0 / length);
    face_planes[face.idx()] = {normal, points[0], dot(normal, points[0]), {}};
  }

  std::size_t reflex_edges = 0;
  for(const auto edge : mesh.edges())
  {
    const auto h0 = mesh.halfedge(edge);
    const auto h1 = mesh.opposite(h0);
    const auto f0 = mesh.face(h0);
    const auto f1 = mesh.face(h1);
    if(f0 == InexactMesh::null_face() || f1 == InexactMesh::null_face())
      continue;
    const auto opposite0 = mesh.target(mesh.next(h0));
    const auto opposite1 = mesh.target(mesh.next(h1));
    const InexactPoint& p0 = mesh.point(opposite0);
    const InexactPoint& p1 = mesh.point(opposite1);
    const HS::Vec3 v0{p0.x(), p0.y(), p0.z()};
    const HS::Vec3 v1{p1.x(), p1.y(), p1.z()};
    if(face_planes[f0.idx()].evaluate(v1) > epsilon ||
       face_planes[f1.idx()].evaluate(v0) > epsilon)
      ++reflex_edges;
  }
  return reflex_edges;
}

static void extractConvexParts(Pipeline& pipeline,
                               NefPolyhedron& decomposed_nef)
{
  int part_id = 0;
  VolumeIterator volume = decomposed_nef.volumes_begin();
  if(volume != decomposed_nef.volumes_end())
    ++volume;  // The first volume is the unbounded exterior.
  for(; volume != decomposed_nef.volumes_end(); ++volume)
  {
    if(!volume->mark() || volume->shells_begin() == volume->shells_end())
      continue;
    ConvexPart part;
    part.id = part_id++;
    decomposed_nef.convert_inner_shell_to_polyhedron(
        volume->shells_begin(), part.mesh);
    if(part.mesh.empty())
      throw std::runtime_error("Convex decomposition produced an empty part");

    part.mesh = triangulatedCopy(part.mesh);
    if(!PMP::is_outward_oriented(part.mesh))
      PMP::reverse_face_orientations(part.mesh);
    part.closed = CGAL::is_closed(part.mesh);
    part.self_intersects = PMP::does_self_intersect(part.mesh);
    const NefPolyhedron part_nef(part.mesh);
    part.convex = part_nef.is_convex();
    part.centroid = vertexAverage(part.mesh);
    part.volume = std::abs(CGAL::to_double(PMP::volume(part.mesh)));
    part.face_count = part.mesh.size_of_facets();
    if(!part.closed || part.self_intersects || !part.convex ||
       part.volume <= 0.0)
      throw std::runtime_error("A decomposed part failed geometric validity");
    pipeline.parts.push_back(std::move(part));
  }
  if(pipeline.parts.empty())
    throw std::runtime_error("Exact convex decomposition returned no parts");
}

static std::vector<Plane> uniquePartPlanes(const ConvexPart& part,
                                           const Options& options)
{
  std::vector<Plane> planes;
  for(auto facet = part.mesh.facets_begin(); facet != part.mesh.facets_end();
      ++facet)
  {
    const Plane candidate = planeFromFacet(facet, &part.centroid, -1);
    bool duplicate = false;
    for(const Plane& plane : planes)
    {
      if(sameDirectedPlane(candidate, plane, options.plane_angle_tolerance,
                           options.plane_distance_tolerance))
      {
        duplicate = true;
        break;
      }
    }
    if(!duplicate)
      planes.push_back(candidate);
  }
  return planes;
}

static void extractHalfspaces(Pipeline& pipeline)
{
  for(ConvexPart& part : pipeline.parts)
  {
    const std::vector<Plane> planes =
        uniquePartPlanes(part, pipeline.options);
    if(planes.size() < 4)
      throw std::runtime_error(
          "A full-dimensional convex part has fewer than four planes");

    part.all_vertices_satisfy = true;
    part.centroid_strictly_inside = true;
    for(const Plane& plane : planes)
    {
      if(plane.evaluate(part.centroid) >=
         -pipeline.options.classification_epsilon)
        part.centroid_strictly_inside = false;
      for(auto vertex = part.mesh.vertices_begin();
          vertex != part.mesh.vertices_end(); ++vertex)
      {
        if(plane.evaluate(toVec3(vertex->point())) >
           pipeline.options.classification_epsilon)
          part.all_vertices_satisfy = false;
      }
    }
    if(!part.all_vertices_satisfy || !part.centroid_strictly_inside)
      throw std::runtime_error(
          "Convex-part vertices or centroid violate extracted halfspaces");

    for(Plane plane : planes)
    {
      std::vector<int> matched_source_faces;
      for(const Plane& original : pipeline.original_planes)
      {
        if(sameDirectedPlane(plane, original,
                             pipeline.options.plane_angle_tolerance,
                             pipeline.options.plane_distance_tolerance))
        {
          matched_source_faces.insert(matched_source_faces.end(),
                                      original.source_face_ids.begin(),
                                      original.source_face_ids.end());
        }
      }
      std::sort(matched_source_faces.begin(), matched_source_faces.end());
      matched_source_faces.erase(
          std::unique(matched_source_faces.begin(),
                      matched_source_faces.end()),
          matched_source_faces.end());

      int leaf_id = -1;
      for(HS::HalfspaceLeaf& leaf : pipeline.leaves)
      {
        if(sameDirectedPlane(leaf, plane,
                             pipeline.options.plane_angle_tolerance,
                             pipeline.options.plane_distance_tolerance))
        {
          leaf_id = leaf.id;
          if(std::find(leaf.convex_part_ids.begin(),
                       leaf.convex_part_ids.end(),
                       part.id) == leaf.convex_part_ids.end())
            leaf.convex_part_ids.push_back(part.id);
          leaf.source_face_ids.insert(leaf.source_face_ids.end(),
                                      matched_source_faces.begin(),
                                      matched_source_faces.end());
          std::sort(leaf.source_face_ids.begin(),
                    leaf.source_face_ids.end());
          leaf.source_face_ids.erase(
              std::unique(leaf.source_face_ids.begin(),
                          leaf.source_face_ids.end()),
              leaf.source_face_ids.end());
          break;
        }
      }
      if(leaf_id < 0)
      {
        HS::HalfspaceLeaf leaf;
        leaf.id = static_cast<int>(pipeline.leaves.size());
        leaf.normal = plane.normal;
        leaf.point = plane.point;
        leaf.offset = plane.offset;
        leaf.convex_part_ids.push_back(part.id);
        leaf.source_face_ids = matched_source_faces;
        leaf.source_type = matched_source_faces.empty()
                               ? "decomposition_cut"
                               : "original_boundary";
        leaf_id = leaf.id;
        pipeline.leaves.push_back(std::move(leaf));
      }
      part.leaf_ids.push_back(leaf_id);
    }
    std::sort(part.leaf_ids.begin(), part.leaf_ids.end());
  }

  for(const auto& leaf : pipeline.leaves)
  {
    if(leaf.source_type == "decomposition_cut")
      ++pipeline.inserted_cut_plane_count;
    else if(leaf.source_type == "original_boundary")
      ++pipeline.original_boundary_plane_count;
  }
}

static void constructTrees(Pipeline& pipeline)
{
  std::vector<HS::ExpressionPtr> part_expressions;
  for(const ConvexPart& part : pipeline.parts)
  {
    if(part.leaf_ids.empty())
      throw std::runtime_error("A convex part has no halfspace leaves");
    std::vector<HS::ExpressionPtr> leaves;
    for(const int leaf_id : part.leaf_ids)
      leaves.push_back(HS::makeLeaf(leaf_id));
    pipeline.raw_leaf_reference_count += leaves.size();
    part_expressions.push_back(
        leaves.size() == 1 ? leaves.front()
                           : HS::makeNode(HS::TreeOperator::MAX, leaves));
  }
  pipeline.raw_tree =
      part_expressions.size() == 1
          ? part_expressions.front()
          : HS::makeNode(HS::TreeOperator::MIN, part_expressions);
}

static std::array<double, 6> meshBoundingBox(const InexactMesh& mesh)
{
  std::array<double, 6> bbox{
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity()};
  for(const auto vertex : mesh.vertices())
  {
    const auto& point = mesh.point(vertex);
    bbox[0] = std::min(bbox[0], point.x());
    bbox[1] = std::min(bbox[1], point.y());
    bbox[2] = std::min(bbox[2], point.z());
    bbox[3] = std::max(bbox[3], point.x());
    bbox[4] = std::max(bbox[4], point.y());
    bbox[5] = std::max(bbox[5], point.z());
  }
  return bbox;
}

static void recordComparison(Pipeline& pipeline,
                             const HS::Vec3& point,
                             const std::string& label,
                             const CGAL::Bounded_side mesh_side,
                             double distance_squared)
{
  const double raw =
      HS::evaluate(pipeline.raw_tree, pipeline.leaves, point);
  const double simplified =
      HS::evaluate(pipeline.simplified_tree, pipeline.leaves, point);
  const double expression_error = std::abs(raw - simplified);
  pipeline.validation.maximum_expression_error =
      std::max(pipeline.validation.maximum_expression_error,
               expression_error);
  if(expression_error > pipeline.options.expression_tolerance ||
     (raw <= pipeline.options.classification_epsilon) !=
         (simplified <= pipeline.options.classification_epsilon))
  {
    std::ostringstream failure;
    failure << "Raw/simplified mismatch at " << label << " ["
            << point.x << ", " << point.y << ", " << point.z
            << "]: raw=" << raw << ", simplified=" << simplified;
    pipeline.validation.failures.push_back(failure.str());
  }

  const double epsilon_squared =
      pipeline.options.classification_epsilon *
      pipeline.options.classification_epsilon;
  if(mesh_side == CGAL::ON_BOUNDARY ||
     distance_squared <= epsilon_squared)
  {
    ++pipeline.validation.confusion.boundary_ignored;
    ++pipeline.validation.tested_points;
    return;
  }
  if(std::abs(raw) <= pipeline.options.classification_epsilon)
  {
    ++pipeline.validation.confusion.ambiguous;
    ++pipeline.validation.tested_points;
    return;
  }

  const bool mesh_inside = mesh_side == CGAL::ON_BOUNDED_SIDE;
  const bool tree_inside = raw <= pipeline.options.classification_epsilon;
  if(mesh_inside && tree_inside)
    ++pipeline.validation.confusion.mesh_inside_tree_inside;
  else if(mesh_inside)
    ++pipeline.validation.confusion.mesh_inside_tree_outside;
  else if(tree_inside)
    ++pipeline.validation.confusion.mesh_outside_tree_inside;
  else
    ++pipeline.validation.confusion.mesh_outside_tree_outside;
  ++pipeline.validation.tested_points;
}

static void validatePipeline(Pipeline& pipeline)
{
  pipeline.validation.all_parts_valid = true;
  for(const ConvexPart& part : pipeline.parts)
  {
    pipeline.validation.all_parts_valid &=
        part.closed && !part.self_intersects && part.convex &&
        part.all_vertices_satisfy && part.centroid_strictly_inside;

    double centroid_max = -std::numeric_limits<double>::infinity();
    for(const int leaf_id : part.leaf_ids)
      centroid_max =
          std::max(centroid_max,
                   HS::evaluate(HS::makeLeaf(leaf_id), pipeline.leaves,
                                part.centroid));
    if(centroid_max >= -pipeline.options.classification_epsilon)
      pipeline.validation.failures.push_back(
          "Convex-part centroid did not evaluate strictly inside");

    for(const int leaf_id : part.leaf_ids)
    {
      const HS::HalfspaceLeaf& leaf = pipeline.leaves.at(leaf_id);
      const double margin =
          -(leaf.normal.x * part.centroid.x +
            leaf.normal.y * part.centroid.y +
            leaf.normal.z * part.centroid.z - leaf.offset);
      const HS::Vec3 outside =
          add(part.centroid,
              scale(leaf.normal,
                    margin + 10.0 * pipeline.options.classification_epsilon));
      double outside_max = -std::numeric_limits<double>::infinity();
      for(const int candidate_id : part.leaf_ids)
        outside_max =
            std::max(outside_max,
                     HS::evaluate(HS::makeLeaf(candidate_id),
                                  pipeline.leaves, outside));
      if(outside_max <= pipeline.options.classification_epsilon)
        pipeline.validation.failures.push_back(
            "Convex-part exterior sample evaluated inside");
    }
  }

  ExactPolyhedron original_tri =
      triangulatedCopy(pipeline.input_polyhedron);
  pipeline.validation.original_volume =
      std::abs(CGAL::to_double(PMP::volume(original_tri)));
  for(const ConvexPart& part : pipeline.parts)
    pipeline.validation.convex_parts_volume_sum += part.volume;
  pipeline.validation.volume_absolute_error =
      std::abs(pipeline.validation.original_volume -
               pipeline.validation.convex_parts_volume_sum);
  const double scaled_volume_tolerance =
      pipeline.options.volume_tolerance *
      std::max(1.0, pipeline.validation.original_volume);
  if(pipeline.validation.volume_absolute_error > scaled_volume_tolerance)
    pipeline.validation.failures.push_back(
        "Convex-part volume sum differs from original volume");

  NefPolyhedron reconstructed;
  for(const ConvexPart& part : pipeline.parts)
    reconstructed += NefPolyhedron(part.mesh);
  pipeline.validation.exact_symmetric_difference_empty =
      pipeline.original_nef.symmetric_difference(reconstructed).is_empty();
  if(!pipeline.validation.exact_symmetric_difference_empty)
    pipeline.validation.failures.push_back(
        "Exact Nef symmetric difference is non-empty");

  CGAL::Side_of_triangle_mesh<InexactMesh, InexactKernel> side(
      pipeline.input_mesh);
  AabbTree distance_tree(faces(pipeline.input_mesh).first,
                         faces(pipeline.input_mesh).second,
                         pipeline.input_mesh);
  distance_tree.accelerate_distance_queries();
  const auto classify = [&](const HS::Vec3& point,
                            const std::string& label) {
    const InexactPoint query(point.x, point.y, point.z);
    recordComparison(pipeline, point, label, side(query),
                     distance_tree.squared_distance(query));
  };

  const auto bbox = meshBoundingBox(pipeline.input_mesh);
  const double dx = bbox[3] - bbox[0];
  const double dy = bbox[4] - bbox[1];
  const double dz = bbox[5] - bbox[2];
  const double diagonal = std::sqrt(dx * dx + dy * dy + dz * dz);
  const double expansion = 0.1 * diagonal;
  std::mt19937_64 generator(20260727);
  std::uniform_real_distribution<double> random_x(bbox[0] - expansion,
                                                   bbox[3] + expansion);
  std::uniform_real_distribution<double> random_y(bbox[1] - expansion,
                                                   bbox[4] + expansion);
  std::uniform_real_distribution<double> random_z(bbox[2] - expansion,
                                                   bbox[5] + expansion);
  for(std::size_t i = 0; i < pipeline.options.validation_samples; ++i)
    classify({random_x(generator), random_y(generator), random_z(generator)},
             "random_" + std::to_string(i));

  int vertex_id = 0;
  for(auto vertex = pipeline.input_polyhedron.vertices_begin();
      vertex != pipeline.input_polyhedron.vertices_end();
      ++vertex, ++vertex_id)
    classify(toVec3(vertex->point()),
             "original_vertex_" + std::to_string(vertex_id));

  for(const ConvexPart& part : pipeline.parts)
    classify(part.centroid,
             "convex_centroid_" + std::to_string(part.id));

  const double normal_step =
      std::max(10.0 * pipeline.options.classification_epsilon,
               diagonal * 1e-7);
  for(std::size_t face_id = 0;
      face_id < pipeline.original_planes.size(); ++face_id)
  {
    const Plane& plane = pipeline.original_planes[face_id];
    classify(plane.point, "original_face_center_" + std::to_string(face_id));
    classify(add(plane.point, scale(plane.normal, normal_step)),
             "original_face_outside_" + std::to_string(face_id));
    classify(add(plane.point, scale(plane.normal, -normal_step)),
             "original_face_inside_" + std::to_string(face_id));
  }
  for(const auto& leaf : pipeline.leaves)
  {
    if(leaf.source_type == "decomposition_cut")
      classify(leaf.point, "cut_plane_" + std::to_string(leaf.id));
  }

  pipeline.validation.probes_passed = true;
  for(std::size_t probe_id = 0; probe_id < pipeline.options.probes.size();
      ++probe_id)
  {
    const auto& probe = pipeline.options.probes[probe_id];
    classify(probe.point, "probe_" + std::to_string(probe_id));
    const bool tree_inside =
        HS::evaluateInside(pipeline.raw_tree, pipeline.leaves, probe.point,
                           pipeline.options.classification_epsilon);
    if(tree_inside != probe.expected_inside)
    {
      pipeline.validation.probes_passed = false;
      pipeline.validation.failures.push_back(
          "A user-specified inside/outside probe failed");
    }
  }

  pipeline.validation.raw_simplified_equivalent =
      pipeline.validation.maximum_expression_error <=
      pipeline.options.expression_tolerance;
  pipeline.validation.special_points_passed =
      pipeline.validation.confusion.mesh_inside_tree_outside == 0 &&
      pipeline.validation.confusion.mesh_outside_tree_inside == 0;
}

static std::string jsonBoolean(bool value)
{
  return value ? "true" : "false";
}

static void writeVec3(std::ostream& out, const HS::Vec3& value)
{
  out << "[" << value.x << ", " << value.y << ", " << value.z << "]";
}

static void writeIntegerArray(std::ostream& out,
                              const std::vector<int>& values)
{
  out << "[";
  for(std::size_t i = 0; i < values.size(); ++i)
  {
    if(i != 0)
      out << ", ";
    out << values[i];
  }
  out << "]";
}

static void writeLeavesArray(std::ostream& out,
                             const std::vector<HS::HalfspaceLeaf>& leaves,
                             int indentation)
{
  const std::string pad(indentation, ' ');
  out << "[\n";
  for(std::size_t i = 0; i < leaves.size(); ++i)
  {
    const auto& leaf = leaves[i];
    out << pad << "  {\"id\": " << leaf.id << ", \"normal\": ";
    writeVec3(out, leaf.normal);
    out << ", \"point\": ";
    writeVec3(out, leaf.point);
    out << ", \"offset\": " << leaf.offset << ", \"convex_part_ids\": ";
    writeIntegerArray(out, leaf.convex_part_ids);
    out << ", \"source_face_ids\": ";
    writeIntegerArray(out, leaf.source_face_ids);
    out << ", \"source_type\": \"" << leaf.source_type << "\"}";
    out << (i + 1 == leaves.size() ? "\n" : ",\n");
  }
  out << pad << "]";
}

static void writeTreeJson(const fs::path& path,
                          const Pipeline& pipeline,
                          const HS::ExpressionPtr& expression,
                          const std::string& tree_kind)
{
  std::ofstream out(path);
  if(!out)
    throw std::runtime_error("Failed to write " + path.string());
  out << std::setprecision(17);
  const HS::SerializedTree tree = HS::serializeTree(expression);
  out << "{\n"
      << "  \"version\": 1,\n"
      << "  \"length_unit\": \"" << rokae_demo::kLengthUnit << "\",\n"
      << "  \"halfspace_offset_unit\": \"" << rokae_demo::kLengthUnit
      << "\",\n"
      << "  \"tree_kind\": \"" << tree_kind << "\",\n"
      << "  \"representation\": \"convex_decomposition_dnf\",\n"
      << "  \"sign_convention\": {\n"
      << "    \"leaf\": \"psi(x) = n^T x - offset = n^T(x-w)\",\n"
      << "    \"inside_leaf\": \"psi(x) <= 0\",\n"
      << "    \"inside_root\": \"value <= 0\",\n"
      << "    \"outside_root\": \"value > 0\"\n"
      << "  },\n"
      << "  \"decomposition\": {\"type\": \"exact\", "
      << "\"convex_part_count\": " << pipeline.parts.size() << "},\n"
      << "  \"leaves\": ";
  writeLeavesArray(out, pipeline.leaves, 2);
  out << ",\n  \"nodes\": [\n";
  for(std::size_t i = 0; i < tree.nodes.size(); ++i)
  {
    const auto& node = tree.nodes[i];
    out << "    {\"id\": " << node.id << ", \"type\": \""
        << HS::operatorName(node.op) << "\", \"boolean_alias\": \""
        << HS::booleanAlias(node.op) << "\", \"leaf_id\": "
        << node.leaf_id << ", \"children\": ";
    writeIntegerArray(out, node.children);
    out << "}" << (i + 1 == tree.nodes.size() ? "\n" : ",\n");
  }
  out << "  ],\n  \"root_id\": " << tree.root_id << "\n}\n";
}

static void writeConvexPartsJson(const fs::path& path,
                                 const Pipeline& pipeline)
{
  std::ofstream out(path);
  if(!out)
    throw std::runtime_error("Failed to write convex_parts.json");
  out << std::setprecision(17) << "{\n  \"units\": {\"length\": \""
      << rokae_demo::kLengthUnit << "\", \"volume\": \""
      << rokae_demo::kVolumeUnit << "\"},\n  \"parts\": [\n";
  for(std::size_t i = 0; i < pipeline.parts.size(); ++i)
  {
    const auto& part = pipeline.parts[i];
    std::ostringstream filename;
    filename << "convex_parts/convex_part_" << std::setw(3)
             << std::setfill('0') << part.id << ".off";
    out << "    {\"part_id\": " << part.id << ", \"mesh_file\": \""
        << filename.str() << "\", \"vertex_count\": "
        << part.mesh.size_of_vertices() << ", \"face_count\": "
        << part.face_count << ", \"unique_plane_count\": "
        << part.leaf_ids.size() << ", \"volume\": " << part.volume
        << ", \"centroid\": ";
    writeVec3(out, part.centroid);
    out << ", \"leaf_ids\": ";
    writeIntegerArray(out, part.leaf_ids);
    out << ", \"convexity_check\": " << jsonBoolean(part.convex)
        << ", \"closed_check\": " << jsonBoolean(part.closed)
        << ", \"self_intersection_check\": "
        << jsonBoolean(part.self_intersects) << "}";
    out << (i + 1 == pipeline.parts.size() ? "\n" : ",\n");
  }
  out << "  ]\n}\n";
}

static void writeLeavesJson(const fs::path& path,
                            const Pipeline& pipeline)
{
  std::ofstream out(path);
  if(!out)
    throw std::runtime_error("Failed to write halfspace_leaves.json");
  out << std::setprecision(17)
      << "{\n  \"length_unit\": \"" << rokae_demo::kLengthUnit << "\",\n"
      << "  \"halfspace_offset_unit\": \"" << rokae_demo::kLengthUnit
      << "\",\n  \"sign_convention\": \"psi(x) = n^T x - offset; "
         "inside iff psi <= 0\",\n"
      << "  \"leaves\": ";
  writeLeavesArray(out, pipeline.leaves, 2);
  out << "\n}\n";
}

static void writeExpressions(const fs::path& path,
                             const HS::ExpressionPtr& expression)
{
  std::ofstream out(path);
  if(!out)
    throw std::runtime_error("Failed to write expression file");
  out << "inside(x) =\n" << HS::booleanExpression(expression, 0)
      << "\n\nphi(x) =\n" << HS::scalarExpression(expression, 0) << "\n";
}

static void writeDot(const fs::path& path,
                     const Pipeline& pipeline,
                     const HS::ExpressionPtr& expression)
{
  const HS::SerializedTree tree = HS::serializeTree(expression);
  std::ofstream out(path);
  if(!out)
    throw std::runtime_error("Failed to write logic_tree.dot");
  out << std::setprecision(8) << "digraph HalfspaceTree {\n"
      << "  rankdir=TB;\n"
      << "  node [fontname=\"Helvetica\"];\n";
  for(const auto& node : tree.nodes)
  {
    if(node.op == HS::TreeOperator::LEAF)
    {
      const auto& leaf = pipeline.leaves.at(node.leaf_id);
      out << "  n" << node.id << " [shape=box,label=\"L" << leaf.id
          << "\\nn=[" << leaf.normal.x << "," << leaf.normal.y << ","
          << leaf.normal.z << "]\\nd=" << leaf.offset << "\"];\n";
    }
    else
    {
      out << "  n" << node.id << " [shape="
          << (node.op == HS::TreeOperator::MAX ? "ellipse" : "diamond")
          << ",label=\"" << HS::operatorName(node.op) << "\\n"
          << HS::booleanAlias(node.op) << "\"];\n";
    }
    for(const int child : node.children)
      out << "  n" << node.id << " -> n" << child << ";\n";
  }
  out << "}\n";
}

static void writeValidationReport(const fs::path& path,
                                  const Pipeline& pipeline)
{
  const auto raw_serialized = HS::serializeTree(pipeline.raw_tree);
  const auto simplified_serialized =
      HS::serializeTree(pipeline.simplified_tree);
  const auto& report = pipeline.validation;
  std::ofstream out(path);
  if(!out)
    throw std::runtime_error("Failed to write validation_report.json");
  out << std::setprecision(17)
      << "{\n"
      << "  \"cgal_version\": \"" << kCgalVersion << "\",\n"
      << "  \"units\": {\"length\": \"" << rokae_demo::kLengthUnit
      << "\", \"area\": \"" << rokae_demo::kAreaUnit
      << "\", \"volume\": \"" << rokae_demo::kVolumeUnit << "\"},\n"
      << "  \"kernel\": \"Exact_predicates_exact_constructions_kernel\",\n"
      << "  \"decomposition_api\": \"CGAL::convex_decomposition_3\",\n"
      << "  \"input\": {\"faces\": "
      << pipeline.input_polyhedron.size_of_facets()
      << ", \"reflex_edges\": " << pipeline.reflex_edge_count
      << ", \"closed\": " << jsonBoolean(report.input_closed)
      << ", \"triangle_mesh\": "
      << jsonBoolean(report.input_triangle_mesh)
      << ", \"self_intersects\": "
      << jsonBoolean(report.input_self_intersects)
      << ", \"bounds_volume\": "
      << jsonBoolean(report.input_bounds_volume)
      << ", \"outward\": " << jsonBoolean(report.input_outward) << "},\n"
      << "  \"decomposition\": {\"convex_part_count\": "
      << pipeline.parts.size() << ", \"inserted_cut_plane_count\": "
      << pipeline.inserted_cut_plane_count
      << ", \"original_boundary_plane_count\": "
      << pipeline.original_boundary_plane_count
      << ", \"unique_leaf_count\": " << pipeline.leaves.size() << "},\n"
      << "  \"tree\": {\"raw_leaf_reference_count\": "
      << pipeline.raw_leaf_reference_count << ", \"raw_node_count\": "
      << raw_serialized.nodes.size()
      << ", \"simplified_node_count\": "
      << simplified_serialized.nodes.size()
      << ", \"raw_node_reference_count\": "
      << HS::treeNodeReferenceCount(pipeline.raw_tree)
      << ", \"simplified_node_reference_count\": "
      << HS::treeNodeReferenceCount(pipeline.simplified_tree)
      << ", \"raw_tree_depth\": "
      << HS::treeDepth(pipeline.raw_tree)
      << ", \"simplified_tree_depth\": "
      << HS::treeDepth(pipeline.simplified_tree)
      << ", \"maximum_expression_error\": "
      << report.maximum_expression_error
      << ", \"raw_simplified_equivalent\": "
      << jsonBoolean(report.raw_simplified_equivalent) << "},\n"
      << "  \"classification\": {\"tested_points\": "
      << report.tested_points
      << ", \"mesh_inside_tree_inside\": "
      << report.confusion.mesh_inside_tree_inside
      << ", \"mesh_inside_tree_outside\": "
      << report.confusion.mesh_inside_tree_outside
      << ", \"mesh_outside_tree_inside\": "
      << report.confusion.mesh_outside_tree_inside
      << ", \"mesh_outside_tree_outside\": "
      << report.confusion.mesh_outside_tree_outside
      << ", \"boundary_ignored\": " << report.confusion.boundary_ignored
      << ", \"ambiguous\": " << report.confusion.ambiguous << "},\n"
      << "  \"volume\": {\"original\": " << report.original_volume
      << ", \"convex_parts_sum\": "
      << report.convex_parts_volume_sum
      << ", \"absolute_error\": " << report.volume_absolute_error << "},\n"
      << "  \"checks\": {\"all_parts_valid\": "
      << jsonBoolean(report.all_parts_valid)
      << ", \"exact_symmetric_difference_empty\": "
      << jsonBoolean(report.exact_symmetric_difference_empty)
      << ", \"special_points_passed\": "
      << jsonBoolean(report.special_points_passed)
      << ", \"probes_passed\": " << jsonBoolean(report.probes_passed)
      << "},\n"
      << "  \"timings_ms\": {\"mesh_conversion\": "
      << pipeline.timings.mesh_conversion_ms
      << ", \"nef_construction\": "
      << pipeline.timings.nef_construction_ms
      << ", \"convex_decomposition\": "
      << pipeline.timings.convex_decomposition_ms
      << ", \"halfspace_extraction\": "
      << pipeline.timings.halfspace_extraction_ms
      << ", \"tree_construction\": "
      << pipeline.timings.tree_construction_ms
      << ", \"simplification\": " << pipeline.timings.simplification_ms
      << ", \"validation\": " << pipeline.timings.validation_ms
      << ", \"total\": " << pipeline.timings.total_ms << "},\n"
      << "  \"failures\": [";
  for(std::size_t i = 0; i < report.failures.size(); ++i)
  {
    if(i != 0)
      out << ", ";
    out << "\"" << report.failures[i] << "\"";
  }
  out << "]\n}\n";
}

static void writeOutputs(const Pipeline& pipeline)
{
  const fs::path output(pipeline.options.output_directory);
  const fs::path parts_directory = output / "convex_parts";
  fs::create_directories(parts_directory);
  for(const ConvexPart& part : pipeline.parts)
  {
    std::ostringstream name;
    name << "convex_part_" << std::setw(3) << std::setfill('0')
         << part.id << ".off";
    writeInexactOff((parts_directory / name.str()).string(), part.mesh);
  }
  writeConvexPartsJson(output / "convex_parts.json", pipeline);
  writeLeavesJson(output / "halfspace_leaves.json", pipeline);
  writeTreeJson(output / "logic_tree_raw.json", pipeline,
                pipeline.raw_tree, "raw");
  writeTreeJson(output / "logic_tree_simplified.json", pipeline,
                pipeline.simplified_tree, "simplified");
  writeExpressions(output / "logic_expression_raw.txt",
                   pipeline.raw_tree);
  writeExpressions(output / "logic_expression_simplified.txt",
                   pipeline.simplified_tree);
  writeDot(output / "logic_tree.dot", pipeline,
           pipeline.simplified_tree);
  writeValidationReport(output / "validation_report.json", pipeline);

  const auto tree = HS::serializeTree(pipeline.simplified_tree);
  std::ofstream benchmark_csv(output / "benchmark.csv");
  benchmark_csv << "pipeline,length_unit,point_count,alpha,offset,wrap_faces,tetra_count,"
                   "convex_parts,unique_planes,tree_nodes,tree_depth,"
                   "alpha_wrap_ms,tet_extract_ms,agglomeration_ms,"
                   "plane_reduction_ms,tree_build_ms,validation_ms,total_ms,"
                   "false_positive,false_negative,preprocess_mode,"
                   "smooth_neighbors,smooth_iterations,preprocessing_ms,"
                   "smoothing_ms,max_smoothing_displacement,requested_offset,"
                   "offset_retries\n"
                << std::setprecision(17) << "exact-nef,"
                << rokae_demo::kLengthUnit << ','
                << pipeline.options.source_point_count << ','
                << pipeline.options.source_alpha << ','
                << pipeline.options.source_offset << ','
                << pipeline.input_polyhedron.size_of_facets()
                << ",0," << pipeline.parts.size() << ','
                << pipeline.leaves.size() << ',' << tree.nodes.size() << ','
                << HS::treeDepth(pipeline.simplified_tree)
                << ",0,0,0," << pipeline.timings.halfspace_extraction_ms
                << ',' << pipeline.timings.tree_construction_ms << ','
                << pipeline.timings.validation_ms << ','
                << pipeline.timings.total_ms << ','
                << pipeline.validation.confusion.mesh_outside_tree_inside
                << ','
                << pipeline.validation.confusion.mesh_inside_tree_outside
                << ",none,0,0,0,0,0," << pipeline.options.source_offset
                << ",0\n";

  std::ofstream benchmark_json(output / "benchmark.json");
  benchmark_json << std::setprecision(17)
                 << "[\n  {\"pipeline\": \"exact-nef\", "
                 << "\"length_unit\": \"" << rokae_demo::kLengthUnit
                 << "\", "
                 << "\"point_count\": "
                 << pipeline.options.source_point_count
                 << ", \"alpha\": " << pipeline.options.source_alpha
                 << ", \"offset\": " << pipeline.options.source_offset
                 << ", \"wrap_faces\": "
                 << pipeline.input_polyhedron.size_of_facets()
                 << ", \"tetra_count\": 0, \"convex_parts\": "
                 << pipeline.parts.size() << ", \"unique_planes\": "
                 << pipeline.leaves.size() << ", \"tree_nodes\": "
                 << tree.nodes.size() << ", \"tree_depth\": "
                 << HS::treeDepth(pipeline.simplified_tree)
                 << ", \"preprocess_mode\": \"none\", "
                    "\"smooth_neighbors\": 0, \"smooth_iterations\": 0, "
                    "\"preprocessing_ms\": 0, \"smoothing_ms\": 0, "
                    "\"maximum_smoothing_displacement\": 0"
                 << ", \"nef_construction_ms\": "
                 << pipeline.timings.nef_construction_ms
                 << ", \"convex_decomposition_ms\": "
                 << pipeline.timings.convex_decomposition_ms
                 << ", \"plane_reduction_ms\": "
                 << pipeline.timings.halfspace_extraction_ms
                 << ", \"tree_build_ms\": "
                 << pipeline.timings.tree_construction_ms
                 << ", \"validation_ms\": "
                 << pipeline.timings.validation_ms
                 << ", \"total_ms\": " << pipeline.timings.total_ms
                 << ", \"false_positive\": "
                 << pipeline.validation.confusion.mesh_outside_tree_inside
                 << ", \"false_negative\": "
                 << pipeline.validation.confusion.mesh_inside_tree_outside
                 << "}\n]\n";
}

static bool pipelinePassed(const Pipeline& pipeline)
{
  const ValidationReport& report = pipeline.validation;
  bool passed =
      report.input_closed && report.input_triangle_mesh &&
      !report.input_self_intersects && report.input_bounds_volume &&
      report.input_outward && report.all_parts_valid &&
      report.exact_symmetric_difference_empty &&
      report.raw_simplified_equivalent && report.special_points_passed &&
      report.probes_passed && report.failures.empty() &&
      report.volume_absolute_error <=
          pipeline.options.volume_tolerance *
              std::max(1.0, report.original_volume);
  if(pipeline.options.expected_parts >= 0)
    passed &= static_cast<int>(pipeline.parts.size()) ==
              pipeline.options.expected_parts;
  if(pipeline.options.expected_planes >= 0)
    passed &= static_cast<int>(pipeline.leaves.size()) ==
              pipeline.options.expected_planes;
  if(pipeline.options.expected_min_parts >= 0)
    passed &= static_cast<int>(pipeline.parts.size()) >=
              pipeline.options.expected_min_parts;
  return passed;
}

int main(int argc, char** argv)
{
  for(int index = 3; index + 1 < argc; ++index)
  {
    if(std::string(argv[index]) == "--pipeline" &&
       std::string(argv[index + 1]) == "octree")
      return rokae_demo::runOctreePipeline(argc, argv);
    if(std::string(argv[index]) == "--pipeline" &&
       (std::string(argv[index + 1]) == "alpha-tet" ||
        std::string(argv[index + 1]) == "alpha_tet"))
      return rokae_demo::runAlphaTetPipeline(argc, argv);
  }

  const auto total_start = std::chrono::steady_clock::now();
  Pipeline pipeline;
  if(!parseOptions(argc, argv, pipeline.options))
  {
    printUsage(argv[0]);
    return EXIT_FAILURE;
  }

  try
  {
    auto start = std::chrono::steady_clock::now();
    validateAndOrientInput(pipeline);
    pipeline.reflex_edge_count =
        countReflexEdges(pipeline.input_mesh,
                         pipeline.options.classification_epsilon);
    pipeline.timings.mesh_conversion_ms = elapsedMilliseconds(start);

    start = std::chrono::steady_clock::now();
    pipeline.original_nef = NefPolyhedron(pipeline.input_polyhedron);
    if(pipeline.original_nef.is_empty() ||
       !pipeline.original_nef.is_valid())
      throw std::runtime_error("Failed to construct a valid bounded Nef");
    pipeline.timings.nef_construction_ms = elapsedMilliseconds(start);

    NefPolyhedron decomposed_nef = pipeline.original_nef;
    start = std::chrono::steady_clock::now();
    CGAL::convex_decomposition_3(decomposed_nef);
    pipeline.timings.convex_decomposition_ms =
        elapsedMilliseconds(start);

    start = std::chrono::steady_clock::now();
    extractConvexParts(pipeline, decomposed_nef);
    extractHalfspaces(pipeline);
    pipeline.timings.halfspace_extraction_ms =
        elapsedMilliseconds(start);

    start = std::chrono::steady_clock::now();
    constructTrees(pipeline);
    pipeline.timings.tree_construction_ms = elapsedMilliseconds(start);

    start = std::chrono::steady_clock::now();
    pipeline.simplified_tree =
        HS::simplifyToFixedPoint(pipeline.raw_tree);
    pipeline.timings.simplification_ms = elapsedMilliseconds(start);

    start = std::chrono::steady_clock::now();
    validatePipeline(pipeline);
    pipeline.timings.validation_ms = elapsedMilliseconds(start);
    pipeline.timings.total_ms = elapsedMilliseconds(total_start);
    writeOutputs(pipeline);

    const auto raw_serialized = HS::serializeTree(pipeline.raw_tree);
    const auto simplified_serialized =
        HS::serializeTree(pipeline.simplified_tree);
    std::cout << std::setprecision(17)
              << "CGAL version: " << kCgalVersion << "\n"
              << "Length unit: " << rokae_demo::kLengthUnit << "\n"
              << "Kernel: Exact_predicates_exact_constructions_kernel\n"
              << "Decomposition: CGAL::convex_decomposition_3\n"
              << "Input faces: "
              << pipeline.input_polyhedron.size_of_facets() << "\n"
              << "Reflex edges: " << pipeline.reflex_edge_count << "\n"
              << "Convex parts: " << pipeline.parts.size() << "\n"
              << "Inserted cut planes: "
              << pipeline.inserted_cut_plane_count << "\n"
              << "Original boundary planes: "
              << pipeline.original_boundary_plane_count << "\n"
              << "Unique leaves: " << pipeline.leaves.size() << "\n"
              << "Raw nodes: " << raw_serialized.nodes.size() << "\n"
              << "Simplified nodes: "
              << simplified_serialized.nodes.size() << "\n"
              << "Raw/simplified max error: "
              << pipeline.validation.maximum_expression_error << "\n"
              << "False negatives: "
              << pipeline.validation.confusion.mesh_inside_tree_outside
              << "\n"
              << "False positives: "
              << pipeline.validation.confusion.mesh_outside_tree_inside
              << "\n"
              << "Original volume (m^3): "
              << pipeline.validation.original_volume << "\n"
              << "Part volume sum (m^3): "
              << pipeline.validation.convex_parts_volume_sum << "\n"
              << "Volume error (m^3): "
              << pipeline.validation.volume_absolute_error << "\n"
              << "Exact symmetric difference empty: "
              << (pipeline.validation.exact_symmetric_difference_empty
                      ? "yes"
                      : "no")
              << "\n"
              << "Total time(ms): " << pipeline.timings.total_ms << "\n"
              << "Output: " << pipeline.options.output_directory << "\n";

    if(!pipelinePassed(pipeline))
    {
      for(const std::string& failure : pipeline.validation.failures)
        std::cerr << "Validation failure: " << failure << "\n";
      if(pipeline.options.expected_parts >= 0 &&
         static_cast<int>(pipeline.parts.size()) !=
             pipeline.options.expected_parts)
        std::cerr << "Expected a different convex-part count.\n";
      if(pipeline.options.expected_planes >= 0 &&
         static_cast<int>(pipeline.leaves.size()) !=
             pipeline.options.expected_planes)
        std::cerr << "Expected a different unique-plane count.\n";
      if(pipeline.options.expected_min_parts >= 0 &&
         static_cast<int>(pipeline.parts.size()) <
             pipeline.options.expected_min_parts)
        std::cerr << "Expected more than one convex part.\n";
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }
  catch(const std::exception& error)
  {
    std::cerr << "build_halfspace_tree failed: " << error.what() << "\n";
    return EXIT_FAILURE;
  }
}
