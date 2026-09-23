#include <rokae_demo/octree_convex_approximation.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/boost/graph/iterator.h>

#include <numeric>
#include <stdexcept>

namespace rokae_demo
{
namespace
{
using K = CGAL::Exact_predicates_exact_constructions_kernel;
using FT = K::FT;
using Point = K::Point_3;
using ExactPlane = K::Plane_3;
using Mesh = CGAL::Surface_mesh<Point>;
using Box = std::array<std::array<FT,3>,2>;

void require(bool condition, const char* message)
{
  if(!condition) throw std::runtime_error(message);
}

void validateOptions(const OctreeConvexOptions& options)
{
  require(std::isfinite(options.max_volume_inflation) && options.max_volume_inflation >= 0 &&
          std::isfinite(options.max_fill_distance) && options.max_fill_distance >= 0 &&
          options.max_candidate_voxels > 0 && options.max_certificate_boxes > 0, "Invalid convex approximation budgets");
}

Point exact(const HS::Vec3& p) { return {p.x,p.y,p.z}; }
// Lazy exact intervals can round differently for algebraically equal volumes
// constructed through different meshes. Canonicalize only reported values;
// acceptance decisions below remain exact rational comparisons.
double rounded(const FT& value) { return CGAL::to_double(CGAL::exact(value)); }
HS::Vec3 approximate(const Point& p) { return {rounded(p.x()),rounded(p.y()),rounded(p.z())}; }
ExactPlane exact(const Plane& p) { return {p.normal.x,p.normal.y,p.normal.z,-FT(p.offset)}; }

Box exactBox(const OctreeBox& box)
{
  Box result;
  for(unsigned k = 0; k < 3; ++k) { result[0][k] = FT(box.lower[k]); result[1][k] = FT(box.upper[k]); }
  return result;
}

FT boxVolume(const OctreeBox& box)
{
  FT result(1);
  for(unsigned k = 0; k < 3; ++k) result *= FT(box.upper[k])-FT(box.lower[k]);
  return result;
}

FT meshVolume(const Mesh& mesh)
{
  require(mesh.number_of_vertices() >= 4 && mesh.number_of_faces() >= 4, "Degenerate convex polyhedron");
  const Point origin = mesh.point(*mesh.vertices().begin());
  FT result(0);
  for(const auto f : mesh.faces())
  {
    std::vector<Point> polygon;
    for(const auto v : CGAL::vertices_around_face(mesh.halfedge(f), mesh)) polygon.push_back(mesh.point(v));
    for(std::size_t i = 1; i+1 < polygon.size(); ++i) result += CGAL::volume(origin,polygon[0],polygon[i],polygon[i+1]);
  }
  require(result > 0, "Convex polyhedron is not outward oriented/full dimensional");
  return result;
}

std::vector<ExactPlane> exactPlanes(const std::vector<Plane>& planes)
{
  std::vector<ExactPlane> result;
  for(const auto& p : planes)
  {
    require(std::isfinite(p.normal.x) && std::isfinite(p.normal.y) && std::isfinite(p.normal.z) &&
            std::isfinite(p.offset), "Nonfinite halfspace coefficient");
    require(p.normal.x != 0 || p.normal.y != 0 || p.normal.z != 0, "Zero plane normal");
    result.push_back(exact(p));
  }
  return result;
}

Mesh intersection(const std::vector<Plane>& planes, const Point& interior)
{
  const auto source = exactPlanes(planes);
  for(const auto& p : source) require(p.has_on_negative_side(interior), "Halfspace intersection lacks strict interior point");
  Mesh mesh;
  CGAL::halfspace_intersection_3(source.begin(), source.end(), mesh, boost::make_optional(interior));
  return mesh;
}

Point center(const std::set<std::size_t>& ids, const std::vector<HS::Vec3>& vertices)
{
  require(!ids.empty(), "Empty convex source vertex set");
  FT x(0),y(0),z(0);
  for(const auto id : ids) { const auto& p = vertices.at(id); x += FT(p.x); y += FT(p.y); z += FT(p.z); }
  const FT count(static_cast<unsigned long>(ids.size()));
  return {x/count,y/count,z/count};
}

void displayMesh(const Mesh& mesh, OctreeConvexCell& cell)
{
  std::map<Mesh::Vertex_index,std::size_t> ids;
  for(const auto v : mesh.vertices()) { ids.emplace(v,cell.mesh_vertices.size()); cell.mesh_vertices.push_back(approximate(mesh.point(v))); }
  for(const auto f : mesh.faces())
  {
    std::vector<std::size_t> polygon;
    for(const auto v : CGAL::vertices_around_face(mesh.halfedge(f),mesh)) polygon.push_back(ids.at(v));
    for(std::size_t i = 1; i+1 < polygon.size(); ++i) cell.mesh_triangles.push_back({polygon[0],polygon[i],polygon[i+1]});
  }
}

void displayBox(const OctreeBox& box, OctreeConvexCell& cell)
{
  for(unsigned mask = 0; mask < 8; ++mask)
    cell.mesh_vertices.push_back({(mask&1) ? box.upper[0] : box.lower[0],
                                  (mask&2) ? box.upper[1] : box.lower[1],
                                  (mask&4) ? box.upper[2] : box.lower[2]});
  for(unsigned axis = 0; axis < 3; ++axis) for(unsigned side = 0; side < 2; ++side)
  {
    const unsigned base = side << axis, u = 1u << ((axis+1)%3), v = 1u << ((axis+2)%3);
    std::array<std::size_t,4> q{base,base|u,base|u|v,base|v};
    if(!side) std::swap(q[1],q[3]);
    cell.mesh_triangles.push_back({q[0],q[1],q[2]});
    cell.mesh_triangles.push_back({q[0],q[2],q[3]});
  }
}

std::vector<Plane> activePlanes(const std::vector<Plane>& planes, const Mesh& mesh)
{
  std::vector<Plane> result;
  for(const auto& plane : planes)
  {
    const auto p = exact(plane);
    std::vector<Point> on;
    for(const auto v : mesh.vertices()) if(p.has_on(mesh.point(v))) on.push_back(mesh.point(v));
    bool facet = false;
    for(std::size_t i = 2; i < on.size() && !facet; ++i) facet = !CGAL::collinear(on[0],on[1],on[i]);
    if(facet && std::none_of(result.begin(),result.end(),[&](const Plane& other) { return exact(other) == p; })) result.push_back(plane);
  }
  require(result.size() >= 4, "Lost essential convex support planes");
  return result;
}

// Each imported double normal defines a DIFFERENT rational plane from the
// ideal hull plane. Recompute its support over all source corners, outwardly
// round the offset and include a conservative floating dot-product guard.
Plane support(const ExactPlane& facet, const std::set<std::size_t>& ids,
              const std::vector<HS::Vec3>& vertices, int face_id)
{
  const FT scale = (std::max)({CGAL::abs(facet.a()),CGAL::abs(facet.b()),CGAL::abs(facet.c())});
  double a = rounded(facet.a()/scale), b = rounded(facet.b()/scale), c = rounded(facet.c()/scale);
  const double length = std::hypot(a,b,c);
  Plane result;
  result.normal = {a/length,b/length,c/length};
  const FT nx(result.normal.x), ny(result.normal.y), nz(result.normal.z);
  FT maximum = nx*FT(vertices.at(*ids.begin()).x) + ny*FT(vertices.at(*ids.begin()).y) + nz*FT(vertices.at(*ids.begin()).z);
  FT magnitude(0);
  for(const auto id : ids)
  {
    const auto& p = vertices.at(id);
    const FT x = nx*FT(p.x), y = ny*FT(p.y), z = nz*FT(p.z);
    maximum = (std::max)(maximum,FT(x+y+z));
    magnitude = (std::max)(magnitude,FT(CGAL::abs(x)+CGAL::abs(y)+CGAL::abs(z)));
  }
  const FT guard = FT(64)*FT(std::numeric_limits<double>::epsilon())*magnitude + FT(16)*FT(std::numeric_limits<double>::denorm_min());
  result.offset = std::nextafter(CGAL::to_interval(maximum+guard).second, std::numeric_limits<double>::infinity());
  require(std::isfinite(result.offset), "Outward plane support is not representable");
  result.point = {result.normal.x*result.offset,result.normal.y*result.offset,result.normal.z*result.offset};
  result.source_face_ids = {face_id};
  return result;
}

FT evaluate(const ExactPlane& p, const std::array<FT,3>& x)
{
  return p.a()*x[0]+p.b()*x[1]+p.c()*x[2]+p.d();
}

bool outside(const Box& box, const std::vector<ExactPlane>& planes)
{
  for(const auto& p : planes)
  {
    const std::array<FT,3> n{p.a(),p.b(),p.c()};
    std::array<FT,3> corner;
    for(unsigned k = 0; k < 3; ++k) corner[k] = box[n[k] < 0 ? 1 : 0][k];
    if(evaluate(p,corner) > 0) return true;
  }
  return false;
}

struct DistanceCertificate { bool passed = false, capped = false; FT bound2 = FT(0); std::size_t boxes = 0; };

DistanceCertificate distanceCertificate(const OctreeBox& domain, const std::vector<Plane>& supports,
                                        const std::vector<OctreeBox>& occupied, double budget, std::size_t limit)
{
  DistanceCertificate result;
  const FT budget2 = FT(budget)*FT(budget);
  const auto planes = exactPlanes(supports);
  std::vector<Box> sources;
  for(const auto& box : occupied) sources.push_back(exactBox(box));
  require(!sources.empty(), "Missing occupied boxes for distance certificate");
  std::vector<Box> pending{exactBox(domain)};
  while(!pending.empty())
  {
    if(result.boxes == limit) { result.capped = true; return result; }
    ++result.boxes;
    const auto box = pending.back(); pending.pop_back();
    if(outside(box,planes)) continue;
    std::array<FT,3> middle;
    for(unsigned k = 0; k < 3; ++k) middle[k] = (box[0][k]+box[1][k])/FT(2);
    FT upper(0), nearest(0);
    bool first = true;
    for(const auto& source : sources)
    {
      FT upper2(0), distance2(0);
      for(unsigned k = 0; k < 3; ++k)
      {
        const FT far = (std::max)({FT(0),FT(source[0][k]-box[0][k]),FT(box[1][k]-source[1][k])});
        const FT near = (std::max)({FT(0),FT(source[0][k]-middle[k]),FT(middle[k]-source[1][k])});
        upper2 += far*far; distance2 += near*near;
      }
      if(first || upper2 < upper) upper = upper2;
      if(first || distance2 < nearest) nearest = distance2;
      first = false;
    }
    if(upper <= budget2) { result.bound2 = (std::max)(result.bound2,upper); continue; }
    if(nearest > budget2 && std::all_of(planes.begin(),planes.end(),[&](const ExactPlane& p) { return evaluate(p,middle) <= 0; }))
      return result; // exact inside witness exceeds the requested distance
    unsigned axis = 0;
    for(unsigned k = 1; k < 3; ++k) if(box[1][k]-box[0][k] > box[1][axis]-box[0][axis]) axis = k;
    Box lower = box, higher = box;
    lower[1][axis] = higher[0][axis] = middle[axis];
    pending.push_back(std::move(higher)); pending.push_back(std::move(lower));
  }
  result.passed = true;
  return result;
}

double distanceBound(const FT& squared, double budget)
{
  if(squared == 0) return 0;
  return (std::min)(budget,std::nextafter(std::sqrt(CGAL::to_interval(CGAL::exact(squared)).second),std::numeric_limits<double>::infinity()));
}

struct Sources
{
  std::vector<HS::Vec3> vertices;
  std::vector<std::array<std::size_t,8>> corners;
};

Sources sourceVertices(const OctreeOccupancy& octree)
{
  Sources result;
  std::map<VoxelIndex,std::size_t> ids;
  for(const auto& voxel : octree.occupied_voxels)
  {
    std::array<std::size_t,8> corners{};
    for(unsigned mask = 0; mask < 8; ++mask)
    {
      auto key = voxel;
      for(unsigned k = 0; k < 3; ++k) key[k] += (mask >> k)&1u;
      const auto entry = ids.emplace(key,result.vertices.size());
      if(entry.second) result.vertices.push_back({static_cast<double>(key[0])*octree.voxel_size,
                                                  static_cast<double>(key[1])*octree.voxel_size,
                                                  static_cast<double>(key[2])*octree.voxel_size});
      corners[mask] = entry.first->second;
    }
    result.corners.push_back(corners);
  }
  return result;
}

std::set<std::size_t> cornerSet(const std::vector<std::size_t>& voxels, const Sources& source)
{
  std::set<std::size_t> result;
  for(const auto id : voxels) for(const auto corner : source.corners.at(id)) result.insert(corner);
  return result;
}

OctreeBox tightBox(const std::set<std::size_t>& ids, const std::vector<HS::Vec3>& source)
{
  const auto& p = source.at(*ids.begin());
  OctreeBox box{{p.x,p.y,p.z},{p.x,p.y,p.z}};
  for(const auto id : ids)
  {
    const auto& v = source.at(id);
    const std::array<double,3> coords{v.x,v.y,v.z};
    for(unsigned k = 0; k < 3; ++k) { box.lower[k] = (std::min)(box.lower[k],coords[k]); box.upper[k] = (std::max)(box.upper[k],coords[k]); }
  }
  return box;
}

std::vector<Plane> boxPlanes(const OctreeBox& box)
{
  std::vector<HS::Vec3> vertices;
  return aabbClusters({box},{0},vertices).front().exact_planes;
}

bool sameCoefficients(const Plane& a, const Plane& b)
{
  return a.normal.x == b.normal.x && a.normal.y == b.normal.y && a.normal.z == b.normal.z && a.offset == b.offset;
}

bool inNode(const VoxelIndex& key, const OctreeNode& node)
{
  for(unsigned k = 0; k < 3; ++k) if(key[k] < node.lower[k] || key[k] >= node.lower[k]+node.width) return false;
  return true;
}

std::vector<OctreeBox> originalBoxes(const std::vector<std::size_t>& voxels, const OctreeCellPartition& baseline)
{
  std::set<std::size_t> ids;
  for(const auto voxel : voxels) ids.insert(baseline.voxel_cell_ids.at(voxel));
  std::vector<OctreeBox> result;
  for(const auto id : ids) result.push_back(baseline.cells.at(id).box);
  return result;
}

FT sourceVolume(const std::vector<std::size_t>& voxels, const OctreeOccupancy& octree)
{
  FT result(0);
  for(const auto voxel : voxels) result += boxVolume(octree.nodes.at(octree.leaf_nodes.at(voxel)).box);
  return result;
}
}  // namespace

OctreeConvexResult approximateOctreeConvex(const OctreeOccupancy& octree, const OctreeCellPartition& baseline,
                                          const OctreeConvexOptions& options)
{
  validateOptions(options);
  validateOctreeCellPartition(octree,baseline);
  OctreeConvexResult result;
  const auto source = sourceVertices(octree);
  result.source_vertices = source.vertices;
  struct Record { ConvexCluster cluster; OctreeConvexCell cell; FT volume; };
  std::vector<Record> pool;
  std::vector<HS::Vec3> ignored;
  auto original = octreeAabbClusters(baseline,ignored);
  std::vector<std::vector<std::size_t>> owned(octree.nodes.size());
  FT original_volume(0);
  for(std::size_t i = 0; i < original.size(); ++i)
  {
    auto& cluster = original[i];
    cluster.vertices = cornerSet(cluster.cells,source);
    OctreeConvexCell cell;
    cell.baseline_cell = static_cast<int>(i);
    cell.source_node = baseline.cells[i].source_node >= 0 ? baseline.cells[i].source_node : baseline.cells[i].parent_node;
    const FT volume = boxVolume(baseline.cells[i].box);
    cell.volume = rounded(volume);
    owned.at(cell.source_node).push_back(i);
    original_volume += volume;
    pool.push_back({std::move(cluster),std::move(cell),volume});
  }
  FT current_volume = original_volume;
  const FT ratio(options.max_volume_inflation);
  int face_id = static_cast<int>(6*original.size());
  struct Subtree { std::vector<std::size_t> voxels, active; };
  std::vector<int> voxel_at_node(octree.nodes.size(),-1);
  for(std::size_t i = 0; i < octree.leaf_nodes.size(); ++i) voxel_at_node[octree.leaf_nodes[i]] = static_cast<int>(i);
  const auto visit = [&](const auto& self, int node_id) -> Subtree {
    const auto& node = octree.nodes.at(node_id);
    Subtree data;
    data.active = owned[node_id];
    if(voxel_at_node[node_id] >= 0) data.voxels.push_back(static_cast<std::size_t>(voxel_at_node[node_id]));
    for(const int child : node.children) if(child >= 0)
    {
      auto sub = self(self,child);
      data.voxels.insert(data.voxels.end(),sub.voxels.begin(),sub.voxels.end());
      data.active.insert(data.active.end(),sub.active.begin(),sub.active.end());
    }
    if(data.active.size() < 2) return data;
    ++result.candidate_nodes;
    if(data.voxels.size() > options.max_candidate_voxels) { ++result.rejected_work_limit; return data; }
    ++result.candidate_hulls;
    std::sort(data.voxels.begin(),data.voxels.end());
    ConvexCluster cluster;
    cluster.cells = data.voxels;
    cluster.vertices = cornerSet(data.voxels,source);
    std::vector<Point> points;
    for(const auto id : cluster.vertices) points.push_back(exact(source.vertices[id]));
    const auto interior = center(cluster.vertices,source.vertices);
    Mesh hull;
    CGAL::convex_hull_3(points.begin(),points.end(),hull);
    const FT source_volume = sourceVolume(data.voxels,octree);
    FT replaced_volume(0);
    std::size_t replaced_planes = 0;
    for(const auto id : data.active) { replaced_volume += pool[id].volume; replaced_planes += pool[id].cluster.reduced_planes.size(); }
    // The ideal hull is a lower bound on the supported output polyhedron.
    // Reject cheaply here, but still measure the rounded output before accept.
    const FT minimum_volume = meshVolume(hull);
    if(minimum_volume-source_volume > ratio*source_volume ||
       current_volume-replaced_volume+minimum_volume-original_volume > ratio*original_volume)
    { ++result.rejected_volume; return data; }
    std::vector<ExactPlane> facets;
    for(const auto face : hull.faces())
    {
      const auto h = hull.halfedge(face);
      ExactPlane plane(hull.point(hull.target(h)),hull.point(hull.target(hull.next(h))),hull.point(hull.target(hull.next(hull.next(h)))));
      if(plane.has_on_positive_side(interior)) plane = plane.opposite();
      if(std::find(facets.begin(),facets.end(),plane) == facets.end()) facets.push_back(plane);
    }
    const auto domain = tightBox(cluster.vertices,source.vertices);
    // Exact axis clamps keep the final polyhedron inside its node and make
    // distinct selected hierarchy regions interior-disjoint even after padding.
    cluster.exact_planes = boxPlanes(domain);
    for(auto& plane : cluster.exact_planes) plane.source_face_ids = {face_id++};
    for(const auto& facet : facets) cluster.exact_planes.push_back(support(facet,cluster.vertices,source.vertices,face_id++));
    Mesh poly = intersection(cluster.exact_planes,interior);
    const FT volume = meshVolume(poly);
    require(volume >= source_volume, "Candidate lost occupied volume");
    if(volume-source_volume > ratio*source_volume || current_volume-replaced_volume+volume-original_volume > ratio*original_volume)
    { ++result.rejected_volume; return data; }
    cluster.reduced_planes = activePlanes(cluster.exact_planes,poly);
    if(cluster.reduced_planes.size() >= replaced_planes) { ++result.rejected_gain; return data; }
    DistanceCertificate certificate;
    if(volume == source_volume) certificate.passed = true; // regular-closed sets + containment + equal volume => equality
    else certificate = distanceCertificate(domain,cluster.reduced_planes,originalBoxes(data.voxels,baseline),
                                            options.max_fill_distance,options.max_certificate_boxes);
    result.certificate_boxes += certificate.boxes;
    if(!certificate.passed)
    {
      if(certificate.capped) ++result.rejected_work_limit;
      else ++result.rejected_distance;
      return data;
    }
    OctreeConvexCell cell;
    cell.source_node = node_id;
    cell.volume = rounded(volume);
    cell.fill_distance_bound = distanceBound(certificate.bound2,options.max_fill_distance);
    for(const auto& p : cluster.reduced_planes)
      if((p.normal.x != 0)+(p.normal.y != 0)+(p.normal.z != 0) > 1) cell.non_aabb = true;
    displayMesh(poly,cell);
    current_volume += volume-replaced_volume;
    data.active = {pool.size()};
    pool.push_back({std::move(cluster),std::move(cell),volume});
    ++result.accepted_merges;
    return data;
  };
  auto final = visit(visit,0).active;
  std::sort(final.begin(),final.end(),[&](auto a, auto b) { return pool[a].cluster.cells.front() < pool[b].cluster.cells.front(); });
  result.voxel_cluster_ids.resize(octree.occupied_voxels.size());
  for(const auto id : final)
  {
    auto& record = pool[id];
    const auto cluster_id = result.clusters.size();
    record.cluster.id = static_cast<int>(cluster_id);
    if(record.cell.mesh_vertices.empty()) displayBox(baseline.cells.at(record.cell.baseline_cell).box,record.cell);
    for(const auto voxel : record.cluster.cells) result.voxel_cluster_ids.at(voxel) = cluster_id;
    result.maximum_fill_distance_bound = (std::max)(result.maximum_fill_distance_bound,record.cell.fill_distance_bound);
    result.clusters.push_back(std::move(record.cluster));
    result.cells.push_back(std::move(record.cell));
  }
  result.occupied_volume = rounded(original_volume);
  result.output_volume = rounded(current_volume);
  result.volume_inflation = rounded((current_volume-original_volume)/original_volume);
  return result;
}

void validateOctreeConvex(const OctreeOccupancy& octree, const OctreeCellPartition& baseline,
                         const OctreeConvexOptions& options, const OctreeConvexResult& result)
{
  validateOptions(options);
  validateOctreeCellPartition(octree,baseline);
  const auto source = sourceVertices(octree);
  require(result.source_vertices.size() == source.vertices.size() && result.clusters.size() == result.cells.size() &&
          !result.clusters.empty() && result.voxel_cluster_ids.size() == octree.occupied_voxels.size(), "Convex result sizes inconsistent");
  for(std::size_t i = 0; i < source.vertices.size(); ++i)
  {
    const auto& a = result.source_vertices[i]; const auto& b = source.vertices[i];
    require(a.x == b.x && a.y == b.y && a.z == b.z, "Convex source vertices changed");
  }
  std::vector<unsigned char> covered(octree.occupied_voxels.size(),0);
  FT volume_sum(0), source_sum(0);
  double maximum_distance = 0;
  std::vector<HS::Vec3> ignored;
  const auto original = octreeAabbClusters(baseline,ignored);
  for(std::size_t i = 0; i < result.clusters.size(); ++i)
  {
    const auto& cluster = result.clusters[i]; const auto& cell = result.cells[i];
    require(cluster.id == static_cast<int>(i) && !cluster.cells.empty() &&
            cluster.vertices == cornerSet(cluster.cells,source), "Invalid convex source cluster");
    for(const auto voxel : cluster.cells)
    {
      require(voxel < covered.size() && !covered[voxel] && result.voxel_cluster_ids[voxel] == i, "Duplicated/missing convex voxel owner");
      covered[voxel] = 1;
    }
    const FT source_volume = sourceVolume(cluster.cells,octree);
    const auto interior = center(cluster.vertices,source.vertices);
    require(cluster.exact_planes.size() >= cluster.reduced_planes.size(), "Invalid reduced support count");
    for(const auto& plane : cluster.reduced_planes)
      require(std::any_of(cluster.exact_planes.begin(),cluster.exact_planes.end(),[&](const Plane& p) { return sameCoefficients(p,plane); }),
              "Reduction introduced or moved a plane");
    for(const auto id : cluster.vertices)
      for(const auto& p : cluster.exact_planes)
        require(!exact(p).has_on_positive_side(exact(source.vertices[id])) && p.evaluate(source.vertices[id]) <= 0,
                "Support plane misses an occupied corner");
    const auto poly = intersection(cluster.reduced_planes,interior);
    const FT volume = meshVolume(poly);
    const auto domain = exactBox(tightBox(cluster.vertices,source.vertices));
    for(const auto v : poly.vertices()) for(unsigned k = 0; k < 3; ++k)
      require(poly.point(v)[k] >= domain[0][k] && poly.point(v)[k] <= domain[1][k], "Polyhedron escapes the distance certificate domain");
    // Check removed constraints on ALL vertices of the reduced intersection:
    // reduction is exact and did not silently enlarge the candidate.
    for(const auto v : poly.vertices()) for(const auto& p : cluster.exact_planes)
      require(!exact(p).has_on_positive_side(poly.point(v)), "Plane reduction changed polyhedron");
    if(cell.baseline_cell >= 0)
    {
      const auto& old = original.at(cell.baseline_cell);
      const auto& old_cell = baseline.cells.at(cell.baseline_cell);
      require(cell.source_node == (old_cell.source_node >= 0 ? old_cell.source_node : old_cell.parent_node), "Retained cell owner changed");
      require(cluster.cells == old.cells && cluster.reduced_planes.size() == 6 && cell.fill_distance_bound == 0 && !cell.non_aabb,
              "Retained exact cell changed");
      for(std::size_t j = 0; j < 6; ++j) require(sameCoefficients(cluster.reduced_planes[j],old.reduced_planes[j]), "Retained AABB plane changed");
    }
    else
    {
      const auto& node = octree.nodes.at(cell.source_node);
      require(node.width > 1 && cluster.cells.size() == node.occupied_count, "Hull owner must cover all occupied voxels in node");
      for(const auto voxel : cluster.cells) require(inNode(octree.occupied_voxels[voxel],node), "Hull source escapes owner node");
      const auto node_box = exactBox(node.box);
      for(const auto v : poly.vertices()) for(unsigned k = 0; k < 3; ++k)
        require(poly.point(v)[k] >= node_box[0][k] && poly.point(v)[k] <= node_box[1][k], "Padded hull escapes node; disjointness lost");
    }
    require(volume >= source_volume && volume-source_volume <= FT(options.max_volume_inflation)*source_volume, "Local volume budget exceeded");
    DistanceCertificate certificate;
    if(volume == source_volume) certificate.passed = true;
    else certificate = distanceCertificate(tightBox(cluster.vertices,source.vertices),cluster.reduced_planes,originalBoxes(cluster.cells,baseline),
                                           options.max_fill_distance,options.max_certificate_boxes);
    require(certificate.passed && cell.fill_distance_bound == distanceBound(certificate.bound2,options.max_fill_distance), "Distance certificate failed");
    require(cell.volume == rounded(volume), "Reported cell volume changed");
    maximum_distance = (std::max)(maximum_distance,cell.fill_distance_bound);
    bool non_aabb = false;
    for(const auto& p : cluster.reduced_planes) non_aabb |= (p.normal.x != 0)+(p.normal.y != 0)+(p.normal.z != 0) > 1;
    require(non_aabb == cell.non_aabb, "Non-AABB metadata inconsistent");
    source_sum += source_volume; volume_sum += volume;
  }
  require(std::all_of(covered.begin(),covered.end(),[](unsigned char v) { return v == 1; }), "Incomplete convex voxel coverage");
  require(volume_sum-source_sum <= FT(options.max_volume_inflation)*source_sum, "Global volume budget exceeded");
  require(result.maximum_fill_distance_bound == maximum_distance, "Reported maximum fill distance inconsistent");
  require(result.occupied_volume == rounded(source_sum) && result.output_volume == rounded(volume_sum) &&
          result.volume_inflation == rounded((volume_sum-source_sum)/source_sum), "Reported convex volumes inconsistent");
}
}  // namespace rokae_demo
