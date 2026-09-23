#include <rokae_demo/convex_cluster.hpp>
#include <rokae_demo/units.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace rokae_demo
{
namespace fs = std::filesystem;
namespace
{
HS::Vec3 add(const HS::Vec3& lhs, const HS::Vec3& rhs)
{
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

HS::Vec3 scale(const HS::Vec3& value, double factor)
{
  return {value.x * factor, value.y * factor, value.z * factor};
}

double dot(const HS::Vec3& lhs, const HS::Vec3& rhs)
{
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

bool samePlane(const Plane& lhs,
               const Plane& rhs,
               double cosine_tolerance,
               double distance_tolerance)
{
  return dot(lhs.normal, rhs.normal) >= cosine_tolerance &&
         std::abs(lhs.offset - rhs.offset) <= distance_tolerance;
}

void writeIntegerArray(std::ostream& out, const std::vector<int>& values)
{
  out << '[';
  for(std::size_t i = 0; i < values.size(); ++i)
  {
    if(i != 0)
      out << ", ";
    out << values[i];
  }
  out << ']';
}

HS::Vec3 vertexAverage(const ConvexCluster& cluster,
                       const std::vector<HS::Vec3>& vertices)
{
  HS::Vec3 result;
  for(const std::size_t vertex : cluster.vertices)
    result = add(result, vertices.at(vertex));
  return scale(result, 1.0 / static_cast<double>(cluster.vertices.size()));
}

}  // namespace

int findOrAddLeaf(const Plane& plane,
                  int cluster_id,
                  double distance_tolerance,
                  std::vector<HS::HalfspaceLeaf>& leaves,
                  const std::string& source_type,
                  bool exact_normal_match)
{
  for(HS::HalfspaceLeaf& leaf : leaves)
  {
    Plane existing;
    existing.normal = leaf.normal;
    existing.offset = leaf.offset;
    const bool normals_match = !exact_normal_match ||
        (plane.normal.x == leaf.normal.x && plane.normal.y == leaf.normal.y && plane.normal.z == leaf.normal.z);
    if(normals_match && samePlane(plane, existing, 1.0 - 1e-12, distance_tolerance))
    {
      if(std::find(leaf.convex_part_ids.begin(),
                   leaf.convex_part_ids.end(), cluster_id) ==
         leaf.convex_part_ids.end())
        leaf.convex_part_ids.push_back(cluster_id);
      leaf.source_face_ids.insert(leaf.source_face_ids.end(),
                                  plane.source_face_ids.begin(),
                                  plane.source_face_ids.end());
      std::sort(leaf.source_face_ids.begin(), leaf.source_face_ids.end());
      leaf.source_face_ids.erase(
          std::unique(leaf.source_face_ids.begin(),
                      leaf.source_face_ids.end()),
          leaf.source_face_ids.end());
      return leaf.id;
    }
  }
  HS::HalfspaceLeaf leaf;
  leaf.id = static_cast<int>(leaves.size());
  leaf.normal = plane.normal;
  leaf.point = plane.point;
  leaf.offset = plane.offset;
  leaf.convex_part_ids.push_back(cluster_id);
  leaf.source_face_ids = plane.source_face_ids;
  leaf.source_type = source_type;
  leaves.push_back(std::move(leaf));
  return leaves.back().id;
}

TreeRepresentation buildClusterTree(std::vector<ConvexCluster>& clusters,
                                    double plane_tolerance,
                                    const std::string& source_type,
                                    bool exact_normal_match)
{
  TreeRepresentation tree;
  std::vector<HS::ExpressionPtr> cluster_expressions;
  for(ConvexCluster& cluster : clusters)
  {
    cluster.leaf_ids.clear();
    std::vector<HS::ExpressionPtr> halfspaces;
    for(const Plane& plane : cluster.reduced_planes)
    {
      const int leaf_id = findOrAddLeaf(plane, cluster.id,
                                        plane_tolerance, tree.leaves, source_type, exact_normal_match);
      cluster.leaf_ids.push_back(leaf_id);
      halfspaces.push_back(HS::makeLeaf(leaf_id));
      ++tree.leaf_references;
    }
    if(halfspaces.empty())
      throw std::runtime_error("A convex cluster has no support planes");
    cluster_expressions.push_back(
        halfspaces.size() == 1
            ? halfspaces.front()
            : HS::makeNode(HS::TreeOperator::MAX, std::move(halfspaces)));
  }
  if(cluster_expressions.empty())
    throw std::runtime_error("No convex clusters for expression tree");
  tree.raw = cluster_expressions.size() == 1
                 ? cluster_expressions.front()
                 : HS::makeNode(HS::TreeOperator::MIN,
                                std::move(cluster_expressions));
  return tree;
}

bool insidePlanes(const std::vector<Plane>& planes,
                  const HS::Vec3& point,
                  double epsilon)
{
  for(const Plane& plane : planes)
  {
    if(plane.evaluate(point) > epsilon)
      return false;
  }
  return true;
}

bool insideClusterUnion(const std::vector<ConvexCluster>& clusters,
                        const HS::Vec3& point,
                        double epsilon)
{
  for(const ConvexCluster& cluster : clusters)
  {
    if(insidePlanes(cluster.exact_planes, point, epsilon))
      return true;
  }
  return false;
}

void validateClusters(std::vector<ConvexCluster>& clusters,
                      const std::vector<HS::Vec3>& vertices,
                      double epsilon,
                      ClusterValidation& validation)
{
  for(ConvexCluster& cluster : clusters)
  {
    cluster.centroid = vertexAverage(cluster, vertices);
    for(const Plane& plane : cluster.exact_planes)
    {
      for(const std::size_t vertex : cluster.vertices)
      {
        if(plane.evaluate(vertices.at(vertex)) > epsilon)
        {
          validation.all_clusters_convex = false;
          validation.failures.push_back(
              "A cluster boundary plane is not supporting");
          break;
        }
      }
    }
    for(const Plane& plane : cluster.reduced_planes)
    {
      for(const std::size_t vertex : cluster.vertices)
      {
        if(plane.evaluate(vertices.at(vertex)) > epsilon)
        {
          validation.all_cluster_vertices_satisfy = false;
          validation.failures.push_back(
              "A cluster vertex violates a reduced halfspace");
          break;
        }
      }
      if(plane.evaluate(cluster.centroid) >= -epsilon)
      {
        validation.all_cluster_centroids_strictly_inside = false;
        validation.failures.push_back(
            "A cluster centroid is not strictly inside its halfspaces");
      }
    }
  }
}

void writeTreeJson(const fs::path& path,
                   const TreeRepresentation& representation,
                   const HS::ExpressionPtr& expression,
                   const std::string& kind,
                   const std::string& pipeline)
{
  std::ofstream out(path);
  if(!out)
    throw std::runtime_error("Failed to write tree JSON");
  out.exceptions(std::ios::badbit | std::ios::failbit);
  const HS::SerializedTree tree = HS::serializeTree(expression);
  out << std::setprecision(17)
      << "{\n  \"pipeline\": \"" << pipeline << "\",\n"
      << "  \"length_unit\": \"" << kLengthUnit << "\",\n"
      << "  \"halfspace_offset_unit\": \"" << kLengthUnit << "\",\n"
      << "  \"tree_kind\": \"" << kind << "\",\n"
      << "  \"sign_convention\": \"psi(x)=n^T x-offset; inside iff root<=0\",\n"
      << "  \"leaves\": [\n";
  for(std::size_t i = 0; i < representation.leaves.size(); ++i)
  {
    const HS::HalfspaceLeaf& leaf = representation.leaves[i];
    out << "    {\"id\": " << leaf.id << ", \"normal\": ["
        << leaf.normal.x << ", " << leaf.normal.y << ", "
        << leaf.normal.z << "], \"point\": [" << leaf.point.x << ", "
        << leaf.point.y << ", " << leaf.point.z << "], \"offset\": "
        << leaf.offset << ", \"convex_part_ids\": ";
    writeIntegerArray(out, leaf.convex_part_ids);
    out << ", \"source_face_ids\": ";
    writeIntegerArray(out, leaf.source_face_ids);
    out << "}" << (i + 1 == representation.leaves.size() ? "\n" : ",\n");
  }
  out << "  ],\n  \"nodes\": [\n";
  for(std::size_t i = 0; i < tree.nodes.size(); ++i)
  {
    const HS::SerializedNode& node = tree.nodes[i];
    out << "    {\"id\": " << node.id << ", \"type\": \""
        << HS::operatorName(node.op) << "\", \"leaf_id\": "
        << node.leaf_id << ", \"children\": ";
    writeIntegerArray(out, node.children);
    out << "}" << (i + 1 == tree.nodes.size() ? "\n" : ",\n");
  }
  out << "  ],\n  \"root_id\": " << tree.root_id << "\n}\n";
  out.close();
}

}  // namespace rokae_demo
