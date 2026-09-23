#pragma once
#include <rokae_demo/halfspace_logic_tree.hpp>
#include <array>
#include <filesystem>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace rokae_demo
{
namespace HS = halfspace;

// Shared representation: cells are source tetrahedra or occupied voxels.
struct FaceKey
{
  std::array<std::size_t, 3> vertices{};

  bool operator<(const FaceKey& other) const
  {
    return vertices < other.vertices;
  }
};

struct EdgeKey
{
  std::array<std::size_t, 2> vertices{};

  bool operator<(const EdgeKey& other) const
  {
    return vertices < other.vertices;
  }
};

struct BoundaryFace
{
  FaceKey key;
  std::size_t opposite_vertex = 0;
  int source_face_id = -1;
};

struct Plane
{
  HS::Vec3 normal;
  HS::Vec3 point;
  double offset = 0.0;
  double area = 0.0;
  std::vector<int> source_face_ids;

  double evaluate(const HS::Vec3& value) const
  {
    return normal.x * value.x + normal.y * value.y +
           normal.z * value.z - offset;
  }
};

struct ConvexCluster
{
  int id = -1;
  bool active = true;
  std::size_t version = 0;
  std::vector<std::size_t> cells;
  std::map<FaceKey, BoundaryFace> boundary_faces;
  std::set<std::size_t> vertices;
  std::set<int> neighbors;
  std::vector<Plane> exact_planes;
  std::vector<Plane> reduced_planes;
  std::vector<int> leaf_ids;
  HS::Vec3 centroid;
};

struct TreeRepresentation
{
  std::vector<HS::HalfspaceLeaf> leaves;
  HS::ExpressionPtr raw;
  HS::ExpressionPtr simplified;
  std::size_t leaf_references = 0;
};

struct ClusterValidation
{
  bool all_clusters_convex = true;
  bool all_cluster_vertices_satisfy = true;
  bool all_cluster_centroids_strictly_inside = true;
  std::vector<std::string> failures;
};

int findOrAddLeaf(const Plane& plane, int cluster_id, double distance_tolerance,
                  std::vector<HS::HalfspaceLeaf>& leaves,
                  const std::string& source_type = "alpha_tet_boundary",
                  bool exact_normal_match = false);
TreeRepresentation buildClusterTree(std::vector<ConvexCluster>& clusters,
                                    double plane_tolerance,
                                    const std::string& source_type = "alpha_tet_boundary",
                                    bool exact_normal_match = false);
bool insidePlanes(const std::vector<Plane>& planes, const HS::Vec3& point, double epsilon);
bool insideClusterUnion(const std::vector<ConvexCluster>& clusters,
                        const HS::Vec3& point, double epsilon);
void validateClusters(std::vector<ConvexCluster>& clusters,
                      const std::vector<HS::Vec3>& vertices, double epsilon,
                      ClusterValidation& validation);
void writeTreeJson(const std::filesystem::path& path,
                   const TreeRepresentation& representation,
                   const HS::ExpressionPtr& expression, const std::string& kind,
                   const std::string& pipeline = "alpha-tet");
}  // namespace rokae_demo
