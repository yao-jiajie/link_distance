#pragma once

#include <rokae_demo/octree_rectangular_pruning.hpp>

#include <filesystem>

namespace rokae_demo
{
struct SparseVoxelOccupancy;
// An oriented grid face. u=(axis+1)%3, v=(axis+2)%3, so u cross v
// points along +axis. Indices describe geometry; meters are converted only
// at export using exactly the same rounded grid coordinates as occupancy.
struct OctreeBoundaryFace
{
  unsigned axis = 0;
  int sign = 1;
  std::int64_t plane = 0, u = 0, v = 0;
  std::size_t voxel_id = 0, convex_cell_id = 0;
};

struct OctreeBoundaryPatch
{
  unsigned axis = 0;
  int sign = 1;
  std::int64_t plane = 0, u0 = 0, u1 = 0, v0 = 0, v1 = 0;
  std::vector<std::size_t> source_face_ids;
};

struct OctreeBoundary
{
  std::size_t voxel_faces_before = 0, internal_faces_removed = 0;
  std::vector<OctreeBoundaryFace> faces;
  std::vector<OctreeBoundaryPatch> patches;
  double extraction_ms = 0, merge_ms = 0;
};

// Exact exterior of the occupied union (including cavity walls), NOT a support
// set for the generally nonconvex union. Never remove inter-cluster constraints
// from ConvexCluster or feed all these boundary planes into a single MAX.
// Requires a valid occupancy and a certified partition of that occupancy.
OctreeBoundary extractOctreeBoundary(const OctreeOccupancy&, const OctreeCellPartition&);
// Same geometry without a hierarchy/rect partition; face owner IDs are voxel IDs.
OctreeBoundary extractVoxelBoundary(const std::vector<VoxelIndex>&);
// Sparse-occupancy convenience path; output matches the key-vector overload.
OctreeBoundary extractSparseVoxelBoundary(const SparseVoxelOccupancy&);
void validateVoxelBoundary(const std::vector<VoxelIndex>&, const OctreeBoundary&);
void writeVoxelBoundary(const std::filesystem::path&, double voxel_size, const OctreeBoundary&);

// Checks every expected exposed face and every patch's exact, disjoint tiling.
// Independent of random validation samples. Throws on any mismatch.
void validateOctreeBoundary(const OctreeOccupancy&, const OctreeCellPartition&, const OctreeBoundary&);

std::array<VoxelIndex, 4> octreeBoundaryCorners(const OctreeBoundaryPatch&);

// boundary_patches.json: merged rectangular polygons, normals and provenance.
// boundary.off: outward triangles of the original exposed unit faces; keeping
// grid subdivisions avoids T-junctions introduced by rectangular patch merges.
// Edge/vertex-only voxel contact may still be nonmanifold. No mesh is fed back
// into the tree, and no convex decomposition is invoked.
void writeOctreeBoundary(const std::filesystem::path&, const OctreeOccupancy&, const OctreeBoundary&);
}  // namespace rokae_demo
