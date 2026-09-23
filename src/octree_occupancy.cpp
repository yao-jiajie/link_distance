#include <rokae_demo/octree_occupancy.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace rokae_demo
{
OctreeBox octreeGridBox(const VoxelIndex& lower, const VoxelIndex& extent, double size)
{
  if(!std::isfinite(size) || size <= 0)
    throw std::invalid_argument("Grid size must be finite and positive");
  OctreeBox box;
  for(std::size_t axis = 0; axis < 3; ++axis)
  {
    if(extent[axis] <= 0 || lower[axis] > std::numeric_limits<std::int64_t>::max() - extent[axis])
      throw std::out_of_range("Invalid or overflowing grid extent");
    box.lower[axis] = static_cast<double>(lower[axis]) * size;
    box.upper[axis] = static_cast<double>(lower[axis] + extent[axis]) * size;
    if(!std::isfinite(box.lower[axis]) || !std::isfinite(box.upper[axis]) ||
       !(box.lower[axis] < box.upper[axis]))
      throw std::out_of_range("Octree grid planes are not representable");
  }
  if(!std::isfinite(box.volume()) || box.volume() <= 0)
    throw std::out_of_range("Octree box volume is not representable");
  return box;
}

namespace
{
OctreeBox gridBox(const VoxelIndex& lower, std::int64_t width, double size)
{
  return octreeGridBox(lower, {width, width, width}, size);
}

int appendNode(OctreeOccupancy& octree, const VoxelIndex& lower,
               std::int64_t width, unsigned depth)
{
  if(octree.nodes.size() >= static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw std::length_error("Too many octree nodes");
  OctreeNode node;
  node.lower = lower;
  node.width = width;
  node.depth = depth;
  node.box = gridBox(lower, width, octree.voxel_size);
  octree.nodes.push_back(std::move(node));
  return static_cast<int>(octree.nodes.size() - 1);
}

void require(bool condition, const char* message)
{
  if(!condition)
    throw std::runtime_error(message);
}

bool fullByCount(const OctreeNode& node)
{
  // Never form width^3 unchecked: sparse roots may span up to 2^51 voxels
  // along an axis. Stop as soon as capacity would exceed occupied_count.
  if(node.width < 1 || node.occupied_count == 0) return false;
  const auto width = static_cast<std::uint64_t>(node.width);
  std::size_t capacity = 1;
  for(unsigned axis = 0; axis < 3; ++axis)
  {
    if(width > node.occupied_count / capacity) return false;
    capacity *= static_cast<std::size_t>(width);
  }
  return capacity == node.occupied_count;
}
}  // namespace

bool OctreeBox::contains(const HS::Vec3& point) const
{
  return point.x >= lower[0] && point.x <= upper[0] &&
         point.y >= lower[1] && point.y <= upper[1] &&
         point.z >= lower[2] && point.z <= upper[2];
}

double OctreeBox::volume() const
{
  return (upper[0] - lower[0]) * (upper[1] - lower[1]) * (upper[2] - lower[2]);
}

bool OctreeOccupancy::contains(const HS::Vec3& point) const
{
  // Closed child boxes can share a boundary; visit all touching occupied
  // children, not just the one selected by a half-open floor index. No division
  // is needed, including at the uppermost supported integer grid plane.
  const auto visit = [&](const auto& self, int id) -> bool {
    const auto& node = nodes[id];
    if(!node.box.contains(point)) return false;
    if(node.width == 1) return true;
    for(const int child : node.children)
      if(child >= 0 && self(self, child)) return true;
    return false;
  };
  return !nodes.empty() && visit(visit, 0);
}

OctreeOccupancy buildOctreeOccupancy(const std::vector<VoxelPoint>& raw_points,
                                     double voxel_size, unsigned max_depth)
{
  if(raw_points.empty())
    throw std::invalid_argument("Cannot build octree from an empty point cloud");
  if(max_depth > 52)
    throw std::invalid_argument("max_depth must be in [0, 52]");
  OctreeOccupancy octree;
  octree.voxel_size = voxel_size;
  std::vector<std::pair<VoxelIndex, std::size_t>> indexed;
  indexed.reserve(raw_points.size());
  for(std::size_t i = 0; i < raw_points.size(); ++i)
    indexed.emplace_back(voxelIndex(raw_points[i], voxel_size), i);
  std::sort(indexed.begin(), indexed.end());
  VoxelIndex lower = indexed.front().first, upper = lower;
  for(const auto& entry : indexed)
    for(std::size_t axis = 0; axis < 3; ++axis)
    {
      lower[axis] = std::min(lower[axis], entry.first[axis]);
      upper[axis] = std::max(upper[axis], entry.first[axis]);
    }
  std::int64_t extent = 1;
  for(std::size_t axis = 0; axis < 3; ++axis)
    extent = std::max(extent, upper[axis] - lower[axis] + 1);
  std::int64_t width = 1;
  while(width < extent)
  {
    width *= 2;  // voxelIndex limits absolute indices to < 2^50.
    ++octree.required_depth;
  }
  if(octree.required_depth > max_depth)
    throw std::invalid_argument("max_depth is too small for voxel_size: required " +
                                std::to_string(octree.required_depth) +
                                "; resolution will not be silently coarsened");
  appendNode(octree, lower, width, 0);
  for(std::size_t begin = 0; begin < indexed.size();)
  {
    std::size_t end = begin + 1;
    const VoxelIndex& key = indexed[begin].first;
    while(end < indexed.size() && indexed[end].first == key)
      ++end;
    octree.occupied_voxels.push_back(key);
    int id = 0;
    for(;;)
    {
      ++octree.nodes[id].occupied_count;
      octree.nodes[id].point_count += end - begin;
      if(octree.nodes[id].width == 1)
        break;
      const auto half = octree.nodes[id].width / 2;
      VoxelIndex child_lower = octree.nodes[id].lower;
      unsigned child = 0;
      for(unsigned axis = 0; axis < 3; ++axis)
        if(key[axis] >= child_lower[axis] + half)
        {
          child |= 1u << axis;
          child_lower[axis] += half;
        }
      if(octree.nodes[id].children[child] < 0)
      {
        // appendNode can reallocate nodes: retain indices, never references.
        const unsigned depth = octree.nodes[id].depth + 1;
        const int next = appendNode(octree, child_lower, half, depth);
        octree.nodes[id].children[child] = next;
      }
      id = octree.nodes[id].children[child];
    }
    auto& point_ids = octree.nodes[id].point_indices;
    point_ids.reserve(end - begin);
    for(std::size_t i = begin; i < end; ++i)
      point_ids.push_back(indexed[i].second);
    octree.leaf_nodes.push_back(id);
    begin = end;
  }
  return octree;
}

OctreePruning pruneOctree(const OctreeOccupancy& octree, bool exact)
{
  require(!octree.nodes.empty(), "Cannot prune an empty octree");
  OctreePruning result;
  result.exact = exact;
  const auto& nodes = octree.nodes;
  std::vector<unsigned char> full(nodes.size(), 0);
  // Children always have larger IDs: reverse creation order is a bottom-up
  // traversal. A full parent replaces eight full children: 8 cells / 48 plane
  // references become 1 cell / 6 references, preserving the exact box union.
  for(std::size_t i = nodes.size(); i-- > 0;)
  {
    const auto& node = nodes[i];
    full[i] = node.width == 1;
    if(!exact || node.width == 1) continue;
    bool all_full = true;
    for(const int child : node.children)
      all_full &= child >= 0 && full.at(static_cast<std::size_t>(child));
    if(all_full)
    {
      full[i] = 1;
      ++result.merge_operations;
    }
  }
  std::vector<int> owner(nodes.size(), -1), cell_id(nodes.size(), -1);
  result.active_nodes.resize(nodes.size(), 0);
  for(std::size_t i = 0; i < nodes.size(); ++i)
  {
    if(owner[i] < 0)
    {
      result.active_nodes[i] = 1;
      if(full[i]) owner[i] = static_cast<int>(i);
    }
    for(const int child : nodes[i].children)
      if(child >= 0) owner.at(child) = owner[i];
  }
  // Stable output ordering without sorting the cut: fine voxels are already
  // lexicographically ordered. Each source voxel is stored in exactly one cell.
  result.voxel_cell_ids.reserve(octree.leaf_nodes.size());
  for(const int leaf : octree.leaf_nodes)
  {
    const int selected = owner.at(leaf);
    require(selected >= 0, "Pruning left an uncovered occupied voxel");
    if(cell_id[selected] < 0)
    {
      cell_id[selected] = static_cast<int>(result.cell_nodes.size());
      result.cell_nodes.push_back(selected);
    }
    result.voxel_cell_ids.push_back(static_cast<std::size_t>(cell_id[selected]));
  }
  return result;
}

void validateOctreePruning(const OctreeOccupancy& octree, const OctreePruning& pruning)
{
  const auto& nodes = octree.nodes;
  const auto count = octree.leaf_nodes.size();
  require(!nodes.empty() && !pruning.cell_nodes.empty(), "Empty pruning cut");
  require(pruning.active_nodes.size() == nodes.size() && pruning.voxel_cell_ids.size() == count,
          "Pruning metadata sizes disagree with occupancy");
  std::vector<int> selected(nodes.size(), -1), owner(nodes.size(), -1);
  for(std::size_t i = 0; i < pruning.cell_nodes.size(); ++i)
  {
    const int id = pruning.cell_nodes[i];
    require(id >= 0 && static_cast<std::size_t>(id) < nodes.size(), "Invalid selected node ID");
    require(selected[id] < 0 && fullByCount(nodes[id]), "Duplicate or non-full selected node");
    require(pruning.exact || nodes[id].width == 1, "No-pruning mode selected an internal node");
    selected[id] = static_cast<int>(i);
  }
  for(std::size_t i = 0; i < nodes.size(); ++i)
  {
    const bool active = owner[i] < 0;
    require(pruning.active_nodes[i] == static_cast<unsigned char>(active), "Invalid pruned active mask");
    if(selected[i] >= 0)
    {
      require(active, "Selected parent and descendant overlap");
      owner[i] = selected[i];
    }
    else if(active)
    {
      require(nodes[i].width > 1, "Pruning missed an occupied leaf");
      require(!pruning.exact || !fullByCount(nodes[i]), "Exact cut left a recoverable full parent");
    }
    for(const int child : nodes[i].children)
      if(child >= 0) owner.at(child) = owner[i];
  }
  std::vector<std::size_t> covered(pruning.cell_nodes.size(), 0);
  std::size_t next_first_cell = 0;
  for(std::size_t i = 0; i < count; ++i)
  {
    const int id = owner.at(octree.leaf_nodes[i]);
    require(id >= 0 && pruning.voxel_cell_ids[i] == static_cast<std::size_t>(id),
            "Fine voxel provenance disagrees with pruning cut");
    if(covered[id]++ == 0)
      require(static_cast<std::size_t>(id) == next_first_cell++, "Pruned cells are not deterministically ordered");
  }
  for(std::size_t i = 0; i < covered.size(); ++i)
    require(covered[i] == nodes[pruning.cell_nodes[i]].occupied_count, "Selected cell coverage count mismatch");
  require(pruning.cell_nodes.size() <= count &&
          (count - pruning.cell_nodes.size()) % 7 == 0 &&
          pruning.merge_operations == (count - pruning.cell_nodes.size()) / 7,
          "Pruning merge/cost accounting mismatch");
}

std::vector<ConvexCluster> octreeAabbClusters(const OctreeOccupancy& octree,
                                             std::vector<HS::Vec3>& vertices)
{
  return octreeAabbClusters(octree, pruneOctree(octree, false), vertices);
}

std::vector<ConvexCluster> octreeAabbClusters(const OctreeOccupancy& octree,
                                             const OctreePruning& pruning,
                                             std::vector<HS::Vec3>& vertices)
{
  std::vector<OctreeBox> boxes;
  boxes.reserve(pruning.cell_nodes.size());
  for(const int id : pruning.cell_nodes) boxes.push_back(octree.nodes.at(id).box);
  return aabbClusters(boxes, pruning.voxel_cell_ids, vertices);
}

std::vector<ConvexCluster> aabbClusters(const std::vector<OctreeBox>& boxes,
                                       const std::vector<std::size_t>& voxel_cell_ids,
                                       std::vector<HS::Vec3>& vertices)
{
  const auto count = boxes.size();
  if(count > static_cast<std::size_t>(std::numeric_limits<int>::max()) / 7)
    throw std::length_error("Too many AABB cells for integer tree/plane IDs");
  std::vector<ConvexCluster> clusters;
  clusters.reserve(count);
  vertices.clear();
  vertices.reserve(8 * count);
  for(std::size_t i = 0; i < count; ++i)
  {
    const auto& box = boxes[i];
    ConvexCluster cluster;
    cluster.id = static_cast<int>(i);
    for(unsigned mask = 0; mask < 8; ++mask)
    {
      cluster.vertices.insert(vertices.size());
      vertices.push_back({(mask & 1) ? box.upper[0] : box.lower[0],
                          (mask & 2) ? box.upper[1] : box.lower[1],
                          (mask & 4) ? box.upper[2] : box.lower[2]});
    }
    for(unsigned axis = 0; axis < 3; ++axis)
      for(unsigned side = 0; side < 2; ++side)
      {
        std::array<double, 3> normal{}, point{};
        normal[axis] = side ? 1 : -1;
        point[axis] = side ? box.upper[axis] : box.lower[axis];
        Plane plane;
        plane.normal = {normal[0], normal[1], normal[2]};
        plane.point = {point[0], point[1], point[2]};
        plane.offset = normal[axis] * point[axis];
        plane.area = (box.upper[(axis + 1) % 3] - box.lower[(axis + 1) % 3]) *
                     (box.upper[(axis + 2) % 3] - box.lower[(axis + 2) % 3]);
        plane.source_face_ids.push_back(static_cast<int>(6 * i + 2 * axis + side));
        cluster.exact_planes.push_back(std::move(plane));
      }
    // Six independent support directions: no local plane is redundant.
    cluster.reduced_planes = cluster.exact_planes;
    clusters.push_back(std::move(cluster));
  }
  for(std::size_t i = 0; i < voxel_cell_ids.size(); ++i)
    clusters.at(voxel_cell_ids[i]).cells.push_back(i);
  return clusters;
}

void validateOctreeOccupancy(const OctreeOccupancy& octree,
                            const std::vector<VoxelPoint>& raw_points)
{
  const auto& nodes = octree.nodes;
  require(!nodes.empty(), "Octree has no root");
  require(octree.leaf_nodes.size() == octree.occupied_voxels.size(), "Leaf/key count mismatch");
  require(std::is_sorted(octree.occupied_voxels.begin(), octree.occupied_voxels.end()) &&
          std::adjacent_find(octree.occupied_voxels.begin(), octree.occupied_voxels.end()) ==
              octree.occupied_voxels.end(), "Occupied keys are not sorted/unique");
  std::vector<unsigned> parents(nodes.size(), 0), seen_points(raw_points.size(), 0);
  std::size_t leaves = 0;
  for(std::size_t id = 0; id < nodes.size(); ++id)
  {
    const auto& node = nodes[id];
    const auto expected = gridBox(node.lower, node.width, octree.voxel_size);
    require(node.box.lower == expected.lower && node.box.upper == expected.upper,
            "Node bounds disagree with integer grid");
    std::size_t occupied = 0, points = 0, children = 0;
    for(unsigned slot = 0; slot < 8; ++slot)
    {
      const int child_id = node.children[slot];
      if(child_id < 0)
        continue;
      require(static_cast<std::size_t>(child_id) > id &&
              static_cast<std::size_t>(child_id) < nodes.size(), "Invalid child ID");
      const auto& child = nodes[child_id];
      ++parents[child_id];
      ++children;
      VoxelIndex lower = node.lower;
      for(unsigned axis = 0; axis < 3; ++axis)
        if(slot & (1u << axis)) lower[axis] += node.width / 2;
      require(node.width > 1 && child.width == node.width / 2 &&
              child.lower == lower && child.depth == node.depth + 1, "Invalid child partition");
      occupied += child.occupied_count;
      points += child.point_count;
    }
    if(node.width == 1)
    {
      ++leaves;
      require(children == 0 && node.occupied_count == 1 &&
              node.point_count == node.point_indices.size() && node.point_count > 0,
              "Invalid occupied leaf");
      for(const auto point_id : node.point_indices)
      {
        require(point_id < raw_points.size(), "Invalid raw point ID");
        ++seen_points[point_id];
        require(voxelIndex(raw_points[point_id], octree.voxel_size) == node.lower,
                "Raw point assigned to wrong leaf");
      }
    }
    else
      require(children > 0 && node.point_indices.empty() &&
              occupied == node.occupied_count && points == node.point_count,
              "Internal octree counts disagree with children");
  }
  require(leaves == octree.leaf_nodes.size(), "Missing occupied leaves");
  require(nodes[0].point_count == raw_points.size() &&
          nodes[0].occupied_count == octree.occupied_voxels.size(), "Root counts mismatch");
  for(std::size_t i = 1; i < nodes.size(); ++i)
    require(parents[i] == 1, "Node is detached or has multiple parents");
  for(const auto count : seen_points) require(count == 1, "Raw point lost or duplicated");
  for(std::size_t i = 0; i < leaves; ++i)
  {
    const auto& node = nodes.at(octree.leaf_nodes[i]);
    require(node.width == 1 && node.lower == octree.occupied_voxels[i], "Leaf/key mismatch");
  }
}
}  // namespace rokae_demo
