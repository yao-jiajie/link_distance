#include <rokae_demo/octree_rectangular_pruning.hpp>

#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace rokae_demo
{
namespace
{
void require(bool valid, const char* message)
{
  if(!valid) throw std::runtime_error(message);
}

unsigned population(unsigned mask)
{
  unsigned result = 0;
  for(; mask; mask >>= 1) result += mask & 1u;
  return result;
}

struct Pattern
{
  unsigned mask = 0;
  VoxelIndex lower{}, extent{};
};

struct PatternTable
{
  std::array<Pattern, 256> patterns{};
  std::array<unsigned, 256> cost{}, choice{};
  PatternTable()
  {
    std::vector<Pattern> candidates;
    for(int x0 = 0; x0 < 2; ++x0) for(int x1 = x0 + 1; x1 <= 2; ++x1)
      for(int y0 = 0; y0 < 2; ++y0) for(int y1 = y0 + 1; y1 <= 2; ++y1)
        for(int z0 = 0; z0 < 2; ++z0) for(int z1 = z0 + 1; z1 <= 2; ++z1)
        {
          Pattern p;
          p.lower = {x0, y0, z0}; p.extent = {x1-x0, y1-y0, z1-z0};
          for(int x = x0; x < x1; ++x) for(int y = y0; y < y1; ++y)
            for(int z = z0; z < z1; ++z) p.mask |= 1u << (x + 2*y + 4*z);
          candidates.push_back(p);
          patterns[p.mask] = p;
        }
    cost.fill(9); cost[0] = 0;
    for(unsigned mask = 1; mask < 256; ++mask)
    {
      const unsigned first = mask & (~mask + 1u);
      for(const auto& p : candidates)
        if((p.mask & first) && (mask & p.mask) == p.mask &&
           1 + cost[mask ^ p.mask] < cost[mask])
        {
          cost[mask] = 1 + cost[mask ^ p.mask];
          choice[mask] = p.mask;
        }
    }
  }
};

const PatternTable& table()
{
  static const PatternTable result;
  return result;
}

unsigned occupiedMask(const OctreeNode& node)
{
  unsigned result = 0;
  for(unsigned i = 0; i < 8; ++i) if(node.children[i] >= 0) result |= 1u << i;
  return result;
}

OctreeRectangularCell cubeCell(const OctreeOccupancy& octree, int id)
{
  const auto& node = octree.nodes.at(id);
  OctreeRectangularCell cell;
  cell.source_node = id;
  cell.lower = node.lower;
  cell.extent = {node.width, node.width, node.width};
  cell.box = node.box;
  cell.occupied_count = node.occupied_count;
  return cell;
}
}  // namespace

std::vector<unsigned> rectangularSiblingPartition(unsigned occupied_mask)
{
  if(occupied_mask > 255) throw std::invalid_argument("Sibling occupancy mask must fit eight bits");
  std::vector<unsigned> result;
  for(unsigned mask = occupied_mask; mask;)
  {
    const unsigned choice = table().choice[mask];
    result.push_back(choice);
    mask ^= choice;
  }
  return result;
}

OctreeCellPartition partitionOctreeCells(const OctreeOccupancy& octree,
                                         OctreePruning base, bool rectangles)
{
  require(!rectangles || base.exact, "Rectangular grouping requires the exact cubic base cut");
  OctreeCellPartition result;
  result.base = std::move(base);
  result.rectangles = rectangles;
  const auto& cut = result.base;
  std::vector<OctreeRectangularCell> provisional;
  std::vector<std::size_t> replacement(cut.cell_nodes.size());
  std::iota(replacement.begin(), replacement.end(), 0);
  std::vector<int> base_at_node(octree.nodes.size(), -1);
  for(std::size_t i = 0; i < cut.cell_nodes.size(); ++i)
  {
    provisional.push_back(cubeCell(octree, cut.cell_nodes[i]));
    base_at_node.at(cut.cell_nodes[i]) = static_cast<int>(i);
  }
  if(rectangles)
    for(std::size_t parent_id = 0; parent_id < octree.nodes.size(); ++parent_id)
    {
      const auto& parent = octree.nodes[parent_id];
      // Only finest sibling groups, outside already recovered full cubes.
      // Larger full cubes remain intact; no cross-parent or pairwise merge.
      if(parent.width != 2 || !cut.active_nodes[parent_id] || base_at_node[parent_id] >= 0) continue;
      const unsigned occupied = occupiedMask(parent);
      const auto masks = rectangularSiblingPartition(occupied);
      if(masks.size() < population(occupied)) ++result.rectangular_groups;
      for(const unsigned mask : masks)
      {
        if(population(mask) == 1) continue;
        const auto& pattern = table().patterns[mask];
        OctreeRectangularCell cell;
        cell.parent_node = static_cast<int>(parent_id);
        cell.sibling_mask = mask;
        cell.lower = parent.lower;
        cell.extent = pattern.extent;
        for(unsigned axis = 0; axis < 3; ++axis) cell.lower[axis] += pattern.lower[axis];
        cell.box = octreeGridBox(cell.lower, cell.extent, octree.voxel_size);
        cell.occupied_count = population(mask);
        const auto id = provisional.size();
        provisional.push_back(std::move(cell));
        ++result.rectangular_merge_count;
        for(unsigned slot = 0; slot < 8; ++slot)
          if(mask & (1u << slot))
          {
            const int base_id = base_at_node.at(parent.children[slot]);
            require(base_id >= 0, "Rectangular sibling is not an active fine cell");
            replacement.at(base_id) = id;
          }
      }
    }
  // Assign final IDs by first source voxel, retaining stable geometry/tree IDs.
  std::vector<int> final_id(provisional.size(), -1);
  result.base_cell_ids.resize(cut.cell_nodes.size());
  result.voxel_cell_ids.reserve(cut.voxel_cell_ids.size());
  for(const auto base_id : cut.voxel_cell_ids)
  {
    const auto candidate = replacement.at(base_id);
    if(final_id[candidate] < 0)
    {
      final_id[candidate] = static_cast<int>(result.cells.size());
      result.cells.push_back(provisional[candidate]);
    }
    const auto id = static_cast<std::size_t>(final_id[candidate]);
    result.base_cell_ids.at(base_id) = id;
    result.voxel_cell_ids.push_back(id);
  }
  return result;
}

void validateOctreeCellPartition(const OctreeOccupancy& octree,
                                const OctreeCellPartition& partition)
{
  validateOctreePruning(octree, partition.base);
  const auto& cut = partition.base;
  require(!partition.cells.empty() && partition.cells.size() <= cut.cell_nodes.size(), "Invalid rectangular cell count");
  require(!partition.rectangles || cut.exact, "Rectangle mode has non-exact base cut");
  require(partition.voxel_cell_ids.size() == octree.occupied_voxels.size() &&
          partition.base_cell_ids.size() == cut.cell_nodes.size(), "Rectangular provenance size mismatch");
  std::vector<int> base_at_node(octree.nodes.size(), -1);
  for(std::size_t i = 0; i < cut.cell_nodes.size(); ++i) base_at_node[cut.cell_nodes[i]] = static_cast<int>(i);
  std::size_t rectangles = 0;
  for(const auto& cell : partition.cells)
  {
    const auto expected = octreeGridBox(cell.lower, cell.extent, octree.voxel_size);
    require(cell.box.lower == expected.lower && cell.box.upper == expected.upper, "Rectangle grid bounds mismatch");
    std::size_t capacity = 1;
    for(const auto width : cell.extent)
    {
      require(static_cast<std::uint64_t>(width) <= cell.occupied_count / capacity, "Rectangle capacity exceeds occupancy");
      capacity *= static_cast<std::size_t>(width);
    }
    require(capacity == cell.occupied_count, "Rectangle is not full");
    if(cell.source_node >= 0)
    {
      const auto& node = octree.nodes.at(cell.source_node);
      require(base_at_node.at(cell.source_node) >= 0 && cell.parent_node == -1 && cell.sibling_mask == 0 &&
              cell.lower == node.lower && cell.extent == VoxelIndex{node.width,node.width,node.width} &&
              cell.occupied_count == node.occupied_count, "Invalid retained cubic cell");
    }
    else
    {
      ++rectangles;
      require(partition.rectangles && cell.source_node == -1 && cell.parent_node >= 0, "Invalid virtual rectangle");
      const auto& parent = octree.nodes.at(cell.parent_node);
      require(parent.width == 2 && cut.active_nodes.at(cell.parent_node) && base_at_node[cell.parent_node] < 0,
              "Rectangle is not within an active finest sibling group");
      require(cell.sibling_mask > 0 && cell.sibling_mask < 255, "Invalid rectangle sibling mask");
      const auto& pattern = table().patterns[cell.sibling_mask];
      require(pattern.mask == cell.sibling_mask && population(pattern.mask) > 1 &&
              (pattern.mask & occupiedMask(parent)) == pattern.mask && cell.extent == pattern.extent,
              "Nonrectangular or unoccupied sibling selection");
      for(unsigned axis = 0; axis < 3; ++axis)
        require(cell.lower[axis] == parent.lower[axis] + pattern.lower[axis], "Rectangle parent bounds mismatch");
    }
  }
  // Each distinct source voxel is assigned once, lies within its assigned box,
  // and each box receives its full integer capacity. Together these prove exact
  // coverage and disjoint interiors without O(cells^2) intersection checks.
  std::vector<std::size_t> covered(partition.cells.size(), 0);
  std::size_t next_id = 0;
  for(std::size_t i = 0; i < octree.occupied_voxels.size(); ++i)
  {
    const auto id = partition.voxel_cell_ids[i];
    require(id < partition.cells.size() && partition.base_cell_ids.at(cut.voxel_cell_ids[i]) == id,
            "Rectangular source ownership mismatch");
    const auto& cell = partition.cells[id];
    const auto& key = octree.occupied_voxels[i];
    for(unsigned axis = 0; axis < 3; ++axis)
      require(key[axis] >= cell.lower[axis] && key[axis] < cell.lower[axis] + cell.extent[axis],
              "Source voxel lies outside its rectangle");
    if(covered[id]++ == 0) require(id == next_id++, "Rectangle ordering is not deterministic");
  }
  for(std::size_t i = 0; i < covered.size(); ++i)
    require(covered[i] == partition.cells[i].occupied_count, "Incomplete rectangular coverage");
  require(rectangles == partition.rectangular_merge_count, "Rectangle merge count mismatch");
  std::size_t groups = 0;
  if(partition.rectangles)
    for(std::size_t i = 0; i < octree.nodes.size(); ++i)
    {
      const auto& parent = octree.nodes[i];
      if(parent.width != 2 || !cut.active_nodes[i] || base_at_node[i] >= 0) continue;
      std::vector<std::size_t> ids;
      for(const int child : parent.children)
        if(child >= 0) ids.push_back(partition.base_cell_ids.at(base_at_node.at(child)));
      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
      const auto mask = occupiedMask(parent);
      require(ids.size() == table().cost[mask], "Sibling partition is not locally minimal");
      groups += ids.size() < population(mask);
    }
  else
    require(partition.cells.size() == cut.cell_nodes.size(), "Disabled rectangle mode changed the cut");
  require(groups == partition.rectangular_groups, "Rectangular group count mismatch");
}

std::vector<ConvexCluster> octreeAabbClusters(const OctreeCellPartition& partition,
                                             std::vector<HS::Vec3>& vertices)
{
  std::vector<OctreeBox> boxes;
  boxes.reserve(partition.cells.size());
  for(const auto& cell : partition.cells) boxes.push_back(cell.box);
  return aabbClusters(boxes, partition.voxel_cell_ids, vertices);
}
}  // namespace rokae_demo
