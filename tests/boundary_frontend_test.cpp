#include "boundary_frontend_checks.hpp"
#include <iostream>
#include <random>

using namespace rokae_demo;
using namespace boundary_frontend_checks;
namespace ref = boundary_frontend_reference;

template<class F> void rejects(F f)
{
  bool rejected = false;
  try { f(); } catch(const std::exception&) { rejected = true; }
  require(rejected, "Invalid frontend input accepted");
}

void exercise(std::vector<VoxelIndex> keys)
{
  std::sort(keys.begin(),keys.end()); keys.erase(std::unique(keys.begin(),keys.end()),keys.end());
  const auto old = ref::extractBoundary(keys,nullptr,keys.size());
  auto current = extractVoxelBoundary(keys);
  sameBoundary(old,current); validateVoxelBoundary(keys,current);
  sameRegistry(ref::makeRegistry(old),buildBoundaryPlaneRegistry(current));

  // Exercise supplied owner metadata independently of geometry.
  OctreeOccupancy occupancy; occupancy.occupied_voxels = keys;
  OctreeCellPartition partition; partition.cells.resize(3);
  for(std::size_t i = 0; i < keys.size(); ++i) partition.voxel_cell_ids.push_back(i%3);
  const auto owned = extractOctreeBoundary(occupancy,partition);
  sameBoundary(ref::extractBoundary(keys,&partition.voxel_cell_ids,partition.cells.size()),owned);
  validateOctreeBoundary(occupancy,partition,owned);

  // Registry callers need not pass extraction order. IDs within each support
  // must still follow the caller's patch order, including both normal signs.
  std::reverse(current.patches.begin(),current.patches.end());
  sameRegistry(ref::makeRegistry(current),buildBoundaryPlaneRegistry(current));
  std::mt19937 random(1147); std::shuffle(current.patches.begin(),current.patches.end(),random);
  sameRegistry(ref::makeRegistry(current),buildBoundaryPlaneRegistry(current));
}

int main()
{
  try
  {
    for(unsigned mask = 1; mask < 256; ++mask)
    {
      std::vector<VoxelIndex> keys;
      for(unsigned bit = 0; bit < 8; ++bit) if(mask & (1u<<bit))
        keys.push_back({int(bit&1)-1,int((bit>>1)&1)-1,int((bit>>2)&1)-1});
      exercise(keys);
    }
    // Runs split, join, change endpoints, disappear, and reappear after a gap.
    for(unsigned axis = 0; axis < 3; ++axis)
    {
      std::vector<VoxelIndex> keys;
      const std::vector<std::vector<int>> rows{{0,1,4,5,8},{0,1,4,5,8},{0,1,2,4,8},{1,2,4,5,8},{},
                                             {1,2,4,5,8},{0,1,2,3,4,5,6,7,8}};
      for(std::size_t v = 0; v < rows.size(); ++v) for(const auto u : rows[v])
      {
        VoxelIndex key{}; key[axis] = -2; key[(axis+1)%3] = u-4; key[(axis+2)%3] = static_cast<std::int64_t>(v)-3;
        keys.push_back(key);
      }
      exercise(keys);
    }
    std::vector<VoxelIndex> shell, slab;
    for(int x = -4; x <= 4; ++x) for(int y = -4; y <= 4; ++y) for(int z = -4; z <= 4; ++z)
    {
      if(std::abs(x)==4 || std::abs(y)==4 || std::abs(z)==4) shell.push_back({x,y,z});
      if(z==0) slab.push_back({x,y,z});
    }
    exercise(shell); exercise(slab);
    const auto low = std::numeric_limits<std::int64_t>::min()+1;
    const auto high = std::numeric_limits<std::int64_t>::max()-1;
    // Use integer geometry directly: converting these coordinates to double
    // would collapse neighboring planes, outside the world-coordinate API.
    exercise({{low,low,low},{low+1,low,low},{high,high,high},{high-1,high,high}});
    exercise({{-1000000000,0,0},{1000000000,0,0}});
    std::mt19937 random(76411);
    for(unsigned trial = 0; trial < 100; ++trial)
    {
      std::vector<VoxelIndex> keys{{0,0,0}};
      for(int x = -3; x <= 3; ++x) for(int y = -3; y <= 3; ++y) for(int z = -3; z <= 3; ++z)
        if(random()%5 <= trial%4) keys.push_back({x,y,z});
      exercise(keys);
    }
    rejects([] { extractVoxelBoundary({}); });
    rejects([] { extractVoxelBoundary({{0,0,0},{0,0,0}}); });
    rejects([] { extractVoxelBoundary({{1,0,0},{0,0,0}}); });
    for(const auto edge : {std::numeric_limits<std::int64_t>::min(),std::numeric_limits<std::int64_t>::max()})
      rejects([&] { extractVoxelBoundary({{edge,0,0}}); });
    auto bad = extractVoxelBoundary({{0,0,0}}); bad.patches[0].axis = 3;
    rejects([&] { buildBoundaryPlaneRegistry(bad); });
    std::cout << "Frontend: exact ordered faces, patches, owners and supports match legacy maps\n";
  }
  catch(const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
