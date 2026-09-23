#include <rokae_demo/sparse_voxel_occupancy.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace rokae_demo
{
SparseVoxelOccupancy buildSparseVoxelOccupancy(const std::vector<VoxelPoint>& points,double h)
{
  if(points.empty()) throw std::invalid_argument("Cannot voxelize an empty cloud");
  SparseVoxelOccupancy o; o.voxel_size=h; o.occupied_voxels.reserve(points.size());
  o.lower=o.upper=voxelIndex(points.front(),h);
  o.occupied_voxels.push_back(o.lower);
  for(std::size_t id=1;id<points.size();++id)
  {
    const auto key=voxelIndex(points[id],h); o.occupied_voxels.push_back(key);
    for(unsigned axis=0;axis<3;++axis)
    {
      o.lower[axis]=std::min(o.lower[axis],key[axis]);
      o.upper[axis]=std::max(o.upper[axis],key[axis]);
    }
  }
  auto& keys=o.occupied_voxels;
  const auto bits=[](std::uint64_t value) {
    unsigned result=0; while(value) { ++result; value>>=1; } return result;
  };
  unsigned widths[3]{};
  for(unsigned axis=0;axis<3;++axis)
    widths[axis]=bits(static_cast<std::uint64_t>(o.upper[axis])-
                     static_cast<std::uint64_t>(o.lower[axis]));
  std::size_t dense_cells=1;
  bool use_dense=true;
  for(unsigned axis=0;axis<3;++axis)
  {
    const auto width=static_cast<std::uint64_t>(o.upper[axis])-
                     static_cast<std::uint64_t>(o.lower[axis])+1;
    if(width>std::numeric_limits<std::size_t>::max() ||
       static_cast<std::size_t>(width)>std::numeric_limits<std::size_t>::max()/dense_cells)
    { use_dense=false; break; }
    o.dense_extent[axis]=static_cast<std::size_t>(width);
    dense_cells*=o.dense_extent[axis];
  }
  constexpr std::size_t absolute_dense_limit=16*1024*1024;
  const auto raw_relative_limit=points.size()<=std::numeric_limits<std::size_t>::max()/64
    ? std::max(std::size_t{4096},points.size()*64) : std::numeric_limits<std::size_t>::max();
  use_dense=use_dense && dense_cells<=absolute_dense_limit && dense_cells<=raw_relative_limit;
  if(use_dense)
  {
    o.dense_bits.assign((dense_cells+63)/64,0);
    for(unsigned axis=0;axis<3;++axis)
      o.boundary_plane_signs[axis].assign(o.dense_extent[axis]+1,0);
    for(const auto& key:keys)
    {
      const auto x=static_cast<std::uint64_t>(key[0])-static_cast<std::uint64_t>(o.lower[0]);
      const auto y=static_cast<std::uint64_t>(key[1])-static_cast<std::uint64_t>(o.lower[1]);
      const auto z=static_cast<std::uint64_t>(key[2])-static_cast<std::uint64_t>(o.lower[2]);
      const auto id=(static_cast<std::size_t>(x)*o.dense_extent[1]+static_cast<std::size_t>(y))*
        o.dense_extent[2]+static_cast<std::size_t>(z);
      o.dense_bits[id/64]|=std::uint64_t{1}<<(id%64);
    }
    keys.clear(); keys.reserve(std::min(points.size(),dense_cells));
    for(std::size_t word_id=0;word_id<o.dense_bits.size();++word_id)
    {
      auto word=o.dense_bits[word_id];
      while(word)
      {
#if defined(__GNUC__) || defined(__clang__)
        const auto bit=static_cast<unsigned>(__builtin_ctzll(word));
#else
        unsigned bit=0; while(!(word&(std::uint64_t{1}<<bit))) ++bit;
#endif
        const auto id=word_id*64+bit;
        if(id<dense_cells)
        {
          const auto x=id/(o.dense_extent[1]*o.dense_extent[2]);
          const auto rest=id%(o.dense_extent[1]*o.dense_extent[2]);
          const auto y=rest/o.dense_extent[2],z=rest%o.dense_extent[2];
          const std::size_t coordinate[3]{x,y,z};
          const std::size_t stride[3]{
            o.dense_extent[1]*o.dense_extent[2],o.dense_extent[2],1};
          for(unsigned axis=0;axis<3;++axis) for(const int side:{-1,1})
          {
            const bool in_bounds=side<0 ? coordinate[axis]>0 :
              coordinate[axis]+1<o.dense_extent[axis];
            const auto neighbor=side<0 ? id-stride[axis] : id+stride[axis];
            const bool neighbor_occupied=in_bounds &&
              ((o.dense_bits[neighbor/64]>>(neighbor%64))&1);
            if(!neighbor_occupied)
              o.boundary_plane_signs[axis][coordinate[axis]+(side>0)]|=
                static_cast<unsigned char>(side<0 ? 1 : 2);
          }
          keys.push_back({o.lower[0]+static_cast<std::int64_t>(x),
                          o.lower[1]+static_cast<std::int64_t>(y),
                          o.lower[2]+static_cast<std::int64_t>(z)});
        }
        word&=word-1;
      }
    }
    const auto unique_relative_limit=keys.size()<=std::numeric_limits<std::size_t>::max()/64
      ? std::max(std::size_t{4096},keys.size()*64) : std::numeric_limits<std::size_t>::max();
    if(dense_cells>unique_relative_limit)
    { o.dense_extent={}; o.dense_bits.clear(); o.boundary_plane_signs={}; }
  }
  else if(widths[0]+widths[1]+widths[2]<=63)
  {
    std::vector<std::uint64_t> packed; packed.reserve(keys.size());
    for(const auto& key:keys)
    {
      const auto x=static_cast<std::uint64_t>(key[0])-static_cast<std::uint64_t>(o.lower[0]);
      const auto y=static_cast<std::uint64_t>(key[1])-static_cast<std::uint64_t>(o.lower[1]);
      const auto z=static_cast<std::uint64_t>(key[2])-static_cast<std::uint64_t>(o.lower[2]);
      packed.push_back(((x<<widths[1])|y)<<widths[2]|z);
    }
    std::sort(packed.begin(),packed.end());
    packed.erase(std::unique(packed.begin(),packed.end()),packed.end());
    const auto z_mask=widths[2] ? (std::uint64_t{1}<<widths[2])-1 : 0;
    const auto y_mask=widths[1] ? (std::uint64_t{1}<<widths[1])-1 : 0;
    keys.clear(); keys.reserve(packed.size());
    for(const auto value:packed)
    {
      const auto z=value&z_mask;
      const auto y=(value>>widths[2])&y_mask;
      const auto x=value>>(widths[1]+widths[2]);
      keys.push_back({o.lower[0]+static_cast<std::int64_t>(x),
                      o.lower[1]+static_cast<std::int64_t>(y),
                      o.lower[2]+static_cast<std::int64_t>(z)});
    }
  }
  else
  {
    std::sort(keys.begin(),keys.end());
    keys.erase(std::unique(keys.begin(),keys.end()),keys.end());
  }
  VoxelIndex extent{}; std::int64_t width=1;
  for(unsigned a=0;a<3;++a) extent[a]=o.upper[a]-o.lower[a]+1;
  while(width<*std::max_element(extent.begin(),extent.end())) { width*=2; ++o.required_depth; }
  o.box=octreeGridBox(o.lower,extent,h);
  o.sampling_box=octreeGridBox(o.lower,{width,width,width},h);
  // Validate rounded fine boxes too; broad sparse extents do not prove this.
  for(const auto& k:keys) (void)octreeGridBox(k,{1,1,1},h);
  return o;
}

bool SparseVoxelOccupancy::contains(const HS::Vec3& p) const
{
  if(!std::isfinite(p.x)||!std::isfinite(p.y)||!std::isfinite(p.z))
    throw std::invalid_argument("Nonfinite voxel occupancy query");
  if(occupied_voxels.empty()||!box.contains(p)) return false;
  const double xyz[3]{p.x,p.y,p.z};
  std::array<std::array<std::int64_t,2>,3> candidates{}; unsigned count[3]{};
  for(unsigned a=0;a<3;++a)
  {
    // Find actual rounded grid planes, without division/floor range problems
    // at the largest supported upper boundary. At most two CLOSED cells/axis.
    auto lo=lower[a], hi=upper[a]+1;
    while(lo<hi)
    {
      const auto mid=lo+(hi-lo+1)/2;
      if(static_cast<double>(mid)*voxel_size<=xyz[a]) lo=mid; else hi=mid-1;
    }
    if(lo<=upper[a]) candidates[a][count[a]++]=lo;
    if(lo>lower[a] && static_cast<double>(lo)*voxel_size==xyz[a]) candidates[a][count[a]++]=lo-1;
  }
  const auto occupied=[&](const VoxelIndex& key) {
    if(dense_bits.empty()) return std::binary_search(occupied_voxels.begin(),occupied_voxels.end(),key);
    for(unsigned axis=0;axis<3;++axis)
      if(key[axis]<lower[axis] || key[axis]>upper[axis]) return false;
    const auto x=static_cast<std::uint64_t>(key[0])-static_cast<std::uint64_t>(lower[0]);
    const auto y=static_cast<std::uint64_t>(key[1])-static_cast<std::uint64_t>(lower[1]);
    const auto z=static_cast<std::uint64_t>(key[2])-static_cast<std::uint64_t>(lower[2]);
    const auto id=(static_cast<std::size_t>(x)*dense_extent[1]+static_cast<std::size_t>(y))*
      dense_extent[2]+static_cast<std::size_t>(z);
    return ((dense_bits[id/64]>>(id%64))&1)!=0;
  };
  for(unsigned x=0;x<count[0];++x) for(unsigned y=0;y<count[1];++y) for(unsigned z=0;z<count[2];++z)
    if(occupied({candidates[0][x],candidates[1][y],candidates[2][z]})) return true;
  return false;
}

void validateSparseVoxelOccupancy(const SparseVoxelOccupancy& o,const std::vector<VoxelPoint>& points)
{
  // Rebuild the raw-point key certificate independently of hierarchy metadata.
  const auto expected=buildSparseVoxelOccupancy(points,o.voxel_size);
  if(o.occupied_voxels!=expected.occupied_voxels || o.lower!=expected.lower || o.upper!=expected.upper ||
     o.dense_extent!=expected.dense_extent || o.dense_bits!=expected.dense_bits ||
     o.boundary_plane_signs!=expected.boundary_plane_signs ||
     o.box.lower!=expected.box.lower || o.box.upper!=expected.box.upper ||
     o.sampling_box.lower!=expected.sampling_box.lower || o.sampling_box.upper!=expected.sampling_box.upper ||
     o.required_depth!=expected.required_depth) throw std::runtime_error("Voxel occupancy certificate mismatch");
  for(const auto& p:points) if(!o.contains({p.x(),p.y(),p.z()}))
    throw std::runtime_error("Raw point outside closed voxel union");
}
} // namespace rokae_demo
