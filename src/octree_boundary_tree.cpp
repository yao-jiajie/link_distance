#include <rokae_demo/octree_boundary_tree.hpp>

#include <chrono>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <unordered_set>

namespace rokae_demo
{
namespace
{
using Clock = std::chrono::steady_clock;
double ms(Clock::time_point start) { return std::chrono::duration<double,std::milli>(Clock::now()-start).count(); }
void require(bool v, const char* message) { if(!v) throw std::runtime_error(message); }
std::size_t at(const std::vector<std::int64_t>& c, std::int64_t x)
{
  const auto i = std::lower_bound(c.begin(),c.end(),x);
  require(i != c.end() && *i == x,"Boundary endpoint has no axis-normal support plane");
  return static_cast<std::size_t>(i-c.begin());
}
BoundaryPlaneRegistry makeRegistry(const OctreeBoundary& b)
{
  BoundaryPlaneRegistry r;
  for(const auto& p : b.patches)
  {
    require(p.axis < 3,"Invalid boundary support axis");
    r.coordinates[p.axis].push_back(p.plane);
  }
  std::array<std::size_t,3> offsets{};
  std::size_t count = 0;
  for(unsigned axis = 0; axis < 3; ++axis)
  {
    auto& c = r.coordinates[axis];
    std::sort(c.begin(),c.end()); c.erase(std::unique(c.begin(),c.end()),c.end());
    require(c.size() >= 2,"Missing boundary extent");
    offsets[axis] = count; count += c.size();
  }
  r.supports.reserve(count);
  for(unsigned axis = 0; axis < 3; ++axis)
    for(const auto coordinate : r.coordinates[axis]) r.supports.push_back({axis,coordinate,0,{}});
  struct DenseLookup
  {
    std::int64_t lower=0;
    std::vector<std::size_t> positions;
  };
  std::array<DenseLookup,3> lookups;
  const auto patch_scaled=b.patches.size()<=std::numeric_limits<std::size_t>::max()/4
    ? b.patches.size()*4 : std::numeric_limits<std::size_t>::max();
  const auto dense_limit=std::max(std::size_t{4096},patch_scaled);
  for(unsigned axis=0;axis<3;++axis)
  {
    const auto& coordinates=r.coordinates[axis];
    const auto span=static_cast<std::uint64_t>(coordinates.back())-
                    static_cast<std::uint64_t>(coordinates.front());
    if(span>dense_limit || span==std::numeric_limits<std::size_t>::max()) continue;
    auto& lookup=lookups[axis]; lookup.lower=coordinates.front();
    lookup.positions.assign(static_cast<std::size_t>(span)+1,
                            std::numeric_limits<std::size_t>::max());
    for(std::size_t id=0;id<coordinates.size();++id)
    {
      const auto position=static_cast<std::uint64_t>(coordinates[id])-
                          static_cast<std::uint64_t>(lookup.lower);
      lookup.positions[static_cast<std::size_t>(position)]=id;
    }
  }
  const auto supportId=[&](const OctreeBoundaryPatch& patch) {
    const auto& coordinates=r.coordinates[patch.axis]; const auto& lookup=lookups[patch.axis];
    std::size_t position;
    if(!lookup.positions.empty())
    {
      const auto dense=static_cast<std::uint64_t>(patch.plane)-
                       static_cast<std::uint64_t>(lookup.lower);
      require(dense<lookup.positions.size(),"Boundary support outside dense registry lookup");
      position=lookup.positions[static_cast<std::size_t>(dense)];
      require(position!=std::numeric_limits<std::size_t>::max(),
              "Boundary support missing from dense registry lookup");
    }
    else
    {
      const auto found=std::lower_bound(coordinates.begin(),coordinates.end(),patch.plane);
      require(found!=coordinates.end() && *found==patch.plane,"Boundary support missing from registry");
      position=static_cast<std::size_t>(found-coordinates.begin());
    }
    return offsets[patch.axis]+position;
  };
  std::vector<std::size_t> support_ids; support_ids.reserve(b.patches.size());
  std::vector<std::size_t> provenance_counts(count,0);
  for(std::size_t id = 0; id < b.patches.size(); ++id)
  {
    const auto support=supportId(b.patches[id]); support_ids.push_back(support);
    ++provenance_counts[support];
  }
  for(std::size_t support=0;support<count;++support)
    r.supports[support].patches.reserve(provenance_counts[support]);
  // Appending in input order retains provenance for shuffled patches and
  // supports that occur with both normal signs.
  for(std::size_t id=0;id<b.patches.size();++id)
  {
    const auto& patch=b.patches[id]; auto& support=r.supports[support_ids[id]];
    support.outward_signs|=patch.sign<0 ? 1 : 2; support.patches.push_back(id);
  }
  return r;
}
struct SparseIndexHash
{
  std::size_t operator()(const VoxelIndex& key) const
  {
    std::size_t hash=0;
    for(const auto coordinate:key)
    {
      const auto value=std::hash<std::int64_t>{}(coordinate);
      hash^=value+std::size_t{0x9e3779b9}+(hash<<6)+(hash>>2);
    }
    return hash;
  }
};
BoundaryPlaneRegistry makeRegistry(const SparseVoxelOccupancy& occupancy)
{
  require(!occupancy.occupied_voxels.empty(),"Cannot register empty sparse occupancy");
  for(std::size_t i=0;i<occupancy.occupied_voxels.size();++i)
  {
    const auto& key=occupancy.occupied_voxels[i];
    require(!i || occupancy.occupied_voxels[i-1]<key,
            "Sparse registry occupancy must be sorted/unique");
    for(unsigned axis=0;axis<3;++axis)
      require(key[axis]>std::numeric_limits<std::int64_t>::min() &&
              key[axis]<std::numeric_limits<std::int64_t>::max(),
              "Sparse registry neighbor coordinate overflow");
  }

  const bool cached_signs=!occupancy.boundary_plane_signs[0].empty();
  for(unsigned axis=1;axis<3;++axis)
    require(occupancy.boundary_plane_signs[axis].empty()==!cached_signs,
            "Incomplete cached sparse registry signs");
  if(cached_signs)
  {
    BoundaryPlaneRegistry result;
    result.patch_provenance_complete=false;
    for(unsigned axis=0;axis<3;++axis)
    {
      const auto extent=static_cast<std::uint64_t>(occupancy.upper[axis])-
                        static_cast<std::uint64_t>(occupancy.lower[axis])+1;
      require(extent<std::numeric_limits<std::size_t>::max() &&
              occupancy.boundary_plane_signs[axis].size()==static_cast<std::size_t>(extent)+1,
              "Invalid cached sparse registry axis span");
      const auto& signs=occupancy.boundary_plane_signs[axis];
      for(std::size_t offset=0;offset<signs.size();++offset) if(signs[offset])
      {
        const auto coordinate=occupancy.lower[axis]+static_cast<std::int64_t>(offset);
        result.coordinates[axis].push_back(coordinate);
        result.supports.push_back({axis,coordinate,signs[offset],{}});
      }
      require(result.coordinates[axis].size()>=2,"Missing cached sparse registry boundary extent");
    }
    return result;
  }

  std::array<std::vector<unsigned char>,3> signs;
  std::array<std::map<std::int64_t,unsigned char>,3> sparse_signs;
  constexpr std::size_t absolute_axis_limit=16*1024*1024;
  const auto relative_axis_limit=occupancy.occupied_voxels.size()<=std::numeric_limits<std::size_t>::max()/64
    ? std::max(std::size_t{4096},occupancy.occupied_voxels.size()*64) : std::numeric_limits<std::size_t>::max();
  for(unsigned axis=0;axis<3;++axis)
  {
    const auto extent=static_cast<std::uint64_t>(occupancy.upper[axis])-
                      static_cast<std::uint64_t>(occupancy.lower[axis])+1;
    require(extent<std::numeric_limits<std::size_t>::max(),
            "Sparse registry axis extent overflow");
    if(extent+1<=absolute_axis_limit && extent+1<=relative_axis_limit)
      signs[axis].assign(static_cast<std::size_t>(extent)+1,0);
  }

  std::unordered_set<VoxelIndex,SparseIndexHash> sparse_lookup;
  if(occupancy.dense_bits.empty())
  {
    require(occupancy.occupied_voxels.size()<=std::numeric_limits<std::size_t>::max()/2,
            "Sparse registry lookup capacity overflow");
    sparse_lookup.reserve(occupancy.occupied_voxels.size()*2);
    sparse_lookup.insert(occupancy.occupied_voxels.begin(),occupancy.occupied_voxels.end());
  }
  const auto occupied=[&](const VoxelIndex& key) {
    for(unsigned axis=0;axis<3;++axis)
      if(key[axis]<occupancy.lower[axis] || key[axis]>occupancy.upper[axis]) return false;
    if(occupancy.dense_bits.empty()) return sparse_lookup.find(key)!=sparse_lookup.end();
    const auto x=static_cast<std::uint64_t>(key[0])-static_cast<std::uint64_t>(occupancy.lower[0]);
    const auto y=static_cast<std::uint64_t>(key[1])-static_cast<std::uint64_t>(occupancy.lower[1]);
    const auto z=static_cast<std::uint64_t>(key[2])-static_cast<std::uint64_t>(occupancy.lower[2]);
    const auto id=(static_cast<std::size_t>(x)*occupancy.dense_extent[1]+
                   static_cast<std::size_t>(y))*occupancy.dense_extent[2]+
                  static_cast<std::size_t>(z);
    return ((occupancy.dense_bits[id/64]>>(id%64))&1)!=0;
  };
  for(const auto& key:occupancy.occupied_voxels)
    for(unsigned axis=0;axis<3;++axis) for(const int side:{-1,1})
    {
      auto neighbor=key; neighbor[axis]+=side;
      if(occupied(neighbor)) continue;
      const auto plane=key[axis]+(side>0);
      const auto sign=static_cast<unsigned char>(side<0 ? 1 : 2);
      if(signs[axis].empty()) sparse_signs[axis][plane]|=sign;
      else
      {
        const auto offset=static_cast<std::uint64_t>(plane)-
                          static_cast<std::uint64_t>(occupancy.lower[axis]);
        require(offset<signs[axis].size(),"Sparse registry plane outside axis span");
        signs[axis][static_cast<std::size_t>(offset)]|=sign;
      }
    }

  BoundaryPlaneRegistry result;
  result.patch_provenance_complete=false;
  for(unsigned axis=0;axis<3;++axis)
  {
    if(signs[axis].empty())
      for(const auto& [coordinate,sign]:sparse_signs[axis])
      {
        result.coordinates[axis].push_back(coordinate);
        result.supports.push_back({axis,coordinate,sign,{}});
      }
    else
      for(std::size_t offset=0;offset<signs[axis].size();++offset)
        if(signs[axis][offset])
        {
          const auto coordinate=occupancy.lower[axis]+static_cast<std::int64_t>(offset);
          result.coordinates[axis].push_back(coordinate);
          result.supports.push_back({axis,coordinate,signs[axis][offset],{}});
        }
    require(result.coordinates[axis].size()>=2,"Missing sparse registry boundary extent");
  }
  return result;
}
void registry(const OctreeBoundary& b, BoundaryTreeResult& r)
{
  auto planes = makeRegistry(b);
  r.coordinates = std::move(planes.coordinates); r.supports = std::move(planes.supports);
}
// y/z denote local transverse coordinates, not necessarily world Y/Z.
struct Rectangle { std::size_t y0,y1,z0,z1; };
std::vector<Rectangle> rectangles(const std::vector<unsigned char>& occupied, std::size_t ny, std::size_t nz, bool transpose)
{
  std::vector<Rectangle> out;
  std::map<std::pair<std::size_t,std::size_t>,std::size_t> previous;
  const auto rows = transpose ? nz : ny, cols = transpose ? ny : nz;
  for(std::size_t row = 0; row < rows; ++row)
  {
    std::map<std::pair<std::size_t,std::size_t>,std::size_t> current;
    const auto blocked = [&](std::size_t col) { return occupied[transpose ? col*nz+row : row*nz+col]; };
    for(std::size_t col = 0; col < cols;)
    {
      if(blocked(col)) { ++col; continue; }
      const auto start = col;
      while(col < cols && !blocked(col)) ++col;
      const auto key = std::make_pair(start,col);
      const auto old = previous.find(key);
      std::size_t id;
      if(old != previous.end())
      { id = old->second; if(transpose) out[id].z1 = row+1; else out[id].y1 = row+1; }
      else
      {
        id = out.size();
        out.push_back(transpose ? Rectangle{start,col,row,row+1} : Rectangle{row,row+1,start,col});
      }
      current.emplace(key,id);
    }
    previous = std::move(current);
  }
  return out;
}
void mergeAxis(std::vector<FreePrism3D>& boxes, unsigned axis)
{
  const unsigned u = (axis+1)%3, v = (axis+2)%3;
  const auto key = [&](const FreePrism3D& b) { return std::make_tuple(b.lower[u],b.upper[u],b.lower[v],b.upper[v],b.lower[axis],b.upper[axis]); };
  std::sort(boxes.begin(),boxes.end(),[&](const auto& a,const auto& b) { return key(a)<key(b); });
  std::size_t n = 0;
  for(const auto b : boxes)
  {
    if(n && boxes[n-1].lower[u] == b.lower[u] && boxes[n-1].upper[u] == b.upper[u] &&
       boxes[n-1].lower[v] == b.lower[v] && boxes[n-1].upper[v] == b.upper[v] && boxes[n-1].upper[axis] == b.lower[axis])
      boxes[n-1].upper[axis] = b.upper[axis];
    else boxes[n++] = b;
  }
  boxes.resize(n);
}
TreeRepresentation closureTree(const BoundaryTreeResult& result, double h,
                               std::vector<std::pair<std::size_t,int>>& provenance)
{
  std::vector<OctreeBox> boxes;
  for(const auto& p : result.prisms)
    boxes.push_back(octreeGridBox(p.lower,{p.upper[0]-p.lower[0],p.upper[1]-p.lower[1],p.upper[2]-p.lower[2]},h));
  std::vector<std::size_t> ids(boxes.size()); std::iota(ids.begin(),ids.end(),0);
  std::vector<HS::Vec3> vertices;
  std::vector<ConvexCluster> clusters;
  if(!boxes.empty()) clusters = aabbClusters(boxes,ids,vertices);
  // Six unbounded exterior halfspaces: no fabricated infinite AABB vertices.
  for(unsigned axis = 0; axis < 3; ++axis) for(unsigned high = 0; high < 2; ++high)
  {
    const int sign = high ? -1 : 1;
    const double coordinate = static_cast<double>(high ? result.coordinates[axis].back() : result.coordinates[axis].front())*h;
    Plane plane; double n[3]{}, p[3]{}; n[axis] = sign; p[axis] = coordinate;
    plane.normal = {n[0],n[1],n[2]}; plane.point = {p[0],p[1],p[2]}; plane.offset = sign*coordinate;
    ConvexCluster c; c.id = static_cast<int>(clusters.size()); c.exact_planes = c.reduced_planes = {plane};
    clusters.push_back(c);
  }
  auto tree = buildClusterTree(clusters,0,"boundary_support",true);
  tree.simplified = HS::simplifyToFixedPoint(tree.raw);
  for(auto& leaf : tree.leaves)
  {
    const double n[3]{leaf.normal.x,leaf.normal.y,leaf.normal.z};
    unsigned axis = 0; while(axis < 3 && n[axis] == 0) ++axis;
    require(axis < 3 && (n[axis] == -1 || n[axis] == 1),"Non-axis boundary leaf");
    const auto s = std::find_if(result.supports.begin(),result.supports.end(),[&](const auto& p) {
      return p.axis == axis && n[axis]*static_cast<double>(p.coordinate)*h == leaf.offset;
    });
    require(s != result.supports.end(),"Artificial support coordinate in boundary tree");
    const auto id = static_cast<std::size_t>(s-result.supports.begin());
    provenance.emplace_back(id,static_cast<int>(n[axis]));
    leaf.source_face_ids.clear();
  }
  return tree;
}
void checkTree(const TreeRepresentation& a, const TreeRepresentation& b)
{
  require(a.leaves.size() == b.leaves.size() && a.leaf_references == b.leaf_references &&
          HS::structuralKey(a.raw) == HS::structuralKey(b.raw) && HS::structuralKey(a.simplified) == HS::structuralKey(b.simplified),
          "Boundary tree structure changed");
  for(std::size_t i = 0; i < a.leaves.size(); ++i)
  {
    const auto& x = a.leaves[i]; const auto& y = b.leaves[i];
    require(x.id == y.id && x.normal.x == y.normal.x && x.normal.y == y.normal.y && x.normal.z == y.normal.z && x.offset == y.offset,
            "Boundary tree coefficients changed");
  }
}
} // namespace

BoundaryPlaneRegistry buildBoundaryPlaneRegistry(const OctreeBoundary& b) { return makeRegistry(b); }
BoundaryPlaneRegistry buildBoundaryPlaneRegistry(const SparseVoxelOccupancy& o) { return makeRegistry(o); }

static BoundaryTreeResult buildBoundaryTreeImpl(double voxel_size, const OctreeBoundary& boundary,
                                                const BoundaryTreeOptions& options,
                                                const std::function<TreeRepresentation()>& fallback_tree)
{
  require(options.max_compressed_cells && options.max_event_updates && options.max_free_prisms && options.max_merge_passes,
          "Boundary work limits must be positive");
  require(options.sweep_axis < 3,"Sweep axis must be X/Y/Z");
  BoundaryTreeResult r;
  r.sweep_axis = options.sweep_axis;
  const unsigned axis = options.sweep_axis, u = (axis+1)%3, v = (axis+2)%3;
  auto start = Clock::now(); registry(boundary,r); r.registry_ms = ms(start);
  const auto fallback = [&](const char* reason) {
    r.fallback = true; r.fallback_reason = reason; r.prisms.clear(); r.leaf_sources.clear();
    const auto t = Clock::now(); r.tree = fallback_tree(); r.tree_ms = ms(t);
    return r;
  };
  const auto ns = r.coordinates[axis].size()-1, nu = r.coordinates[u].size()-1, nv = r.coordinates[v].size()-1;
  std::size_t product = 1;
  for(const auto n : {ns,nu,nv})
  {
    if(n > std::numeric_limits<std::size_t>::max()/product) return fallback("compressed_cell_overflow");
    product *= n;
  }
  r.compressed_cells = product;
  if(product > options.max_compressed_cells) return fallback("compressed_cell_limit");
  start = Clock::now();
  std::vector<std::vector<std::size_t>> events(ns+1);
  for(std::size_t i = 0; i < boundary.patches.size(); ++i)
    if(boundary.patches[i].axis == axis) events[at(r.coordinates[axis],boundary.patches[i].plane)].push_back(i);
  std::vector<unsigned char> occupied(nu*nv,0);
  for(std::size_t slab = 0; slab <= ns; ++slab)
  {
    for(const auto id : events[slab])
    {
      const auto& p = boundary.patches[id];
      const auto y0 = at(r.coordinates[u],p.u0), y1 = at(r.coordinates[u],p.u1);
      const auto z0 = at(r.coordinates[v],p.v0), z1 = at(r.coordinates[v],p.v1);
      const auto updates = (y1-y0)*(z1-z0);
      if(updates > options.max_event_updates-r.event_updates)
      { r.sweep_ms = ms(start)-r.partition_ms; return fallback("event_update_limit"); }
      r.event_updates += updates;
      for(auto y = y0; y < y1; ++y) for(auto z = z0; z < z1; ++z)
      {
        auto& state = occupied[y*nv+z];
        if((p.sign < 0 && state) || (p.sign > 0 && !state))
          throw std::runtime_error("Boundary event mismatch: patch="+std::to_string(id)+" axis="+std::to_string(axis)+
                                   " plane="+std::to_string(p.plane)+" u_index="+std::to_string(y)+
                                   " v_index="+std::to_string(z)+" sign="+std::to_string(p.sign));
        state = p.sign < 0;
      }
    }
    if(slab == ns) break;
    const auto t = Clock::now();
    auto a = rectangles(occupied,nu,nv,false), b = rectangles(occupied,nu,nv,true);
    if(b.size() < a.size()) a = std::move(b); // stable tie: runs along local v
    r.partition_ms += ms(t);
    if(a.size() > options.max_free_prisms-r.prisms.size())
    { r.sweep_ms = ms(start)-r.partition_ms; return fallback("free_prism_limit"); }
    for(const auto& q : a)
    {
      FreePrism3D prism;
      prism.lower[axis] = r.coordinates[axis][slab]; prism.upper[axis] = r.coordinates[axis][slab+1];
      prism.lower[u] = r.coordinates[u][q.y0]; prism.upper[u] = r.coordinates[u][q.y1];
      prism.lower[v] = r.coordinates[v][q.z0]; prism.upper[v] = r.coordinates[v][q.z1];
      r.prisms.push_back(prism);
    }
  }
  require(std::none_of(occupied.begin(),occupied.end(),[](auto v) { return v != 0; }),"Boundary sweep did not exit all obstacles");
  r.sweep_ms = ms(start)-r.partition_ms; r.prisms_before_merge = r.prisms.size();
  start = Clock::now();
  for(unsigned pass = 0; pass < options.max_merge_passes; ++pass)
  {
    const auto before = r.prisms.size();
    for(unsigned axis = 0; axis < 3; ++axis) mergeAxis(r.prisms,axis);
    ++r.merge_passes;
    if(r.prisms.size() == before) { r.merge_fixed_point = true; break; }
  }
  std::sort(r.prisms.begin(),r.prisms.end(),[](const auto& a,const auto& b) { return std::tie(a.lower,a.upper)<std::tie(b.lower,b.upper); });
  r.merge_ms = ms(start);
  start = Clock::now(); r.tree = closureTree(r,voxel_size,r.leaf_sources); r.tree_ms = ms(start);
  return r;
}

static BoundaryTreeBatch buildBatch(const std::function<BoundaryTreeResult(const BoundaryTreeOptions&)>& build,
                                    const BoundaryTreeOptions& options, bool best_axis)
{
  BoundaryTreeBatch batch;
  const auto start = Clock::now();
  for(unsigned i = 0; i < (best_axis ? 3u : 1u); ++i)
  {
    auto candidate_options = options;
    if(best_axis) candidate_options.sweep_axis = i;
    const auto t = Clock::now();
    batch.candidates.push_back(build(candidate_options));
    batch.candidate_build_ms.push_back(ms(t));
  }
  batch.all_candidates_ms = ms(start);
  const auto selection_start = Clock::now();
  const auto key = [](const BoundaryTreeResult& r) {
    return std::make_tuple(r.fallback,HS::serializeTree(r.tree.simplified).nodes.size(),r.tree.leaf_references,r.prisms.size(),r.sweep_axis);
  };
  if(best_axis)
  {
    auto best = key(batch.candidates.front());
    for(std::size_t i = 1; i < batch.candidates.size(); ++i)
    {
      const auto candidate = key(batch.candidates[i]);
      if(candidate < best) { batch.selected = i; best = candidate; }
    }
  }
  batch.selection_ms = ms(selection_start);
  return batch;
}

BoundaryTreeResult buildBoundaryTree(const OctreeOccupancy& o, const OctreeCellPartition& p,
                                    const OctreeBoundary& b, const BoundaryTreeOptions& options)
{
  return buildBoundaryTreeImpl(o.voxel_size,b,options,[&] {
    std::vector<HS::Vec3> vertices; auto clusters=octreeAabbClusters(p,vertices);
    auto tree=buildClusterTree(clusters,0,"octree_aabb");
    tree.simplified=HS::simplifyToFixedPoint(tree.raw); return tree;
  });
}
BoundaryTreeBatch buildBoundaryTreeBatch(const OctreeOccupancy& o, const OctreeCellPartition& p,
                                        const OctreeBoundary& b, const BoundaryTreeOptions& options, bool best)
{ return buildBatch([&](const auto& opts) { return buildBoundaryTree(o,p,b,opts); },options,best); }
BoundaryTreeBatch buildVoxelBoundaryTreeBatch(double h, const OctreeBoundary& b, const BoundaryTreeOptions& options,
                                             bool best, const std::function<TreeRepresentation()>& fallback)
{ return buildBatch([&](const auto& opts) { return buildBoundaryTreeImpl(h,b,opts,fallback); },options,best); }

static void validateTreeData(const std::vector<VoxelIndex>& occupied, double h,
                             const OctreeBoundary& boundary, const BoundaryTreeResult& r,
                             const TreeRepresentation* fallback_tree)
{
  require(r.sweep_axis < 3,"Invalid result sweep axis");
  BoundaryTreeResult expected; registry(boundary,expected);
  require(r.coordinates == expected.coordinates && r.supports.size() == expected.supports.size(),"Boundary support registry changed");
  for(std::size_t i = 0; i < r.supports.size(); ++i)
  {
    const auto& a = r.supports[i]; const auto& b = expected.supports[i];
    require(a.axis == b.axis && a.coordinate == b.coordinate && a.outward_signs == b.outward_signs && a.patches == b.patches,
            "Boundary support provenance changed");
  }
  if(r.fallback)
  {
    require(!r.fallback_reason.empty() && r.prisms.empty() && r.leaf_sources.empty(),"Invalid fallback metadata");
    require(fallback_tree != nullptr,"Missing independently validated fallback tree");
    checkTree(r.tree,*fallback_tree); return;
  }
  const auto nx = r.coordinates[0].size()-1, ny = r.coordinates[1].size()-1, nz = r.coordinates[2].size()-1;
  require(r.compressed_cells/ny/nz == nx && r.compressed_cells == nx*ny*nz,"Invalid compressed grid size");
  std::vector<std::size_t> count(r.compressed_cells,0);
  const auto index = [&](std::size_t x,std::size_t y,std::size_t z) { return (x*ny+y)*nz+z; };
  for(const auto& q : occupied)
  {
    std::size_t ix[3];
    for(unsigned k = 0; k < 3; ++k)
    {
      const auto& c = r.coordinates[k]; const auto i = std::upper_bound(c.begin(),c.end(),q[k]);
      require(i != c.begin() && i != c.end(),"Source voxel outside compressed bbox"); ix[k] = static_cast<std::size_t>(i-c.begin()-1);
    }
    ++count[index(ix[0],ix[1],ix[2])];
  }
  // Every occupied compressed open cell must be entirely occupied: cardinality
  // proves this without floating centers, dense fine grids, or overflowing volumes.
  for(std::size_t x = 0; x < nx; ++x) for(std::size_t y = 0; y < ny; ++y) for(std::size_t z = 0; z < nz; ++z)
  {
    const auto n = count[index(x,y,z)]; if(!n) continue;
    const std::size_t ix[3]{x,y,z}; std::size_t capacity = 1;
    for(unsigned k = 0; k < 3; ++k)
    {
      const auto width = static_cast<std::uint64_t>(r.coordinates[k][ix[k]+1]-r.coordinates[k][ix[k]]);
      require(width <= n/capacity,"Compressed cell is not homogeneous"); capacity *= static_cast<std::size_t>(width);
    }
    require(capacity == n,"Compressed occupied-cell volume mismatch");
  }
  std::vector<unsigned char> covered(r.compressed_cells,0);
  for(const auto& box : r.prisms)
  {
    std::size_t lo[3],hi[3];
    for(unsigned k = 0; k < 3; ++k)
    { lo[k] = at(r.coordinates[k],box.lower[k]); hi[k] = at(r.coordinates[k],box.upper[k]); require(lo[k] < hi[k],"Empty free prism"); }
    for(auto x = lo[0]; x < hi[0]; ++x) for(auto y = lo[1]; y < hi[1]; ++y) for(auto z = lo[2]; z < hi[2]; ++z)
    {
      const auto id = index(x,y,z);
      require(!count[id] && !covered[id],"Free prism overlaps occupied cell/another prism"); covered[id] = 1;
    }
  }
  for(std::size_t i = 0; i < count.size(); ++i) require(bool(covered[i]) == (count[i] == 0),"Incomplete free-space coverage");
  std::vector<std::pair<std::size_t,int>> sources;
  const auto tree = closureTree(r,h,sources);
  require(sources == r.leaf_sources,"Tree support polarity/provenance changed"); checkTree(r.tree,tree);
}

void validateBoundaryTree(const OctreeOccupancy& o, const OctreeCellPartition& p,
                          const OctreeBoundary& b, const BoundaryTreeResult& r)
{
  validateOctreeCellPartition(o,p); validateOctreeBoundary(o,p,b);
  TreeRepresentation fallback;
  if(r.fallback)
  {
    std::vector<HS::Vec3> vertices; auto clusters=octreeAabbClusters(p,vertices);
    fallback=buildClusterTree(clusters,0,"octree_aabb");
    fallback.simplified=HS::simplifyToFixedPoint(fallback.raw);
  }
  validateTreeData(o.occupied_voxels,o.voxel_size,b,r,r.fallback ? &fallback : nullptr);
}
void validateVoxelBoundaryTree(const SparseVoxelOccupancy& o, const OctreeBoundary& b, const BoundaryTreeResult& r)
{
  require(!r.fallback,"Voxel fallback requires the legacy partition certificate");
  validateVoxelBoundary(o.occupied_voxels,b);
  validateTreeData(o.occupied_voxels,o.voxel_size,b,r,nullptr);
}
template<class Occupancy>
static bool strictFree(const BoundaryTreeResult& r, const Occupancy& octree, const HS::Vec3& p, bool* used_zero_guard)
{
  require(std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z),"Nonfinite free-space query");
  if(used_zero_guard) *used_zero_guard = false;
  const double value = HS::evaluate(r.tree.simplified,r.tree.leaves,p);
  if(r.fallback) return value > 0;
  if(value < 0) return true;
  if(value > 0) return false;
  if(used_zero_guard) *used_zero_guard = true;
  return !octree.contains(p);
}
bool boundaryStrictFree(const BoundaryTreeResult& r,const OctreeOccupancy& o,const HS::Vec3& p,bool* zero)
{ return strictFree(r,o,p,zero); }
bool boundaryStrictFree(const BoundaryTreeResult& r,const SparseVoxelOccupancy& o,const HS::Vec3& p,bool* zero)
{ return strictFree(r,o,p,zero); }
} // namespace rokae_demo
