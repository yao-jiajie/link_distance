#include <rokae_demo/octree_boundary.hpp>
#include <rokae_demo/sparse_voxel_occupancy.hpp>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>
#ifdef ROKAE_USE_ABSEIL_BOUNDARY_HASH
#include <absl/container/flat_hash_set.h>
#else
#include <unordered_set>
#endif

namespace rokae_demo
{
namespace
{
using Clock = std::chrono::steady_clock;
using GridFaceKey = std::tuple<unsigned, int, std::int64_t, std::int64_t, std::int64_t>;

struct VoxelIndexHash
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

double elapsed(Clock::time_point start)
{
  return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

void require(bool valid, const char* message)
{
  if(!valid) throw std::runtime_error(message);
}

void checkInput(const std::vector<VoxelIndex>& occupied, const std::vector<std::size_t>* owners, std::size_t cells)
{
  require(!occupied.empty() &&
          (!owners || owners->size() == occupied.size()), "Boundary input size mismatch");
  require(occupied.size() <= std::numeric_limits<std::size_t>::max() / 6,
          "Boundary face count overflow");
  for(std::size_t i = 0; i < occupied.size(); ++i)
  {
    require(!i || occupied[i-1] < occupied[i], "Boundary occupancy must be sorted/unique");
    require(!owners || (*owners)[i] < cells, "Boundary cell ID out of range");
    for(const auto c : occupied[i])
      require(c > std::numeric_limits<std::int64_t>::min() && c < std::numeric_limits<std::int64_t>::max(),
              "Boundary neighbor index overflow");
  }
}

GridFaceKey key(const OctreeBoundaryFace& f) { return {f.axis, f.sign, f.plane, f.u, f.v}; }

OctreeBoundaryPatch unitPatch(const OctreeBoundaryFace& f)
{
  return {f.axis, f.sign, f.plane, f.u, f.u + 1, f.v, f.v + 1, {}};
}

std::ofstream outputFile(const std::filesystem::path& path)
{
  std::ofstream out(path);
  if(!out) throw std::runtime_error("Cannot write boundary output: " + path.string());
  out.exceptions(std::ios::badbit | std::ios::failbit);
  out << std::setprecision(17);
  return out;
}

template<class Range> void arrayJson(std::ostream& out, const Range& values)
{
  out << '[';
  bool first = true;
  for(const auto& value : values)
  {
    if(!first) out << ',';
    first = false;
    out << value;
  }
  out << ']';
}
}  // namespace

static OctreeBoundary extractBoundary(const std::vector<VoxelIndex>& occupied,
                                     const std::vector<std::size_t>* owners, std::size_t cells)
{
  const auto start = Clock::now();
  checkInput(occupied, owners, cells);
  OctreeBoundary result;
  result.voxel_faces_before = 6 * occupied.size();
  result.faces.reserve(result.voxel_faces_before);
  VoxelIndex lower=occupied.front(),upper=occupied.front();
  for(const auto& voxel:occupied)
    for(unsigned axis=0;axis<3;++axis)
    {
      lower[axis]=std::min(lower[axis],voxel[axis]);
      upper[axis]=std::max(upper[axis],voxel[axis]);
    }
  std::array<std::size_t,3> extent{};
  std::size_t dense_cells=1;
  bool use_dense=true;
  for(unsigned axis=0;axis<3;++axis)
  {
    const auto width=static_cast<std::uint64_t>(upper[axis])-static_cast<std::uint64_t>(lower[axis])+1;
    if(width>std::numeric_limits<std::size_t>::max() ||
       static_cast<std::size_t>(width)>std::numeric_limits<std::size_t>::max()/dense_cells)
    { use_dense=false; break; }
    extent[axis]=static_cast<std::size_t>(width); dense_cells*=extent[axis];
  }
  constexpr std::size_t absolute_dense_limit=16*1024*1024;
  const auto relative_dense_limit=occupied.size()<=std::numeric_limits<std::size_t>::max()/64
    ? std::max(std::size_t{4096},occupied.size()*64) : std::numeric_limits<std::size_t>::max();
  use_dense=use_dense && dense_cells<=absolute_dense_limit && dense_cells<=relative_dense_limit;
  const auto denseIndex=[&](const VoxelIndex& voxel) {
    const auto x=static_cast<std::uint64_t>(voxel[0])-static_cast<std::uint64_t>(lower[0]);
    const auto y=static_cast<std::uint64_t>(voxel[1])-static_cast<std::uint64_t>(lower[1]);
    const auto z=static_cast<std::uint64_t>(voxel[2])-static_cast<std::uint64_t>(lower[2]);
    return (static_cast<std::size_t>(x)*extent[1]+static_cast<std::size_t>(y))*extent[2]+
           static_cast<std::size_t>(z);
  };
  std::vector<unsigned char> dense_lookup;
#ifdef ROKAE_USE_ABSEIL_BOUNDARY_HASH
  absl::flat_hash_set<VoxelIndex,VoxelIndexHash> occupied_lookup;
#else
  std::unordered_set<VoxelIndex,VoxelIndexHash> occupied_lookup;
#endif
  if(use_dense)
  {
    dense_lookup.assign(dense_cells,0);
    for(const auto& voxel:occupied) dense_lookup[denseIndex(voxel)]=1;
  }
  else
  {
    occupied_lookup.reserve(occupied.size()*2);
    occupied_lookup.insert(occupied.begin(),occupied.end());
  }
  const auto isOccupied=[&](const VoxelIndex& voxel) {
    if(!use_dense) return occupied_lookup.find(voxel)!=occupied_lookup.end();
    for(unsigned axis=0;axis<3;++axis)
      if(voxel[axis]<lower[axis] || voxel[axis]>upper[axis]) return false;
    return dense_lookup[denseIndex(voxel)]!=0;
  };
  for(std::size_t i = 0; i < occupied.size(); ++i)
    for(unsigned axis = 0; axis < 3; ++axis)
      for(const int sign : {-1, 1})
      {
        auto neighbor = occupied[i];
        neighbor[axis] += sign;
        if(isOccupied(neighbor))
        {
          ++result.internal_faces_removed; // BOTH directed sides are removed.
          continue;
        }
        result.faces.push_back({axis, sign, occupied[i][axis] + (sign > 0),
                                occupied[i][(axis+1)%3], occupied[i][(axis+2)%3],
                                i, owners ? (*owners)[i] : i});
      }
  result.extraction_ms = elapsed(start);
  const auto merge_start = Clock::now();
  // Keep public faces in place. Compact sort records avoid repeated indirect
  // face loads, and the six orientation buckets avoid sorting axis/sign keys.
  // Their iteration order matches the former nested ordered maps exactly.
  struct Position { std::int64_t plane, v, u; std::size_t face; };
  struct PackedPosition { std::uint64_t key; std::size_t face; };
  struct SideEncoding
  {
    std::int64_t min_plane=0,min_v=0,min_u=0,max_plane=0,max_v=0,max_u=0;
    unsigned plane_bits=0,v_bits=0,u_bits=0;
    bool initialized=false,packed=false;
  };
  std::array<std::vector<Position>,6> expanded_sides;
  std::array<std::vector<PackedPosition>,6> packed_sides;
  std::array<SideEncoding,6> encodings;
  std::array<std::size_t,6> counts{};
  for(const auto& f : result.faces)
  {
    const auto side=2*f.axis+(f.sign>0); ++counts[side]; auto& e=encodings[side];
    if(!e.initialized)
    {
      e.min_plane=e.max_plane=f.plane; e.min_v=e.max_v=f.v; e.min_u=e.max_u=f.u;
      e.initialized=true;
    }
    else
    {
      e.min_plane=std::min(e.min_plane,f.plane); e.max_plane=std::max(e.max_plane,f.plane);
      e.min_v=std::min(e.min_v,f.v); e.max_v=std::max(e.max_v,f.v);
      e.min_u=std::min(e.min_u,f.u); e.max_u=std::max(e.max_u,f.u);
    }
  }
  const auto bits=[](std::uint64_t value) {
    unsigned result=0; while(value) { ++result; value>>=1; } return result;
  };
  for(unsigned side=0;side<6;++side)
  {
    auto& e=encodings[side]; if(!e.initialized) continue;
    e.plane_bits=bits(static_cast<std::uint64_t>(e.max_plane)-static_cast<std::uint64_t>(e.min_plane));
    e.v_bits=bits(static_cast<std::uint64_t>(e.max_v)-static_cast<std::uint64_t>(e.min_v));
    e.u_bits=bits(static_cast<std::uint64_t>(e.max_u)-static_cast<std::uint64_t>(e.min_u));
    // Keeping the total below 64 also makes every decoded signed delta safe.
    e.packed=e.plane_bits+e.v_bits+e.u_bits<=63;
    if(e.packed) packed_sides[side].reserve(counts[side]);
    else expanded_sides[side].reserve(counts[side]);
  }
  for(std::size_t id = 0; id < result.faces.size(); ++id)
  {
    const auto& f=result.faces[id]; const auto side=2*f.axis+(f.sign>0); const auto& e=encodings[side];
    if(e.packed)
    {
      const auto plane=static_cast<std::uint64_t>(f.plane)-static_cast<std::uint64_t>(e.min_plane);
      const auto v=static_cast<std::uint64_t>(f.v)-static_cast<std::uint64_t>(e.min_v);
      const auto u=static_cast<std::uint64_t>(f.u)-static_cast<std::uint64_t>(e.min_u);
      packed_sides[side].push_back({((plane<<e.v_bits)|v)<<e.u_bits|u,id});
    }
    else expanded_sides[side].push_back({f.plane,f.v,f.u,id});
  }
  struct Run { std::int64_t u0, u1; std::size_t patch; };
  std::vector<Run> previous, current;
  const auto mergeSide=[&](unsigned side,auto& order,auto planeOf,auto vOf,auto uOf) {
    for(std::size_t begin = 0; begin < order.size();)
    {
      const auto plane=planeOf(order[begin]);
      auto end = begin + 1;
      while(end<order.size() && planeOf(order[end])==plane) ++end;
      previous.clear();
      while(begin < end)
      {
        const auto row=vOf(order[begin]);
        current.clear();
        std::size_t old = 0;
        while(begin<end && vOf(order[begin])==row)
        {
          const auto first = begin;
          const auto u0=uOf(order[begin]);
          auto u1 = u0;
          do { ++u1; ++begin; }
          while(begin<end && vOf(order[begin])==row && uOf(order[begin])==u1);
          // Both rows' disjoint runs are ordered. Merge-join them in linear time;
          // only an identical interval in the immediately adjacent row extends.
          while(old < previous.size() && std::tie(previous[old].u0, previous[old].u1) < std::tie(u0, u1)) ++old;
          std::size_t id;
          if(old < previous.size() && previous[old].u0 == u0 && previous[old].u1 == u1 &&
             result.patches[previous[old].patch].v1 == row)
          {
            id = previous[old].patch;
            result.patches[id].v1 = row + 1;
          }
          else
          {
            id = result.patches.size();
            result.patches.push_back({side/2, side%2 ? 1 : -1, plane, u0, u1, row, row + 1, {}});
          }
          auto& ids = result.patches[id].source_face_ids;
          if(ids.empty()) ids.reserve(begin-first);
          for(auto i = first; i < begin; ++i) ids.push_back(order[i].face);
          current.push_back({u0, u1, id});
        }
        previous.swap(current); // reuse both row buffers across rows and planes
      }
    }
  };
  for(unsigned side = 0; side < 6; ++side)
  {
    const auto& encoding=encodings[side];
    if(encoding.packed)
    {
      auto& order=packed_sides[side];
      std::sort(order.begin(),order.end(),[](const auto& a,const auto& b) { return a.key<b.key; });
      const auto u_mask=encoding.u_bits ? (std::uint64_t{1}<<encoding.u_bits)-1 : 0;
      const auto v_mask=encoding.v_bits ? (std::uint64_t{1}<<encoding.v_bits)-1 : 0;
      mergeSide(side,order,
        [&](const auto& position) { return encoding.min_plane+static_cast<std::int64_t>(position.key>>(encoding.v_bits+encoding.u_bits)); },
        [&](const auto& position) { return encoding.min_v+static_cast<std::int64_t>((position.key>>encoding.u_bits)&v_mask); },
        [&](const auto& position) { return encoding.min_u+static_cast<std::int64_t>(position.key&u_mask); });
    }
    else
    {
      auto& order=expanded_sides[side];
      std::sort(order.begin(),order.end(),[](const auto& a,const auto& b) {
        return std::tie(a.plane,a.v,a.u)<std::tie(b.plane,b.v,b.u);
      });
      mergeSide(side,order,[](const auto& position) { return position.plane; },
        [](const auto& position) { return position.v; },[](const auto& position) { return position.u; });
    }
  }
  result.merge_ms = elapsed(merge_start);
  return result;
}

OctreeBoundary extractOctreeBoundary(const OctreeOccupancy& o, const OctreeCellPartition& p)
{ return extractBoundary(o.occupied_voxels, &p.voxel_cell_ids, p.cells.size()); }
OctreeBoundary extractVoxelBoundary(const std::vector<VoxelIndex>& keys)
{ return extractBoundary(keys, nullptr, keys.size()); }
OctreeBoundary extractSparseVoxelBoundary(const SparseVoxelOccupancy& occupancy)
{ return extractBoundary(occupancy.occupied_voxels,nullptr,occupancy.occupied_voxels.size()); }


std::array<VoxelIndex, 4> octreeBoundaryCorners(const OctreeBoundaryPatch& patch)
{
  require(patch.axis < 3 && (patch.sign == -1 || patch.sign == 1) &&
          patch.u0 < patch.u1 && patch.v0 < patch.v1, "Invalid boundary patch geometry");
  std::array<VoxelIndex, 4> corners{};
  for(unsigned i = 0; i < 4; ++i)
  {
    corners[i][patch.axis] = patch.plane;
    corners[i][(patch.axis+1)%3] = i == 1 || i == 2 ? patch.u1 : patch.u0;
    corners[i][(patch.axis+2)%3] = i >= 2 ? patch.v1 : patch.v0;
  }
  if(patch.sign < 0) std::swap(corners[1], corners[3]);
  return corners;
}

static void validateBoundary(const std::vector<VoxelIndex>& occupied, const std::vector<std::size_t>* owners,
                             std::size_t cells, const OctreeBoundary& boundary)
{
  checkInput(occupied, owners, cells);
  // Independently reconstruct the expected oriented faces and their owners.
  std::map<GridFaceKey, std::size_t> expected;
  for(std::size_t i = 0; i < occupied.size(); ++i)
    for(unsigned side = 0; side < 6; ++side)
    {
      const auto axis = side / 2;
      const int sign = side % 2 ? 1 : -1;
      auto other = occupied[i]; other[axis] += sign;
      if(std::binary_search(occupied.begin(), occupied.end(), other)) continue;
      expected.emplace(GridFaceKey{axis, sign, occupied[i][axis] + (sign > 0),
                              occupied[i][(axis+1)%3], occupied[i][(axis+2)%3]}, i);
    }
  require(boundary.voxel_faces_before == 6 * occupied.size() &&
          boundary.faces.size() == expected.size() &&
          boundary.internal_faces_removed == boundary.voxel_faces_before - expected.size(),
          "Boundary face count certificate failed");
  for(const auto& f : boundary.faces)
  {
    const auto found = expected.find(key(f));
    require(found != expected.end() && f.voxel_id == found->second &&
            f.convex_cell_id == (owners ? (*owners)[found->second] : found->second), "Missing, duplicate or invalid boundary face");
    expected.erase(found);
  }
  require(expected.empty(), "Boundary misses exposed faces");
  std::vector<unsigned char> covered(boundary.faces.size(), 0);
  for(const auto& p : boundary.patches)
  {
    require(p.axis < 3 && (p.sign == 1 || p.sign == -1) && p.u0 < p.u1 && p.v0 < p.v1,
            "Invalid boundary rectangle");
    // Unsigned subtraction avoids signed overflow, even for malicious metadata.
    const auto width = static_cast<std::uint64_t>(p.u1) - static_cast<std::uint64_t>(p.u0);
    const auto height = static_cast<std::uint64_t>(p.v1) - static_cast<std::uint64_t>(p.v0);
    require(!p.source_face_ids.empty() && width <= p.source_face_ids.size() &&
            p.source_face_ids.size() % width == 0 && p.source_face_ids.size() / width == height,
            "Boundary patch is not a full rectangle");
    for(const auto id : p.source_face_ids)
    {
      require(id < covered.size() && !covered[id], "Duplicate or invalid patch source face");
      covered[id] = 1;
      const auto& f = boundary.faces[id];
      require(f.axis == p.axis && f.sign == p.sign && f.plane == p.plane &&
              p.u0 <= f.u && f.u < p.u1 && p.v0 <= f.v && f.v < p.v1,
              "Patch crosses a gap or changes orientation");
    }
  }
  require(std::all_of(covered.begin(), covered.end(), [](unsigned char c) { return c == 1; }),
          "Boundary patch tiling is incomplete");
}

void validateOctreeBoundary(const OctreeOccupancy& o, const OctreeCellPartition& p, const OctreeBoundary& b)
{ validateBoundary(o.occupied_voxels, &p.voxel_cell_ids, p.cells.size(), b); }
void validateVoxelBoundary(const std::vector<VoxelIndex>& keys, const OctreeBoundary& b)
{ validateBoundary(keys, nullptr, keys.size(), b); }

static void writeBoundary(const std::filesystem::path& output, double h,
                          const OctreeBoundary& boundary, bool voxel_owners)
{
  auto json = outputFile(output / "boundary_patches.json");
  json << "{\"pipeline\":\"" << (voxel_owners ? "voxels" : "octree")
       << (voxel_owners ? "\",\"face_owner_kind\":\"voxel" : "")
       << "\",\"phase\":3,\"length_unit\":\"m\",\"voxel_size\":" << h
       << ",\"representation\":\"exact_oriented_rectangular_boundary\",\"is_cluster_support_set\":false,"
          "\"mesh_manifold_guaranteed\":false,\"voxel_faces_before\":" << boundary.voxel_faces_before
       << ",\"internal_faces_removed\":" << boundary.internal_faces_removed << ",\"faces\":[\n";
  for(std::size_t i = 0; i < boundary.faces.size(); ++i)
  {
    const auto& f = boundary.faces[i];
    json << (i ? ",\n" : "") << "{\"id\":" << i << ",\"axis\":" << f.axis << ",\"sign\":" << f.sign
         << ",\"plane_index\":" << f.plane << ",\"u\":" << f.u << ",\"v\":" << f.v
         << ",\"source_voxel_id\":" << f.voxel_id << ",\"convex_cell_id\":" << f.convex_cell_id << '}';
  }
  json << "\n],\"patches\":[\n";
  for(std::size_t i = 0; i < boundary.patches.size(); ++i)
  {
    const auto& p = boundary.patches[i];
    const auto corners = octreeBoundaryCorners(p);
    std::array<int,3> normal{}; normal[p.axis] = p.sign;
    json << (i ? ",\n" : "") << "{\"id\":" << i << ",\"axis\":" << p.axis << ",\"sign\":" << p.sign
         << ",\"plane_index\":" << p.plane << ",\"uv_bounds\": [" << p.u0 << ',' << p.u1 << ',' << p.v0 << ',' << p.v1
         << "],\"normal\":"; arrayJson(json, normal);
    json << ",\"offset\":" << p.sign * (static_cast<double>(p.plane) * h)
         << ",\"corner_indices\":[";
    for(unsigned j = 0; j < 4; ++j) { if(j) json << ','; arrayJson(json, corners[j]); }
    json << "],\"corners_m\":[";
    for(unsigned j = 0; j < 4; ++j)
    {
      if(j) json << ',';
      std::array<double,3> point{};
      for(unsigned axis = 0; axis < 3; ++axis) point[axis] = static_cast<double>(corners[j][axis]) * h;
      arrayJson(json, point);
    }
    json << "],\"source_face_ids\":"; arrayJson(json, p.source_face_ids);
    std::set<std::size_t> owners;
    for(const auto face : p.source_face_ids) owners.insert(boundary.faces[face].convex_cell_id);
    json << ",\"convex_cell_ids\":"; arrayJson(json, owners);
    json << '}';
  }
  json << "\n]}\n";
  json.close();

  std::map<VoxelIndex, std::size_t> ids;
  std::vector<VoxelIndex> vertices;
  std::vector<std::array<std::size_t,4>> quads;
  for(const auto& f : boundary.faces)
  {
    const auto corners = octreeBoundaryCorners(unitPatch(f));
    std::array<std::size_t,4> quad{};
    for(unsigned j = 0; j < 4; ++j)
    {
      const auto added = ids.emplace(corners[j], vertices.size());
      if(added.second) vertices.push_back(corners[j]);
      quad[j] = added.first->second;
    }
    quads.push_back(quad);
  }
  require(quads.size() <= std::numeric_limits<std::size_t>::max()/2, "Boundary triangle count overflow");
  auto off = outputFile(output / "boundary.off");
  off << "OFF\n" << vertices.size() << ' ' << 2 * quads.size() << " 0\n";
  for(const auto& v : vertices)
    off << static_cast<double>(v[0])*h << ' ' << static_cast<double>(v[1])*h
        << ' ' << static_cast<double>(v[2])*h << '\n';
  for(const auto& q : quads)
    off << "3 " << q[0] << ' ' << q[1] << ' ' << q[2] << "\n3 " << q[0] << ' ' << q[2] << ' ' << q[3] << '\n';
  off.close();
}
void writeOctreeBoundary(const std::filesystem::path& p, const OctreeOccupancy& o, const OctreeBoundary& b)
{ writeBoundary(p, o.voxel_size, b, false); }
void writeVoxelBoundary(const std::filesystem::path& p, double h, const OctreeBoundary& b)
{ writeBoundary(p, h, b, true); }
}  // namespace rokae_demo
