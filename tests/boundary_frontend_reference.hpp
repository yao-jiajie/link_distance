#pragma once
// Frozen pre-optimization ordered-map implementation, only for differential
// tests and same-process benchmarks. Keep independent of production merging.
#include <rokae_demo/octree_boundary_tree.hpp>
#include <algorithm>
#include <chrono>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>

namespace boundary_frontend_reference
{
using namespace rokae_demo;
using Clock = std::chrono::steady_clock;
using PlaneKey = std::tuple<unsigned, int, std::int64_t>;
inline double elapsed(Clock::time_point start)
{
  return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

inline void require(bool valid, const char* message)
{
  if(!valid) throw std::runtime_error(message);
}

inline void checkInput(const std::vector<VoxelIndex>& occupied, const std::vector<std::size_t>* owners, std::size_t cells)
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

inline OctreeBoundary extractBoundary(const std::vector<VoxelIndex>& occupied,
                                     const std::vector<std::size_t>* owners, std::size_t cells)
{
  const auto start = Clock::now();
  checkInput(occupied, owners, cells);
  OctreeBoundary result;
  result.voxel_faces_before = 6 * occupied.size();
  for(std::size_t i = 0; i < occupied.size(); ++i)
    for(unsigned axis = 0; axis < 3; ++axis)
      for(const int sign : {-1, 1})
      {
        auto neighbor = occupied[i];
        neighbor[axis] += sign;
        if(std::binary_search(occupied.begin(), occupied.end(), neighbor))
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
  // Sparse row runs on each oriented plane. Extend only an identical interval
  // in the immediately previous row: never span holes, gaps, or opposite normals.
  // O(F log F) time/O(F) storage, without scanning the (possibly huge) bbox or
  // comparing all pairs. Deterministic rectangular tiling, not a minimum tiling.
  using RowPosition = std::pair<std::int64_t, std::int64_t>; // v, u
  std::map<PlaneKey, std::map<RowPosition, std::size_t>> planes;
  for(std::size_t i = 0; i < result.faces.size(); ++i)
  {
    const auto& f = result.faces[i];
    planes[{f.axis, f.sign, f.plane}].emplace(RowPosition{f.v, f.u}, i);
  }
  for(const auto& plane : planes)
  {
    using Interval = std::pair<std::int64_t, std::int64_t>;
    std::map<Interval, std::size_t> previous;
    auto it = plane.second.begin();
    while(it != plane.second.end())
    {
      const auto row = it->first.first;
      std::map<Interval, std::size_t> current;
      while(it != plane.second.end() && it->first.first == row)
      {
        const auto u0 = it->first.second;
        auto u1 = u0;
        std::vector<std::size_t> ids;
        do
        {
          ids.push_back(it->second);
          ++u1; ++it;
        } while(it != plane.second.end() && it->first.first == row && it->first.second == u1);
        const Interval interval{u0, u1};
        const auto old = previous.find(interval);
        std::size_t id;
        if(old != previous.end() && result.patches[old->second].v1 == row)
        {
          id = old->second;
          auto& patch = result.patches[id];
          patch.v1 = row + 1;
          patch.source_face_ids.insert(patch.source_face_ids.end(), ids.begin(), ids.end());
        }
        else
        {
          id = result.patches.size();
          result.patches.push_back({std::get<0>(plane.first), std::get<1>(plane.first), std::get<2>(plane.first),
                                    u0, u1, row, row + 1, std::move(ids)});
        }
        current.emplace(interval, id);
      }
      previous = std::move(current);
    }
  }
  result.merge_ms = elapsed(merge_start);
  return result;
}

inline BoundaryPlaneRegistry makeRegistry(const OctreeBoundary& b)
{
  BoundaryPlaneRegistry r;
  std::map<std::pair<unsigned,std::int64_t>,BoundarySupportPlane> planes;
  for(std::size_t i = 0; i < b.patches.size(); ++i)
  {
    const auto& p = b.patches[i]; auto& s = planes[{p.axis,p.plane}];
    s.axis = p.axis; s.coordinate = p.plane; s.outward_signs |= p.sign < 0 ? 1 : 2;
    s.patches.push_back(i);
  }
  for(const auto& item : planes)
  { r.supports.push_back(item.second); r.coordinates[item.first.first].push_back(item.first.second); }
  for(const auto& c : r.coordinates) require(c.size() >= 2,"Missing boundary extent");
  return r;
}
} // namespace boundary_frontend_reference
