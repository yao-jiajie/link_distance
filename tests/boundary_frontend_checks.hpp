#pragma once
#include "boundary_frontend_reference.hpp"

namespace boundary_frontend_checks
{
using namespace rokae_demo;
inline void require(bool ok, const char* why) { if(!ok) throw std::runtime_error(why); }
inline void sameBoundary(const OctreeBoundary& a, const OctreeBoundary& b)
{
  require(a.voxel_faces_before == b.voxel_faces_before && a.internal_faces_removed == b.internal_faces_removed &&
          a.faces.size() == b.faces.size() && a.patches.size() == b.patches.size(), "Frontend face/patch counts changed");
  for(std::size_t i = 0; i < a.faces.size(); ++i)
  {
    const auto& x = a.faces[i]; const auto& y = b.faces[i];
    require(std::tie(x.axis,x.sign,x.plane,x.u,x.v,x.voxel_id,x.convex_cell_id) ==
            std::tie(y.axis,y.sign,y.plane,y.u,y.v,y.voxel_id,y.convex_cell_id), "Frontend face order/owner changed");
  }
  for(std::size_t i = 0; i < a.patches.size(); ++i)
  {
    const auto& x = a.patches[i]; const auto& y = b.patches[i];
    require(std::tie(x.axis,x.sign,x.plane,x.u0,x.u1,x.v0,x.v1,x.source_face_ids) ==
            std::tie(y.axis,y.sign,y.plane,y.u0,y.u1,y.v0,y.v1,y.source_face_ids), "Frontend patch tiling/provenance changed");
    require(octreeBoundaryCorners(x) == octreeBoundaryCorners(y), "Frontend oriented corners changed");
  }
}
inline void sameRegistry(const BoundaryPlaneRegistry& a, const BoundaryPlaneRegistry& b)
{
  require(a.coordinates == b.coordinates && a.supports.size() == b.supports.size(), "Frontend support coordinates changed");
  for(std::size_t i = 0; i < a.supports.size(); ++i)
  {
    const auto& x = a.supports[i]; const auto& y = b.supports[i];
    require(std::tie(x.axis,x.coordinate,x.outward_signs,x.patches) ==
            std::tie(y.axis,y.coordinate,y.outward_signs,y.patches), "Frontend support order/provenance changed");
  }
}
} // namespace boundary_frontend_checks
