#pragma once
#include <rokae_demo/octree_boundary.hpp>
#include <rokae_demo/sparse_voxel_occupancy.hpp>

namespace rokae_demo
{
// Independent Euclidean distance to FINITE exposed rectangular patches.
// Validation/reference only, not used inside the smooth function or its gradient.
// Caller must certify boundary against the SAME occupancy snapshot first.
class BoundaryDistanceReference
{
  std::vector<OctreeBox> patches_;
  double h_;
public:
  BoundaryDistanceReference(double h,const OctreeBoundary& boundary): h_(h)
  {
    if(!std::isfinite(h) || h<=0 || boundary.patches.empty()) throw std::invalid_argument("Invalid boundary distance input");
    patches_.reserve(boundary.patches.size());
    for(const auto& p:boundary.patches)
    {
      if(p.axis>=3 || p.u0>=p.u1 || p.v0>=p.v1) throw std::invalid_argument("Invalid finite boundary patch");
      const auto corners=octreeBoundaryCorners(p); OctreeBox b;
      for(unsigned a=0;a<3;++a)
      {
        b.lower[a]=static_cast<double>(std::min(corners[0][a],corners[2][a]))*h;
        b.upper[a]=static_cast<double>(std::max(corners[0][a],corners[2][a]))*h;
        if(!std::isfinite(b.lower[a]) || !std::isfinite(b.upper[a])) throw std::invalid_argument("Nonfinite patch coordinates");
      }
      patches_.push_back(b);
    }
  }
  std::size_t patchCount() const { return patches_.size(); }
  double unsignedDistance(const HS::Vec3& p) const
  {
    if(!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z)) throw std::invalid_argument("Nonfinite distance query");
    const double x[3]{p.x,p.y,p.z}; double result=INFINITY;
    for(const auto& b:patches_)
    {
      double d[3]{};
      for(unsigned a=0;a<3;++a) d[a]=x[a]<b.lower[a] ? b.lower[a]-x[a] : x[a]>b.upper[a] ? x[a]-b.upper[a] : 0;
      result=std::min(result,std::hypot(d[0],d[1],d[2]));
    }
    if(!std::isfinite(result)) throw std::overflow_error("Unrepresentable boundary distance");
    return result;
  }
  double signedDistance(const HS::Vec3& p,const SparseVoxelOccupancy& occupancy) const
  {
    if(occupancy.voxel_size!=h_) throw std::invalid_argument("Distance reference voxel size mismatch");
    const double distance=unsignedDistance(p);
    return distance==0 ? 0 : occupancy.contains(p) ? -distance : distance;
  }
};
} // namespace rokae_demo
