#include <rokae_demo/boundary_direct_frame_cache.hpp>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace rokae_demo;
namespace
{
using Clock=std::chrono::steady_clock;
double elapsed(Clock::time_point start)
{ return std::chrono::duration<double,std::milli>(Clock::now()-start).count(); }
double median(std::vector<double> values)
{ std::sort(values.begin(),values.end()); return values[values.size()/2]; }
void require(bool ok,const char* why) { if(!ok) throw std::runtime_error(why); }

std::vector<VoxelPoint> readCloud(const char* path)
{
  std::ifstream in(path); require(bool(in),"Cannot read benchmark XYZ");
  std::vector<VoxelPoint> points; std::string line;
  while(std::getline(in,line))
  {
    const auto comment=line.find('#'); if(comment!=std::string::npos) line.resize(comment);
    std::istringstream row(line); row>>std::ws; if(row.eof()) continue;
    double x,y,z; require(bool(row>>x>>y>>z),"Invalid benchmark XYZ"); row>>std::ws;
    require(row.eof() && std::isfinite(x) && std::isfinite(y) && std::isfinite(z),"Invalid benchmark point");
    points.emplace_back(x,y,z);
  }
  require(!in.bad() && !points.empty(),"Empty benchmark XYZ"); return points;
}

std::vector<VoxelPoint> moveWithinVoxels(const std::vector<VoxelPoint>& raw,double h)
{
  std::vector<VoxelPoint> result; result.reserve(raw.size());
  for(const auto& point:raw)
  {
    const auto key=voxelIndex(point,h); double p[3]{};
    for(unsigned a=0;a<3;++a)
    {
      const double lo=static_cast<double>(key[a])*h,hi=static_cast<double>(key[a]+1)*h;
      p[a]=lo+(hi-lo)/2;
    }
    result.emplace_back(p[0],p[1],p[2]);
  }
  return result;
}
} // namespace

int main(int argc,char** argv)
{
  try
  {
    require(argc==4,"Usage: boundary_direct_frame_benchmark INPUT.xyz VOXEL_SIZE REPEATS");
    std::size_t used=0; const double h=std::stod(argv[2],&used);
    require(used==std::string(argv[2]).size() && std::isfinite(h) && h>0,"Invalid voxel size");
    used=0; const auto repeats=std::stoull(argv[3],&used);
    require(used==std::string(argv[3]).size() && repeats>0,"Invalid repeat count");
    const auto raw=readCloud(argv[1]),moved=moveWithinVoxels(raw,h);
    require(buildSparseVoxelOccupancy(raw,h).occupied_voxels==buildSparseVoxelOccupancy(moved,h).occupied_voxels,
            "Equivalent benchmark frame changed occupancy");

    BoundaryDirectFrameOptions options; options.smooth=true; options.smoothing={0,.001,true,true};
    options.cache_validation_values=true;
    std::vector<double> cold_times,hit_times,cold_build,hit_occupancy;
    cold_times.reserve(repeats); hit_times.reserve(repeats);
    double expected=0; std::size_t nodes=0;
    for(std::size_t i=0;i<repeats;++i)
    {
      BoundaryDirectFrameCache cold; const auto start=Clock::now();
      auto frame=cold.prepare(i%2 ? raw : moved,h,options); cold_times.push_back(elapsed(start));
      require(!frame.cache_hit && frame.snapshot->diagnostics().strict_sign_certified,"Cold frame was not certified");
      cold_build.push_back(frame.snapshot_build_ms); auto& snapshot=*frame.snapshot;
      const auto p=HS::Vec3{raw[0].x(),raw[0].y(),raw[0].z()};
      const auto value=snapshot.compiledProgram()->evaluate(p,snapshot.queryWorkspace());
      if(!i) { expected=value; nodes=snapshot.compiledProgram()->nodeCount(); }
      else require(value==expected && snapshot.compiledProgram()->nodeCount()==nodes,"Cold snapshots changed function");
    }

    BoundaryDirectFrameCache cache; auto seed=cache.prepare(raw,h,options);
    require(!seed.cache_hit,"Cache seed unexpectedly hit"); const auto snapshot=seed.snapshot;
    auto* query_data=snapshot->queryWorkspace().data(); auto* smooth_data=snapshot->smoothWorkspace().data();
    for(std::size_t i=0;i<repeats;++i)
    {
      const auto start=Clock::now(); auto frame=cache.prepare(i%2 ? raw : moved,h,options);
      hit_times.push_back(elapsed(start)); hit_occupancy.push_back(frame.occupancy_ms);
      require(frame.cache_hit && frame.snapshot==snapshot,"Equivalent frame missed reusable snapshot");
      require(frame.snapshot->queryWorkspace().data()==query_data && frame.snapshot->smoothWorkspace().data()==smooth_data,
              "Frame hit replaced a workspace");
      const auto p=HS::Vec3{raw[0].x(),raw[0].y(),raw[0].z()};
      require(frame.snapshot->compiledProgram()->evaluate(p,frame.snapshot->queryWorkspace())==expected,
              "Frame hit changed function value");
    }
    const double cold=median(cold_times),hit=median(hit_times);
    std::cout<<std::setprecision(17)<<"{\"raw_points\":"<<raw.size()
      <<",\"occupied_voxels\":"<<snapshot->occupancy().occupied_voxels.size()
      <<",\"execution_nodes\":"<<nodes<<",\"repeats\":"<<repeats
      <<",\"cold_prepare_ms\":"<<cold<<",\"cold_snapshot_build_ms\":"<<median(cold_build)
      <<",\"hit_prepare_ms\":"<<hit<<",\"hit_occupancy_ms\":"<<median(hit_occupancy)
      <<",\"reduction_fraction\":"<<(cold>0 ? 1-hit/cold : 0)
      <<",\"cache_hits\":"<<cache.hits()<<",\"cache_misses\":"<<cache.misses()
      <<",\"same_occupancy\":true,\"same_function\":true,\"same_snapshot\":true,\"workspaces_reused\":true}\n";
  }
  catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
