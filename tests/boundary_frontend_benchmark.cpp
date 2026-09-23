#include "boundary_frontend_checks.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace rokae_demo;
using namespace boundary_frontend_checks;
namespace ref = boundary_frontend_reference;
namespace
{
using Clock = std::chrono::steady_clock;
double ms(Clock::time_point start)
{ return std::chrono::duration<double,std::milli>(Clock::now()-start).count(); }
double median(std::vector<double> values)
{
  std::sort(values.begin(),values.end());
  return (values[(values.size()-1)/2]+values[values.size()/2])/2;
}
std::vector<VoxelPoint> readCloud(const char* path)
{
  std::ifstream in(path); require(bool(in),"Cannot read frontend benchmark XYZ");
  std::vector<VoxelPoint> points; std::string line;
  while(std::getline(in,line))
  {
    const auto comment = line.find('#'); if(comment != std::string::npos) line.resize(comment);
    std::istringstream row(line); row >> std::ws; if(row.eof()) continue;
    double x,y,z; require(bool(row>>x>>y>>z),"Invalid frontend benchmark XYZ"); row >> std::ws;
    require(row.eof() && std::isfinite(x) && std::isfinite(y) && std::isfinite(z),"Invalid frontend benchmark point");
    points.emplace_back(x,y,z);
  }
  require(!in.bad() && !points.empty(),"Empty frontend benchmark XYZ"); return points;
}
struct Sample
{
  OctreeBoundary boundary;
  BoundaryPlaneRegistry registry;
  double registry_ms = 0, frontend_ms = 0;
};
Sample run(const std::vector<VoxelIndex>& keys, bool legacy)
{
  const auto start = Clock::now(); Sample result;
  result.boundary = legacy ? ref::extractBoundary(keys,nullptr,keys.size()) : extractVoxelBoundary(keys);
  const auto registry_start = Clock::now();
  result.registry = legacy ? ref::makeRegistry(result.boundary) : buildBoundaryPlaneRegistry(result.boundary);
  result.registry_ms = ms(registry_start); result.frontend_ms = ms(start); return result;
}
} // namespace

int main(int argc, char** argv)
{
  try
  {
    require(argc==4,"Usage: boundary_frontend_benchmark INPUT.xyz VOXEL_SIZE REPEATS");
    std::size_t used = 0; const double h = std::stod(argv[2],&used);
    require(used==std::string(argv[2]).size() && std::isfinite(h) && h>0,"Invalid voxel size");
    const std::string count = argv[3];
    require(!count.empty() && count.find_first_not_of("0123456789")==std::string::npos,"Invalid repeat count");
    const auto repeats = std::stoull(count,&used);
    require(used==count.size() && repeats>0 && repeats<=10000,"Invalid repeat count");
    const auto raw = readCloud(argv[1]); const auto occupancy = buildSparseVoxelOccupancy(raw,h);
    const auto& keys = occupancy.occupied_voxels;
    // The final results stay alive during each pair. Checks and their
    // allocations are outside measured intervals; both modes get a warmup.
    std::array<std::array<std::vector<double>,4>,2> trials;
    std::size_t faces = 0, patches = 0, supports = 0;
    for(std::size_t trial = 0; trial < repeats+2; ++trial)
    {
      std::array<Sample,2> samples;
      for(unsigned j = 0; j < 2; ++j)
      {
        const unsigned mode = (trial+j)%2;
        samples[mode] = run(keys,mode==0);
      }
      sameBoundary(samples[0].boundary,samples[1].boundary);
      sameRegistry(samples[0].registry,samples[1].registry);
      validateVoxelBoundary(keys,samples[0].boundary); validateVoxelBoundary(keys,samples[1].boundary);
      for(unsigned mode = 0; mode < 2; ++mode) if(trial>=2)
      {
        const auto& s = samples[mode];
        const double values[]{s.boundary.extraction_ms,s.boundary.merge_ms,s.registry_ms,s.frontend_ms};
        for(unsigned metric = 0; metric < 4; ++metric) trials[mode][metric].push_back(values[metric]);
      }
      faces = samples[1].boundary.faces.size(); patches = samples[1].boundary.patches.size();
      supports = samples[1].registry.supports.size();
    }
    const char* names[]{"extraction_ms","merge_ms","registry_ms","frontend_ms"};
    const char* modes[]{"legacy","contiguous"};
    std::cout << std::setprecision(17) << "{\"raw_points\":" << raw.size() << ",\"occupied_voxels\":" << keys.size()
              << ",\"faces\":" << faces << ",\"patches\":" << patches << ",\"supports\":" << supports
              << ",\"repeats\":" << repeats << ",\"warmups\":2,\"exact_boundary\":true,\"exact_registry\":true";
    for(unsigned mode = 0; mode < 2; ++mode)
    {
      std::cout << ",\"" << modes[mode] << "\":{";
      for(unsigned metric = 0; metric < 4; ++metric)
      {
        if(metric) std::cout << ',';
        std::cout << '\"' << names[metric] << "\":" << median(trials[mode][metric]);
        std::cout << ",\"" << names[metric] << "_trials\":[";
        for(std::size_t i = 0; i < repeats; ++i) { if(i) std::cout << ','; std::cout << trials[mode][metric][i]; }
        std::cout << ']';
      }
      std::cout << '}';
    }
    std::cout << "}\n";
  }
  catch(const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
