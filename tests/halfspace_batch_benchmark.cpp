#include <rokae_demo/boundary_direct_frame_cache.hpp>
#include <algorithm>
#include <array>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>

using namespace rokae_demo;
namespace
{
using Clock=std::chrono::steady_clock;
void require(bool ok,const char* why) { if(!ok) throw std::runtime_error(why); }
bool same(double a,double b) { return std::memcmp(&a,&b,sizeof(double))==0; }
double elapsed(Clock::time_point start)
{ return std::chrono::duration<double,std::micro>(Clock::now()-start).count(); }
double median(std::vector<double> values)
{
  std::sort(values.begin(),values.end()); const auto n=values.size();
  return n%2 ? values[n/2] : (values[n/2-1]+values[n/2])/2;
}
struct Measurement
{
  std::array<double,3> medians;
  std::array<std::vector<double>,3> trials;
  double operator[](std::size_t mode) const { return medians[mode]; }
};
void writeTrials(const char* name,const Measurement& measured)
{
  std::cout<<",\""<<name<<"_trials_us\":[";
  for(std::size_t mode=0;mode<3;++mode)
  {
    std::cout<<(mode ? ",[" : "[");
    for(std::size_t i=0;i<measured.trials[mode].size();++i)
      std::cout<<(i ? "," : "")<<measured.trials[mode][i];
    std::cout<<']';
  }
  std::cout<<']';
}

std::vector<VoxelPoint> readCloud(const char* path)
{
  std::ifstream in(path); require(bool(in),"Cannot read batch benchmark XYZ");
  std::vector<VoxelPoint> points; std::string line;
  while(std::getline(in,line))
  {
    const auto comment=line.find('#'); if(comment!=std::string::npos) line.resize(comment);
    std::istringstream row(line); row>>std::ws; if(row.eof()) continue;
    double x,y,z; require(bool(row>>x>>y>>z),"Invalid batch benchmark XYZ"); row>>std::ws;
    require(row.eof() && std::isfinite(x) && std::isfinite(y) && std::isfinite(z),"Invalid batch benchmark point");
    points.emplace_back(x,y,z);
  }
  require(!in.bad() && !points.empty(),"Empty batch benchmark XYZ"); return points;
}

template<class First,class Second,class Third,class Verify>
Measurement measureModes(std::size_t repeats,std::size_t point_count,First&& first,Second&& second,Third&& third,Verify&& verify)
{
  std::array<std::vector<double>,3> times; for(auto& values:times) values.reserve(repeats);
  const auto run=[&](std::size_t mode) { if(mode==0) first(); else if(mode==1) second(); else third(); };
  for(std::size_t mode=0;mode<3;++mode) run(mode);
  for(std::size_t trial=0;trial<repeats;++trial)
  {
    const auto offset=trial%3;
    for(std::size_t step=0;step<3;++step)
    {
      const auto mode=(offset+step)%3; const auto start=Clock::now(); run(mode);
      times[mode].push_back(elapsed(start)/point_count);
    }
    verify(); // Consume and compare all timed outputs outside the measured interval.
  }
  return {{median(times[0]),median(times[1]),median(times[2])},std::move(times)};
}
} // namespace

int main(int argc,char** argv)
{
  try
  {
    require(argc==6,"Usage: halfspace_batch_benchmark INPUT.xyz VOXEL_SIZE POINTS REPEATS WORKERS");
    std::size_t used=0; const double h=std::stod(argv[2],&used);
    require(used==std::string(argv[2]).size() && std::isfinite(h) && h>0,"Invalid voxel size");
    const auto parseSize=[&](int index) {
      const std::string argument=argv[index];
      require(!argument.empty() && argument.find_first_not_of("0123456789")==std::string::npos,"Invalid positive integer");
      std::size_t parsed=0; const auto value=std::stoull(argv[index],&parsed);
      require(parsed==argument.size() && value>0 && value<=std::numeric_limits<std::size_t>::max(),"Invalid positive integer");
      return static_cast<std::size_t>(value);
    };
    const auto point_count=parseSize(3),repeats=parseSize(4),workers=parseSize(5);
    const auto raw=readCloud(argv[1]); BoundaryDirectFrameOptions options;
    options.direct.max_canonical_cells=131072;
    options.direct.max_expanded_nodes=16000000;
    options.direct.signs.max_strata=1000000;
    options.direct.signs.max_sign_operations=20000000000ULL;
    options.smooth=true; options.smoothing={0,.001,true,true};
    BoundaryDirectFrameCache cache; auto prepared=cache.prepare(raw,h,options); auto& snapshot=*prepared.snapshot;
    require(snapshot.diagnostics().strict_sign_certified && snapshot.smoothTree(),"Batch benchmark snapshot not certified/smoothed");
    const auto& program=*snapshot.compiledProgram(); const auto& smooth=*snapshot.smoothTree();

    std::mt19937_64 random(20260909); std::uniform_real_distribution<double> unit(-.05,1.05);
    const auto& box=snapshot.occupancy().sampling_box; std::vector<HS::Vec3> points; points.reserve(point_count);
    for(std::size_t i=0;i<point_count;++i)
      points.push_back({box.lower[0]+unit(random)*(box.upper[0]-box.lower[0]),
        box.lower[1]+unit(random)*(box.upper[1]-box.lower[1]),box.lower[2]+unit(random)*(box.upper[2]-box.lower[2])});

    std::vector<double> hard_expected(point_count),hard_serial(point_count),hard_parallel(point_count);
    std::vector<double> smooth_expected(point_count),smooth_serial(point_count),smooth_parallel(point_count);
    std::vector<HS::ValueGradient> gradient_expected(point_count),gradient_serial(point_count),gradient_parallel(point_count);
    auto hard_work=program.makeWorkspace(),smooth_work=smooth.makeValueWorkspace(); auto gradient_work=smooth.makeGradientWorkspace();
    for(std::size_t i=0;i<point_count;++i)
    {
      hard_expected[i]=program.evaluate(points[i],hard_work);
      smooth_expected[i]=smooth.evaluate(points[i],smooth_work);
      gradient_expected[i]=smooth.evaluateWithGradient(points[i],gradient_work);
    }
    auto hard_one=program.makeBatchWorkspace(1),hard_many=program.makeBatchWorkspace(workers);
    auto smooth_one=smooth.makeValueBatchWorkspace(1),smooth_many=smooth.makeValueBatchWorkspace(workers);
    auto gradient_one=smooth.makeGradientBatchWorkspace(1),gradient_many=smooth.makeGradientBatchWorkspace(workers);
    program.evaluateBatch(points,hard_serial,hard_one); program.evaluateBatch(points,hard_parallel,hard_many);
    smooth.evaluateBatch(points,smooth_serial,smooth_one); smooth.evaluateBatch(points,smooth_parallel,smooth_many);
    smooth.evaluateWithGradientBatch(points,gradient_serial,gradient_one);
    smooth.evaluateWithGradientBatch(points,gradient_parallel,gradient_many);
    const auto verify=[&] {
      for(std::size_t i=0;i<point_count;++i)
      {
        require(same(hard_expected[i],hard_serial[i]) && same(hard_expected[i],hard_parallel[i]),"Hard batch result changed");
        require(same(smooth_expected[i],smooth_serial[i]) && same(smooth_expected[i],smooth_parallel[i]),"Smooth batch result changed");
        const auto compare=[&](const auto& value) { return same(value.value,gradient_expected[i].value) &&
          same(value.gradient.x,gradient_expected[i].gradient.x) && same(value.gradient.y,gradient_expected[i].gradient.y) &&
          same(value.gradient.z,gradient_expected[i].gradient.z); };
        require(compare(gradient_serial[i]) && compare(gradient_parallel[i]),"Gradient batch result changed");
      }
    };
    verify();

    const auto hard_times=measureModes(repeats,point_count,[&] {
      for(std::size_t i=0;i<point_count;++i) hard_expected[i]=program.evaluate(points[i],hard_work);
    },[&] { program.evaluateBatch(points,hard_serial,hard_one); },
      [&] { program.evaluateBatch(points,hard_parallel,hard_many); },verify);
    const auto smooth_times=measureModes(repeats,point_count,[&] {
      for(std::size_t i=0;i<point_count;++i) smooth_expected[i]=smooth.evaluate(points[i],smooth_work);
    },[&] { smooth.evaluateBatch(points,smooth_serial,smooth_one); },
      [&] { smooth.evaluateBatch(points,smooth_parallel,smooth_many); },verify);
    const auto gradient_times=measureModes(repeats,point_count,[&] {
      for(std::size_t i=0;i<point_count;++i) gradient_expected[i]=smooth.evaluateWithGradient(points[i],gradient_work);
    },[&] {
      smooth.evaluateWithGradientBatch(points,gradient_serial,gradient_one);
    },[&] {
      smooth.evaluateWithGradientBatch(points,gradient_parallel,gradient_many);
    },verify);
    const auto hard_scalar_us=hard_times[0],hard_batch_one_us=hard_times[1],hard_parallel_us=hard_times[2];
    const auto smooth_scalar_us=smooth_times[0],smooth_batch_one_us=smooth_times[1],smooth_parallel_us=smooth_times[2];
    const auto gradient_scalar_us=gradient_times[0],gradient_batch_one_us=gradient_times[1],gradient_parallel_us=gradient_times[2];

    std::cout<<std::setprecision(17)<<"{\"points\":"<<point_count<<",\"workers\":"<<workers
      <<",\"hardware_threads\":"<<std::thread::hardware_concurrency()<<",\"execution_nodes\":"<<program.nodeCount()
      <<",\"repeats\":"<<repeats<<",\"hard_scalar_us\":"<<hard_scalar_us
      <<",\"hard_batch_one_us\":"<<hard_batch_one_us<<",\"hard_parallel_us\":"<<hard_parallel_us
      <<",\"hard_speedup\":"<<hard_scalar_us/hard_parallel_us
      <<",\"smooth_scalar_us\":"<<smooth_scalar_us<<",\"smooth_batch_one_us\":"<<smooth_batch_one_us
      <<",\"smooth_parallel_us\":"<<smooth_parallel_us<<",\"smooth_speedup\":"<<smooth_scalar_us/smooth_parallel_us
      <<",\"gradient_scalar_us\":"<<gradient_scalar_us<<",\"gradient_batch_one_us\":"<<gradient_batch_one_us
      <<",\"gradient_parallel_us\":"<<gradient_parallel_us<<",\"gradient_speedup\":"<<gradient_scalar_us/gradient_parallel_us
      <<",\"hard_batch_workspace_bytes\":"<<hard_many.storageBytes()
      <<",\"smooth_batch_workspace_bytes\":"<<smooth_many.storageBytes()
      <<",\"gradient_batch_workspace_bytes\":"<<gradient_many.storageBytes()
      <<",\"exact_hard\":true,\"exact_smooth\":true,\"exact_gradient\":true";
    writeTrials("hard",hard_times); writeTrials("smooth",smooth_times); writeTrials("gradient",gradient_times);
    std::cout<<"}\n";
  }
  catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
