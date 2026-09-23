#include <rokae_demo/boundary_direct_tree.hpp>
#include <rokae_demo/link_distance_evaluator.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
using Clock=std::chrono::steady_clock;
using namespace rokae_demo;

void require(bool condition,const char* message)
{ if(!condition) throw std::runtime_error(message); }

double milliseconds(Clock::time_point start)
{ return std::chrono::duration<double,std::milli>(Clock::now()-start).count(); }

double median(std::vector<double> values)
{
  std::sort(values.begin(),values.end());
  const auto middle=values.size()/2;
  return values.size()%2 ? values[middle] : (values[middle-1]+values[middle])/2;
}

double p95(std::vector<double> values)
{
  std::sort(values.begin(),values.end());
  return values[static_cast<std::size_t>(std::ceil(.95*values.size()))-1];
}

std::vector<VoxelPoint> readCloud(const char* path)
{
  std::ifstream input(path); require(bool(input),"Cannot open frame benchmark cloud");
  std::vector<VoxelPoint> points; std::string line;
  while(std::getline(input,line))
  {
    const auto comment=line.find('#');
    if(comment!=std::string::npos) line.resize(comment);
    std::istringstream row(line); row>>std::ws; if(row.eof()) continue;
    double x,y,z;
    require(bool(row>>x>>y>>z) && std::isfinite(x) && std::isfinite(y) && std::isfinite(z),
            "Invalid frame benchmark point");
    row>>std::ws; require(row.eof(),"Frame benchmark cloud requires exactly three columns");
    points.emplace_back(x,y,z);
  }
  require(!input.bad() && !points.empty(),"Empty frame benchmark cloud");
  return points;
}

struct FrameTimes
{
  double occupancy=0,boundary=0,direct=0,field=0,workspace=0,query=0,total=0;
  double boundary_extract=0,boundary_merge=0;
  double direct_registry=0,direct_canonical=0,direct_partition=0;
  double direct_prefix=0,direct_recursion=0,direct_serialization=0;
  double program_compile=0,smooth_compile=0,evaluator_setup=0;
};
} // namespace

int main(int argc,char** argv)
{
  try
  {
    require(argc==5 || argc==6 || argc==7,
            "Usage: link_distance_frame_benchmark ROBOT.urdf CLOUD.xyz VOXEL_SIZE REPEATS [CSE:intern|legacy|none] [MATERIALIZE:true|false]");
    const double h=std::stod(argv[3]);
    const int repeats=std::stoi(argv[4]);
    const std::string cse=argc==5 ? "intern" : argv[5];
    require(cse=="intern" || cse=="legacy" || cse=="none" ||
            cse=="true" || cse=="false",
            "CSE must be intern, legacy, or none");
    const bool share=cse!="none" && cse!="false";
    const bool intern=share && cse!="legacy";
    const bool materialize=argc<7 ? false : std::string(argv[6])=="true";
    require(argc<7 || std::string(argv[6])=="true" || std::string(argv[6])=="false",
            "MATERIALIZE must be true or false");
    require(std::isfinite(h) && h>0 && repeats>0,"Invalid frame benchmark parameter");
    const auto robot=SphericalRobotModel::fromUrdf(argv[1]);
    require(robot.dof()==7,"Frame benchmark requires the supplied seven-DOF Panda");
    const auto points=readCloud(argv[2]);
    Eigen::VectorXd q(7); q<<0.0,-0.4,0.0,-2.0,0.0,1.6,0.8;

    BoundaryDirectOptions direct_options;
    direct_options.max_canonical_cells=262144;
    direct_options.max_expanded_nodes=32000000;
    direct_options.signs.max_strata=2000000;
    direct_options.signs.max_sign_operations=40000000000ULL;
    direct_options.intern_subexpressions=intern;
    direct_options.materialize_expression=materialize;
    direct_options.retain_partition_diagnostics=false;
    HS::SmoothTreeOptions smooth_options;
    smooth_options.max_error=.001;
    smooth_options.share_subexpressions=share;
    smooth_options.binary_kernel=true;
    smooth_options.predecode_kernel=true;
    smooth_options.paired_batch=true;
    LinkDistanceOptions link_options;
    link_options.max_smoothing_error=.001;
    link_options.global_margin=.005;
    link_options.field_workers=std::min(6u,std::max(1u,std::thread::hardware_concurrency()));
    link_options.dynamic_field_batch=link_options.field_workers>4;
    halfspace::BatchExecutor field_executor(link_options.field_workers);
    auto reusable_field_workspace=
        field_executor.makeWorkspace<halfspace::ValueGradient>(
            direct_options.max_nodes,true);

    std::vector<FrameTimes> trials; trials.reserve(repeats);
    SparseVoxelOccupancy last_occupancy;
    OctreeBoundary last_boundary;
    BoundaryDirectResult last_direct;
    std::vector<LinkDistanceResult> last_results;
    std::size_t execution_nodes=0;
    for(int trial=-1;trial<repeats;++trial)
    {
      FrameTimes time;
      const auto total_start=Clock::now();
      auto start=Clock::now();
      auto occupancy=buildSparseVoxelOccupancy(points,h);
      time.occupancy=milliseconds(start);
      start=Clock::now();
      auto direct=buildBoundaryDirectTree(occupancy,direct_options);
      time.direct=milliseconds(start);
      time.direct_registry=direct.registry_ms;
      time.direct_canonical=direct.canonical_ms;
      time.direct_partition=direct.tree_ms;
      time.direct_prefix=direct.prefix_ms;
      time.direct_recursion=direct.recursion_ms;
      time.direct_serialization=direct.serialization_ms;
      start=Clock::now();
      std::shared_ptr<const HS::CompiledTree> program;
      if(intern)
        program=std::make_shared<const HS::CompiledTree>(
            direct.serialized,direct.tree.leaves,HS::preinterned_program);
      else
        program=std::make_shared<const HS::CompiledTree>(
            direct.serialized,direct.tree.leaves,share);
      time.program_compile=milliseconds(start);
      start=Clock::now();
      HS::SmoothTree field(program,smooth_options);
      time.smooth_compile=milliseconds(start);
      start=Clock::now();
      LinkDistanceEvaluator evaluator(robot,field,link_options);
      time.evaluator_setup=milliseconds(start);
      time.field=time.program_compile+time.smooth_compile+time.evaluator_setup;
      start=Clock::now();
      auto workspace=evaluator.makeWorkspace(std::move(reusable_field_workspace));
      time.workspace=milliseconds(start);
      start=Clock::now();
      const auto& results=evaluator.evaluate(q,Eigen::Isometry3d::Identity(),workspace);
      time.query=milliseconds(start);
      reusable_field_workspace=std::move(*workspace.field_batch_workspace);
      time.total=milliseconds(total_start);
      require(results.size()==7,"Frame benchmark did not produce seven link results");
      if(trial>=0) trials.push_back(time);
      if(trial==repeats-1)
      {
        execution_nodes=field.nodeCount();
        last_results=results;
        last_occupancy=std::move(occupancy);
        last_direct=std::move(direct);
      }
    }

    // Certification is intentionally outside every measured frame.
    const auto validation_start=Clock::now();
    last_boundary=extractSparseVoxelBoundary(last_occupancy);
    validateSparseVoxelOccupancy(last_occupancy,points);
    const auto diagnostics=validateBoundaryDirectTree(
        last_occupancy,last_boundary,last_direct,direct_options);
    const double validation_ms=milliseconds(validation_start);
    require(diagnostics.applied && diagnostics.strict_sign_certified,
            "Frame benchmark result failed out-of-band certification");
    for(const auto& result:last_results)
      require(std::isfinite(result.distance_proxy) && result.distance_proxy>0,
              "Frame benchmark link result is not collision-free");

    const auto collect=[&](auto member) {
      std::vector<double> values; values.reserve(trials.size());
      for(const auto& trial:trials) values.push_back(trial.*member);
      return std::pair<double,double>{median(values),p95(values)};
    };
    const auto total=collect(&FrameTimes::total);
    std::cout<<std::setprecision(17)
      <<"{\"points\":"<<points.size()
      <<",\"voxel_size_m\":"<<h
      <<",\"occupied_voxels\":"<<last_occupancy.occupied_voxels.size()
      <<",\"boundary_patches\":"<<last_boundary.patches.size()
      <<",\"execution_nodes\":"<<execution_nodes
      <<",\"split_checks\":"<<last_direct.split_checks
      <<",\"expanded_nodes\":"<<last_direct.expanded_nodes
      <<",\"constructed_nodes\":"<<last_direct.constructed_nodes
      <<",\"repeats\":"<<repeats
      <<",\"workers\":"<<link_options.field_workers
      <<",\"structural_cse\":"<<(share ? "true" : "false")
      <<",\"construction_interning\":"<<(intern ? "true" : "false")
      <<",\"materialize_expression\":"<<(materialize ? "true" : "false")
      <<",\"frame_median_ms\":"<<total.first
      <<",\"frame_p95_ms\":"<<total.second;
    const auto emit=[&](const char* name,auto member) {
      const auto value=collect(member);
      std::cout<<",\""<<name<<"_median_ms\":"<<value.first
               <<",\""<<name<<"_p95_ms\":"<<value.second;
    };
    emit("occupancy",&FrameTimes::occupancy);
    emit("boundary",&FrameTimes::boundary);
    emit("boundary_extract",&FrameTimes::boundary_extract);
    emit("boundary_merge",&FrameTimes::boundary_merge);
    emit("direct_tree",&FrameTimes::direct);
    emit("direct_registry",&FrameTimes::direct_registry);
    emit("direct_canonical",&FrameTimes::direct_canonical);
    emit("direct_partition",&FrameTimes::direct_partition);
    emit("direct_prefix",&FrameTimes::direct_prefix);
    emit("direct_recursion",&FrameTimes::direct_recursion);
    emit("direct_serialization",&FrameTimes::direct_serialization);
    emit("field_compile",&FrameTimes::field);
    emit("program_compile",&FrameTimes::program_compile);
    emit("smooth_compile",&FrameTimes::smooth_compile);
    emit("evaluator_setup",&FrameTimes::evaluator_setup);
    emit("workspace",&FrameTimes::workspace);
    emit("query",&FrameTimes::query);
    std::cout<<",\"validation_outside_timing_ms\":"<<validation_ms
             <<",\"strict_sign_certified\":true}\n";
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr<<"link_distance_frame_benchmark failed: "<<error.what()<<'\n';
    return 1;
  }
}
