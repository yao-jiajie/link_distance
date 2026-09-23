#include <rokae_demo/boundary_direct_frame_cache.hpp>
#include <rokae_demo/link_distance_evaluator.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace
{
using Clock = std::chrono::steady_clock;
using namespace rokae_demo;

void require(bool ok, const char* why)
{
  if(!ok) throw std::runtime_error(why);
}

bool same(double a, double b)
{ return std::memcmp(&a, &b, sizeof(double)) == 0; }

std::vector<VoxelPoint> readCloud(const char* path)
{
  std::ifstream input(path);
  require(bool(input), "Cannot open benchmark cloud");
  std::vector<VoxelPoint> points;
  std::string line;
  while(std::getline(input, line))
  {
    const auto comment = line.find('#');
    if(comment != std::string::npos) line.resize(comment);
    std::istringstream row(line);
    row >> std::ws;
    if(row.eof()) continue;
    double x, y, z;
    require(bool(row >> x >> y >> z) && std::isfinite(x) &&
                std::isfinite(y) && std::isfinite(z),
            "Invalid benchmark cloud point");
    row >> std::ws;
    require(row.eof(), "Benchmark cloud requires exactly three columns");
    points.emplace_back(x, y, z);
  }
  require(!input.bad() && !points.empty(), "Empty benchmark cloud");
  return points;
}

double median(std::vector<double> values)
{
  std::sort(values.begin(), values.end());
  const auto middle = values.size() / 2;
  return values.size() % 2 ? values[middle]
    : (values[middle - 1] + values[middle]) / 2.0;
}

double p95(std::vector<double> values)
{
  std::sort(values.begin(), values.end());
  return values[static_cast<std::size_t>(std::ceil(0.95 * values.size())) - 1];
}

void compare(const LinkDistanceWorkspace& a, const LinkDistanceWorkspace& b)
{
  require(a.field_samples.size() == b.field_samples.size() &&
              a.results.size() == b.results.size(), "Kernel result size changed");
  for(std::size_t i = 0; i < a.field_samples.size(); ++i)
  {
    const auto& x = a.field_samples[i];
    const auto& y = b.field_samples[i];
    require(same(x.value, y.value) && same(x.gradient.x, y.gradient.x) &&
                same(x.gradient.y, y.gradient.y) &&
                same(x.gradient.z, y.gradient.z) &&
                same(a.sphere_clearances[i], b.sphere_clearances[i]),
            "Kernel changed sphere field value, gradient, or clearance bits");
  }
  for(std::size_t i = 0; i < a.results.size(); ++i)
  {
    const auto& x = a.results[i];
    const auto& y = b.results[i];
    require(x.nearest_sphere == y.nearest_sphere &&
                same(x.distance_proxy, y.distance_proxy) &&
                same(x.hard_min_proxy, y.hard_min_proxy) &&
                x.gradient_q.size() == y.gradient_q.size(),
            "Kernel changed link result bits");
    for(Eigen::Index joint = 0; joint < x.gradient_q.size(); ++joint)
      require(same(x.gradient_q[joint], y.gradient_q[joint]),
              "Kernel changed link gradient bits");
  }
}

} // namespace

int main(int argc, char** argv)
{
  try
  {
    require(argc >= 6 && argc <= 8,
            "Usage: link_distance_kernel_benchmark ROBOT.urdf CLOUD.xyz VOXEL_SIZE FIELD_Y REPEATS [WORKERS [DYNAMIC_WORKERS]]");
    const double h = std::stod(argv[3]);
    const double field_y = std::stod(argv[4]);
    const int repeats = std::stoi(argv[5]);
    const int workers = argc >= 7 ? std::stoi(argv[6]) : 4;
    const int dynamic_workers = argc == 8 ? std::stoi(argv[7]) : workers;
    require(std::isfinite(h) && h > 0.0 && std::isfinite(field_y) && repeats > 0 &&
                workers > 0 && dynamic_workers > 0,
            "Invalid benchmark parameter");

    const auto robot = SphericalRobotModel::fromUrdf(argv[1]);
    require(robot.dof() == 7 && robot.spheres().size() == 59,
            "Benchmark requires the supplied Panda sphere model");
    const auto points = readCloud(argv[2]);
    BoundaryDirectFrameOptions options;
    // Also admits the six-shape h=0.02 stress cloud. These are rejection
    // limits only and do not add work to smaller benchmark inputs.
    options.direct.max_canonical_cells = 262144;
    options.direct.max_expanded_nodes = 32000000;
    options.direct.signs.max_strata = 2000000;
    options.direct.signs.max_sign_operations = 40000000000ULL;
    options.direct.intern_subexpressions = true;
    options.smooth = true;
    options.share_subexpressions = true;
    options.smoothing = {0.0, 0.001, true, true, false};
    BoundaryDirectFrameCache cache(1);
    const auto prepared = cache.prepare(points, h, options);
    require(prepared.snapshot && prepared.snapshot->diagnostics().strict_sign_certified,
            "Benchmark environment is not strictly certified");
    const auto& program = prepared.snapshot->compiledProgram();
    halfspace::SmoothTree legacy(program, {0.0, 0.001, true, true, false});
    halfspace::SmoothTree decoded(program, {0.0, 0.001, true, true, true});
    halfspace::SmoothTree paired(program, {0.0, 0.001, true, true, true, true});
    LinkDistanceOptions link_options;
    link_options.max_smoothing_error = 0.001;
    link_options.global_margin = 0.005;
    link_options.field_workers = static_cast<std::size_t>(workers);
    LinkDistanceEvaluator old_evaluator(robot, legacy, link_options);
    LinkDistanceEvaluator new_evaluator(robot, decoded, link_options);
    link_options.dynamic_field_batch = true;
    link_options.field_workers = static_cast<std::size_t>(dynamic_workers);
    LinkDistanceEvaluator dynamic_evaluator(robot, decoded, link_options);
    LinkDistanceEvaluator paired_evaluator(robot, paired, link_options);
    auto old_workspace = old_evaluator.makeWorkspace();
    auto new_workspace = new_evaluator.makeWorkspace();
    auto dynamic_workspace = dynamic_evaluator.makeWorkspace();
    auto paired_workspace = paired_evaluator.makeWorkspace();
    Eigen::VectorXd q(7);
    q << 0.0, -0.4, 0.0, -2.0, 0.0, 1.6, 0.8;
    Eigen::Isometry3d field_T_base = Eigen::Isometry3d::Identity();
    field_T_base.translation().y() = field_y;
    const auto run = [&](int mode) {
      if(mode == 0) old_evaluator.evaluate(q, field_T_base, old_workspace);
      else if(mode == 1) new_evaluator.evaluate(q, field_T_base, new_workspace);
      else if(mode == 2) dynamic_evaluator.evaluate(q, field_T_base, dynamic_workspace);
      else paired_evaluator.evaluate(q, field_T_base, paired_workspace);
    };
    for(int warmup = 0; warmup < 6; ++warmup)
      for(int mode = 0; mode < 4; ++mode) run((warmup + mode) % 4);
    compare(old_workspace, new_workspace);
    compare(new_workspace, dynamic_workspace);
    compare(dynamic_workspace, paired_workspace);
    std::vector<double> timings[4];
    for(auto& timing : timings) timing.reserve(repeats);
    for(int trial = 0; trial < repeats; ++trial)
    {
      for(int step = 0; step < 4; ++step)
      {
        const int mode = (trial + step) % 4;
        const auto start = Clock::now();
        run(mode);
        timings[mode].push_back(
            std::chrono::duration<double, std::milli>(Clock::now() - start).count());
      }
      compare(old_workspace, new_workspace); // outside all timed intervals
      compare(new_workspace, dynamic_workspace);
      compare(dynamic_workspace, paired_workspace);
    }
    // Production-only repetition avoids cache interference from rotating the
    // legacy comparison fields. Keep it separate from the paired A/B results.
    std::vector<double> isolated_paired;
    isolated_paired.reserve(repeats);
    for(int warmup=0;warmup<6;++warmup) run(3);
    for(int trial=0;trial<repeats;++trial)
    {
      const auto start=Clock::now();
      run(3);
      isolated_paired.push_back(
          std::chrono::duration<double,std::milli>(Clock::now()-start).count());
    }
    std::cout << std::setprecision(17)
              << "{\"voxel_size_m\":" << h
              << ",\"points\":" << points.size()
              << ",\"spheres\":" << robot.spheres().size()
              << ",\"execution_nodes\":" << decoded.nodeCount()
              << ",\"repeats\":" << repeats
              << ",\"workers\":" << workers
              << ",\"dynamic_workers\":" << dynamic_workers
              << ",\"legacy_median_ms\":" << median(timings[0])
              << ",\"legacy_p95_ms\":" << p95(timings[0])
              << ",\"predecoded_median_ms\":" << median(timings[1])
              << ",\"predecoded_p95_ms\":" << p95(timings[1])
              << ",\"dynamic_median_ms\":" << median(timings[2])
              << ",\"dynamic_p95_ms\":" << p95(timings[2])
              << ",\"paired_median_ms\":" << median(timings[3])
              << ",\"paired_p95_ms\":" << p95(timings[3])
              << ",\"isolated_paired_median_ms\":" << median(isolated_paired)
              << ",\"isolated_paired_p95_ms\":" << p95(isolated_paired)
              << ",\"exact_outputs\":true}\n";
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr << "link_distance_kernel_benchmark failed: " << error.what() << '\n';
    return 1;
  }
}
