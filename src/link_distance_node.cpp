#include <rokae_demo/boundary_direct_frame_cache.hpp>
#include <rokae_demo/link_distance_evaluator.hpp>
#include <rokae_demo/spherical_robot_model.hpp>

#include <Eigen/Geometry>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace
{
namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;

struct Options
{
  fs::path urdf;
  fs::path cloud;
  fs::path output;
  double voxel_size = 0.0;
  double environment_beta = 0.0;
  double environment_error = 0.0;
  double link_alpha = 0.0;
  double link_error = 0.0;
  double margin = 0.0;
  std::size_t query_workers = 1;
  bool profile_query = false;
  bool stream = false;
  bool reuse_environment = true;
  std::vector<double> configuration;
  Eigen::Isometry3d field_T_base = Eigen::Isometry3d::Identity();
  rokae_demo::BoundaryDirectOptions direct_limits;
};

double milliseconds(Clock::time_point start)
{
  return std::chrono::duration<double, std::milli>(Clock::now() - start)
      .count();
}

double finiteDouble(const std::string& text, const char* description)
{
  std::size_t parsed = 0;
  const double result = std::stod(text, &parsed);
  if(parsed != text.size() || !std::isfinite(result))
    throw std::invalid_argument(std::string("Invalid ") + description +
                                ": " + text);
  return result;
}

std::size_t positiveSize(const std::string& text, const char* description)
{
  std::size_t parsed = 0;
  unsigned long long value = 0;
  try
  {
    value = std::stoull(text, &parsed);
  }
  catch(const std::exception&)
  {
    throw std::invalid_argument(std::string("Invalid ") + description +
                                ": " + text);
  }
  if(parsed != text.size() || value == 0 ||
     value > std::numeric_limits<std::size_t>::max())
    throw std::invalid_argument(std::string("Invalid ") + description +
                                ": " + text);
  return static_cast<std::size_t>(value);
}

std::vector<double> commaList(const std::string& text,
                              const char* description)
{
  std::vector<double> result;
  std::istringstream input(text);
  std::string item;
  while(std::getline(input, item, ','))
  {
    if(item.empty())
      throw std::invalid_argument(std::string("Empty value in ") +
                                  description);
    result.push_back(finiteDouble(item, description));
  }
  if(result.empty() || (!text.empty() && text.back() == ','))
    throw std::invalid_argument(std::string("Invalid ") + description);
  return result;
}

Eigen::Isometry3d transform(const std::vector<double>& values)
{
  if(values.size() != 6)
    throw std::invalid_argument(
        "--field-from-base requires x,y,z,roll,pitch,yaw");
  Eigen::Isometry3d result = Eigen::Isometry3d::Identity();
  result.translation() = Eigen::Vector3d(values[0], values[1], values[2]);
  result.linear() =
      (Eigen::AngleAxisd(values[5], Eigen::Vector3d::UnitZ()) *
       Eigen::AngleAxisd(values[4], Eigen::Vector3d::UnitY()) *
       Eigen::AngleAxisd(values[3], Eigen::Vector3d::UnitX()))
          .toRotationMatrix();
  return result;
}

void usage(const char* program)
{
  std::cerr
      << "Usage: " << program
      << " robot.urdf environment.xyz output.json [options]\n"
      << "Required:\n"
      << "  --voxel-size METERS\n"
      << "  --q q1,q2,...,qn\n"
      << "  exactly one of --env-lse-beta 1/M or --env-lse-error M\n"
      << "  exactly one of --link-lse-alpha 1/M or --link-lse-error M\n"
      << "Optional:\n"
      << "  --margin METERS\n"
      << "  --field-from-base x,y,z,roll,pitch,yaw (default identity)\n"
      << "  --query-workers COUNT (default min(6, hardware concurrency))\n"
      << "  --profile-query true|false (optional phase times in JSON)\n"
      << "  --stream true|false (with cloud/output set to -; stdin rows: cloud<TAB>output[<TAB>q_csv])\n"
      << "  --reuse-environment true|false (stream A/B switch; default true)\n"
      << "Direct-tree work limits:\n"
      << "  --boundary-max-cells COUNT\n"
      << "  --boundary-max-direct-nodes COUNT\n"
      << "  --boundary-max-expanded-nodes COUNT\n"
      << "  --boundary-max-split-checks COUNT\n"
      << "  --boundary-max-direct-depth COUNT (<=256)\n"
      << "  --boundary-max-strata COUNT\n"
      << "  --boundary-max-sign-ops COUNT\n";
}

Options parse(int argc, char** argv)
{
  if(argc < 4)
    throw std::invalid_argument("Expected URDF, XYZ cloud, and output JSON");
  Options options;
  // The generic library defaults intentionally target small scenes. These
  // bounded node defaults also cover the supplied camera cloud at h=0.01 m.
  options.direct_limits.max_canonical_cells = 131072;
  options.direct_limits.max_expanded_nodes = 16000000;
  options.direct_limits.signs.max_strata = 1000000;
  options.direct_limits.signs.max_sign_operations = 20000000000ULL;
  options.query_workers = std::min(
      6u, std::max(1u, std::thread::hardware_concurrency()));
  options.urdf = argv[1];
  options.cloud = argv[2];
  options.output = argv[3];
  std::set<std::string> seen;
  for(int index = 4; index < argc; index += 2)
  {
    if(index + 1 >= argc)
      throw std::invalid_argument("Missing option value");
    const std::string flag = argv[index];
    const std::string value = argv[index + 1];
    if(!seen.insert(flag).second)
      throw std::invalid_argument("Repeated option: " + flag);
    if(flag == "--voxel-size")
      options.voxel_size = finiteDouble(value, "voxel size");
    else if(flag == "--q")
      options.configuration = commaList(value, "configuration");
    else if(flag == "--env-lse-beta")
      options.environment_beta = finiteDouble(value, "environment beta");
    else if(flag == "--env-lse-error")
      options.environment_error = finiteDouble(value, "environment error");
    else if(flag == "--link-lse-alpha")
      options.link_alpha = finiteDouble(value, "link alpha");
    else if(flag == "--link-lse-error")
      options.link_error = finiteDouble(value, "link error");
    else if(flag == "--margin")
      options.margin = finiteDouble(value, "margin");
    else if(flag == "--field-from-base")
      options.field_T_base = transform(commaList(value, "base transform"));
    else if(flag == "--query-workers")
      options.query_workers = positiveSize(value, "query workers");
    else if(flag == "--profile-query")
    {
      if(value != "true" && value != "false")
        throw std::invalid_argument("--profile-query requires true or false");
      options.profile_query = value == "true";
    }
    else if(flag == "--stream" || flag == "--reuse-environment")
    {
      if(value != "true" && value != "false")
        throw std::invalid_argument(flag + " requires true or false");
      if(flag == "--stream") options.stream = value == "true";
      else options.reuse_environment = value == "true";
    }
    else if(flag == "--boundary-max-cells")
      options.direct_limits.max_canonical_cells =
          positiveSize(value, "boundary max cells");
    else if(flag == "--boundary-max-direct-nodes")
      options.direct_limits.max_nodes =
          positiveSize(value, "boundary max direct nodes");
    else if(flag == "--boundary-max-expanded-nodes")
      options.direct_limits.max_expanded_nodes =
          positiveSize(value, "boundary max expanded nodes");
    else if(flag == "--boundary-max-split-checks")
      options.direct_limits.max_split_checks =
          positiveSize(value, "boundary max split checks");
    else if(flag == "--boundary-max-direct-depth")
    {
      const std::size_t depth =
          positiveSize(value, "boundary max direct depth");
      if(depth > 256)
        throw std::invalid_argument(
            "Boundary max direct depth must not exceed 256");
      options.direct_limits.max_depth = static_cast<unsigned>(depth);
    }
    else if(flag == "--boundary-max-strata")
      options.direct_limits.signs.max_strata =
          positiveSize(value, "boundary max strata");
    else if(flag == "--boundary-max-sign-ops")
      options.direct_limits.signs.max_sign_operations =
          positiveSize(value, "boundary max sign operations");
    else
      throw std::invalid_argument("Unknown option: " + flag);
  }
  if(!(options.voxel_size > 0.0) || options.configuration.empty() ||
     options.margin < 0.0 ||
     ((options.environment_beta > 0.0) ==
      (options.environment_error > 0.0)) ||
     ((options.link_alpha > 0.0) == (options.link_error > 0.0)))
    throw std::invalid_argument(
        "Missing/invalid voxel, q, LSE, or margin option");
  if(options.stream && (options.cloud != "-" || options.output != "-"))
    throw std::invalid_argument("--stream true requires cloud and output positional arguments to be -");
  if(!options.stream && !options.reuse_environment)
    throw std::invalid_argument("--reuse-environment false requires --stream true");
  return options;
}

Eigen::VectorXd jointConfiguration(const std::vector<double>& values,
                                   std::size_t dof)
{
  if(values.size() != dof)
    throw std::invalid_argument("--q has " + std::to_string(values.size()) +
                                " values, but URDF requires " +
                                std::to_string(dof));
  Eigen::VectorXd result(static_cast<Eigen::Index>(dof));
  for(std::size_t index = 0; index < dof; ++index)
    result[static_cast<Eigen::Index>(index)] = values[index];
  return result;
}

std::vector<rokae_demo::VoxelPoint> readCloud(const fs::path& path)
{
  std::ifstream input(path);
  if(!input)
    throw std::runtime_error("Cannot read XYZ point cloud: " + path.string());
  std::vector<rokae_demo::VoxelPoint> points;
  std::string line;
  while(std::getline(input, line))
  {
    const std::size_t comment = line.find('#');
    if(comment != std::string::npos)
      line.resize(comment);
    std::istringstream row(line);
    row >> std::ws;
    if(row.eof())
      continue;
    double x = 0.0, y = 0.0, z = 0.0;
    if(!(row >> x >> y >> z) || !std::isfinite(x) || !std::isfinite(y) ||
       !std::isfinite(z))
      throw std::runtime_error("Invalid XYZ point in " + path.string());
    row >> std::ws;
    if(!row.eof())
      throw std::runtime_error("XYZ input requires exactly three columns");
    points.emplace_back(x, y, z);
  }
  if(input.bad() || points.empty())
    throw std::runtime_error("Empty or unreadable XYZ point cloud");
  return points;
}

std::string jsonString(const std::string& value)
{
  std::ostringstream output;
  output << '"';
  for(const unsigned char c : value)
  {
    switch(c)
    {
      case '"': output << "\\\""; break;
      case '\\': output << "\\\\"; break;
      case '\b': output << "\\b"; break;
      case '\f': output << "\\f"; break;
      case '\n': output << "\\n"; break;
      case '\r': output << "\\r"; break;
      case '\t': output << "\\t"; break;
      default:
        if(c < 0x20)
          output << "\\u" << std::hex << std::setw(4)
                 << std::setfill('0') << static_cast<unsigned>(c)
                 << std::dec << std::setfill(' ');
        else
          output << static_cast<char>(c);
    }
  }
  output << '"';
  return output.str();
}

template<class Derived>
void jsonVector(std::ostream& output,
                const Eigen::MatrixBase<Derived>& values)
{
  output << '[';
  for(Eigen::Index index = 0; index < values.size(); ++index)
    output << (index ? "," : "") << values[index];
  output << ']';
}

void writeOutput(const Options& options,
                 const rokae_demo::SphericalRobotModel& robot,
                 const rokae_demo::SparseVoxelOccupancy& occupancy,
                 const rokae_demo::OctreeBoundary& boundary,
                 const rokae_demo::halfspace::SmoothTree& field,
                 const rokae_demo::LinkDistanceEvaluator& evaluator,
                 const std::vector<rokae_demo::LinkDistanceResult>& results,
                 const rokae_demo::LinkDistanceWorkspace& workspace,
                 const Eigen::VectorXd& configuration,
                 std::size_t point_count,
                 double environment_build_ms,
                 double query_ms,
                 const rokae_demo::LinkDistancePhaseTimes* phase_times)
{
  if(!options.output.parent_path().empty())
    fs::create_directories(options.output.parent_path());
  std::ofstream output(options.output);
  if(!output)
    throw std::runtime_error("Cannot write output JSON: " +
                             options.output.string());
  output << std::setprecision(17);
  output << "{\n  \"schema\":1,\n"
         << "  \"value_kind\":\"smooth_clearance_proxy_not_euclidean_sdf\",\n"
         << "  \"length_unit\":\"m\",\n"
         << "  \"gradient_order\":\"joint_order\",\n"
         << "  \"robot\":" << jsonString(robot.name()) << ",\n"
         << "  \"urdf\":" << jsonString(options.urdf.string()) << ",\n"
         << "  \"cloud\":" << jsonString(options.cloud.string()) << ",\n"
         << "  \"point_count\":" << point_count << ",\n"
         << "  \"voxel_size_m\":" << options.voxel_size << ",\n"
         << "  \"environment_beta_per_m\":" << field.beta() << ",\n"
         << "  \"environment_scalar_smoothing_bound_m\":"
         << field.errorBound() << ",\n"
         << "  \"link_alpha_per_m\":" << evaluator.alpha() << ",\n"
         << "  \"maximum_link_smoothing_bound_m\":"
         << evaluator.maximumLinkSmoothingBound() << ",\n"
         << "  \"global_margin_m\":" << options.margin << ",\n"
         << "  \"environment_direct_limits\":{"
         << "\"max_canonical_cells\":"
         << options.direct_limits.max_canonical_cells
         << ",\"max_direct_nodes\":" << options.direct_limits.max_nodes
         << ",\"max_expanded_nodes\":"
         << options.direct_limits.max_expanded_nodes
         << ",\"max_split_checks\":"
         << options.direct_limits.max_split_checks
         << ",\"max_direct_depth\":" << options.direct_limits.max_depth
         << ",\"max_strata\":"
         << options.direct_limits.signs.max_strata
         << ",\"max_sign_operations\":"
         << options.direct_limits.signs.max_sign_operations << "},\n"
         << "  \"environment_build_ms\":" << environment_build_ms << ",\n"
         << "  \"query_workers\":" << evaluator.fieldWorkerCount() << ",\n"
         << "  \"query_ms\":" << query_ms << ",\n";
  if(phase_times)
    output << "  \"query_phase_ms\":{\"kinematics\":"
           << phase_times->kinematics_ms << ",\"field\":"
           << phase_times->field_ms << ",\"link\":"
           << phase_times->link_ms << "},\n";
  output << "  \"joint_order\":[";
  for(std::size_t index = 0; index < robot.jointNames().size(); ++index)
    output << (index ? "," : "") << jsonString(robot.jointNames()[index]);
  output << "],\n  \"configuration\":";
  jsonVector(output, configuration);
  output << ",\n  \"field_from_base\":[";
  for(Eigen::Index row = 0; row < 4; ++row)
    for(Eigen::Index column = 0; column < 4; ++column)
      output << ((row || column) ? "," : "")
             << options.field_T_base.matrix()(row, column);
  output << "],\n  \"voxel_occupancy\":{\n"
         << "    \"occupied_voxels\":" << occupancy.occupied_voxels.size()
         << ",\n    \"faces_before_common_face_removal\":"
         << boundary.voxel_faces_before
         << ",\n    \"internal_common_faces_removed\":"
         << boundary.internal_faces_removed / 2
         << ",\n    \"internal_directed_faces_removed\":"
         << boundary.internal_faces_removed
         << ",\n    \"exposed_unit_faces\":" << boundary.faces.size()
         << ",\n    \"merged_boundary_patches\":"
         << boundary.patches.size()
         << ",\n    \"coplanar_exposed_faces_merged\":"
         << boundary.faces.size() - boundary.patches.size()
         << ",\n    \"boundary_patches\":[\n";
  for(std::size_t index = 0; index < boundary.patches.size(); ++index)
  {
    const auto& patch = boundary.patches[index];
    const auto corners = rokae_demo::octreeBoundaryCorners(patch);
    output << (index ? ",\n" : "") << "      {\"axis\":" << patch.axis
           << ",\"sign\":" << patch.sign << ",\"corners_field_m\":[";
    for(std::size_t corner = 0; corner < corners.size(); ++corner)
    {
      output << (corner ? "," : "") << '[';
      for(std::size_t axis = 0; axis < 3; ++axis)
        output << (axis ? "," : "")
               << static_cast<double>(corners[corner][axis]) *
                      occupancy.voxel_size;
      output << ']';
    }
    output << "]}";
  }
  output << "\n    ]\n  },\n  \"links\":[\n";
  for(std::size_t index = 0; index < results.size(); ++index)
  {
    const auto& result = results[index];
    output << (index ? ",\n" : "") << "    {\"name\":"
           << jsonString(result.link_name)
           << ",\"joint\":" << jsonString(result.joint_name)
           << ",\"joint_configuration_index\":"
           << result.joint_configuration_index
           << ",\"source_links\":[";
    for(std::size_t source = 0; source < result.source_links.size(); ++source)
      output << (source ? "," : "")
             << jsonString(result.source_links[source]);
    output << ']'
           << ",\"sphere_count\":" << result.sphere_count
           << ",\"distance_proxy_m\":" << result.distance_proxy
           << ",\"gradient_dq\":";
    jsonVector(output, result.gradient_q);
    output << ",\"hard_min_proxy_m\":" << result.hard_min_proxy
           << ",\"nearest_sphere\":" << result.nearest_sphere
           << ",\"link_smoothing_bound_m\":"
           << result.link_smoothing_bound << '}';
  }
  output << "\n  ],\n  \"spheres\":[\n";
  std::map<std::string, std::size_t> link_indices;
  for(std::size_t index = 0; index < robot.links().size(); ++index)
    link_indices.emplace(robot.links()[index].name, index);
  std::map<std::string, const rokae_demo::LinkDistanceResult*>
      group_by_source_link;
  for(const auto& result : results)
    for(const auto& source_link : result.source_links)
      if(!group_by_source_link.emplace(source_link, &result).second)
        throw std::logic_error(
            "A source link belongs to multiple joint collision groups");
  for(std::size_t index = 0; index < robot.spheres().size(); ++index)
  {
    const auto& sphere = robot.spheres()[index];
    const auto& kinematics = workspace.robot.spheres[index];
    const auto& field_sample = workspace.field_samples[index];
    output << (index ? ",\n" : "")
           << "    {\"id\":" << index
           << ",\"link\":" << jsonString(sphere.link_name)
           << ",\"center_urdf_link_m\":";
    jsonVector(output, sphere.center_link);
    const auto group = group_by_source_link.find(sphere.link_name);
    if(group == group_by_source_link.end())
      output << ",\"joint_group\":null,\"group_link\":null,"
                "\"joint_configuration_index\":null,"
                "\"center_group_link_m\":null";
    else
    {
      const auto group_link = link_indices.find(group->second->link_name);
      if(group_link == link_indices.end())
        throw std::logic_error("Joint collision group link is unknown");
      const Eigen::Vector3d center_group_link =
          workspace.robot.field_T_link[group_link->second].inverse() *
          kinematics.center_field;
      output << ",\"joint_group\":"
             << jsonString(group->second->joint_name)
             << ",\"group_link\":"
             << jsonString(group->second->link_name)
             << ",\"joint_configuration_index\":"
             << group->second->joint_configuration_index
             << ",\"center_group_link_m\":";
      jsonVector(output, center_group_link);
    }
    output
           << ",\"center_field_m\":";
    jsonVector(output, kinematics.center_field);
    output << ",\"radius_m\":" << sphere.radius
           << ",\"environment_value_m\":" << field_sample.value
           << ",\"clearance_proxy_m\":"
           << workspace.sphere_clearances[index]
           << ",\"spatial_gradient\":[" << field_sample.gradient.x << ','
           << field_sample.gradient.y << ',' << field_sample.gradient.z
           << "]}";
  }
  output << "\n  ]\n}\n";
}

}  // namespace

int main(int argc, char** argv)
{
  try
  {
    const Options options = parse(argc, argv);
    const auto robot =
        rokae_demo::SphericalRobotModel::fromUrdf(options.urdf);
    const Eigen::VectorXd default_configuration =
        jointConfiguration(options.configuration, robot.dof());
    rokae_demo::BoundaryDirectFrameOptions environment_options;
    environment_options.direct = options.direct_limits;
    environment_options.direct.intern_subexpressions = true;
    environment_options.direct.materialize_expression = false;
    environment_options.direct.retain_partition_diagnostics = false;
    environment_options.smooth = true;
    environment_options.share_subexpressions = true;
    environment_options.smoothing.share_subexpressions = true;
    environment_options.smoothing.binary_kernel = true;
    environment_options.smoothing.predecode_kernel = true;
    environment_options.smoothing.paired_batch = true;
    environment_options.smoothing.beta = options.environment_beta;
    environment_options.smoothing.max_error = options.environment_error;

    rokae_demo::BoundaryDirectFrameCache cache(1);
    rokae_demo::LinkDistanceOptions distance_options;
    distance_options.alpha = options.link_alpha;
    distance_options.max_smoothing_error = options.link_error;
    distance_options.global_margin = options.margin;
    distance_options.field_workers = options.query_workers;
    distance_options.dynamic_field_batch = options.query_workers > 4;
    rokae_demo::halfspace::BatchExecutor field_executor(
        std::min(distance_options.field_workers,robot.spheres().size()));
    std::optional<rokae_demo::halfspace::GradientBatchWorkspace>
        reusable_field_workspace;
    reusable_field_workspace.emplace(
        field_executor.makeWorkspace<rokae_demo::halfspace::ValueGradient>(
            std::min(options.direct_limits.max_nodes,std::size_t{16384}),
            true));
    std::shared_ptr<rokae_demo::BoundaryDirectFrameSnapshot> active_snapshot;
    std::optional<rokae_demo::LinkDistanceEvaluator> evaluator;
    std::optional<rokae_demo::LinkDistanceWorkspace> workspace;
    const auto run_frame = [&](const Options& frame_options,
                               const Eigen::VectorXd& configuration,
                               std::size_t frame_index) {
      const auto points = readCloud(frame_options.cloud); // I/O outside core timing
      const auto core_start = Clock::now();
      if(!options.reuse_environment) cache.clear();
      const auto environment_start = Clock::now();
      const auto prepared = cache.prepare(points, options.voxel_size,
                                          environment_options);
      const double environment_build_ms = milliseconds(environment_start);
      if(!prepared.snapshot || !prepared.snapshot->smoothTree() ||
         !prepared.snapshot->diagnostics().strict_sign_certified)
        throw std::runtime_error(
            "Environment construction did not produce a certified smooth field");

      const bool evaluator_reused = active_snapshot == prepared.snapshot;
      if(!evaluator_reused)
      {
        if(workspace && workspace->field_batch_workspace)
        {
          reusable_field_workspace.emplace(
              std::move(*workspace->field_batch_workspace));
        }
        workspace.reset();
        evaluator.reset();
        active_snapshot = prepared.snapshot;
        evaluator.emplace(robot, *active_snapshot->smoothTree(), distance_options);
        workspace.emplace(evaluator->makeWorkspace(
            std::move(*reusable_field_workspace)));
        reusable_field_workspace.reset();
      }
      rokae_demo::LinkDistancePhaseTimes phase_times;
      const auto query_start = Clock::now();
      const auto& results = options.profile_query
          ? evaluator->evaluateProfiled(configuration, options.field_T_base,
                                       *workspace, phase_times)
          : evaluator->evaluate(configuration, options.field_T_base, *workspace);
      const double query_ms = milliseconds(query_start);
      const double core_ms = milliseconds(core_start); // excludes JSON writing
      writeOutput(frame_options, robot, prepared.snapshot->occupancy(),
                  prepared.snapshot->boundary(),
                  *prepared.snapshot->smoothTree(), *evaluator, results,
                  *workspace, configuration, points.size(),
                  environment_build_ms, query_ms,
                  options.profile_query ? &phase_times : nullptr);
      if(options.stream)
      {
        std::cout << std::setprecision(17)
                  << "{\"frame\":" << frame_index
                  << ",\"cloud\":" << jsonString(frame_options.cloud.string())
                  << ",\"output\":" << jsonString(frame_options.output.string())
                  << ",\"cache_hit\":" << (prepared.cache_hit ? "true" : "false")
                  << ",\"evaluator_reused\":" << (evaluator_reused ? "true" : "false")
                  << ",\"point_count\":" << points.size()
                  << ",\"occupied_voxels\":" << prepared.snapshot->occupancy().occupied_voxels.size()
                  << ",\"occupancy_ms\":" << prepared.occupancy_ms
                  << ",\"lookup_ms\":" << prepared.lookup_ms
                  << ",\"raw_validation_ms\":" << prepared.raw_validation_ms
                  << ",\"boundary_ms\":" << prepared.boundary_ms
                  << ",\"direct_ms\":" << prepared.direct_ms
                  << ",\"environment_ms\":" << environment_build_ms
                  << ",\"query_ms\":" << query_ms
                  << ",\"core_ms\":" << core_ms
                  << ",\"core_no_validation_ms\":";
        if(prepared.cache_hit)
          std::cout << core_ms - prepared.raw_validation_ms;
        else
          std::cout << "null"; // cold/miss path has other certification work
        std::cout << "}" << std::endl;
      }
      else
        std::cout << "Wrote " << results.size() << " link distance results to "
                  << frame_options.output << "\n";
    };

    if(!options.stream)
      run_frame(options, default_configuration, 0);
    else
    {
      std::size_t line_number = 0, frame_index = 0;
      std::string line;
      while(std::getline(std::cin, line))
      {
        ++line_number;
        if(!line.empty() && line.back() == '\r') line.pop_back();
        if(line.empty() || line[0] == '#') continue;
        const auto first_tab = line.find('\t');
        const auto second_tab = first_tab == std::string::npos
            ? std::string::npos : line.find('\t', first_tab + 1);
        if(first_tab == std::string::npos || first_tab == 0 ||
           first_tab + 1 == line.size() ||
           (second_tab != std::string::npos &&
            (second_tab == first_tab + 1 || second_tab + 1 == line.size() ||
             line.find('\t', second_tab + 1) != std::string::npos)))
          throw std::invalid_argument("Invalid stream row at line " +
                                      std::to_string(line_number));
        Options frame_options = options;
        frame_options.cloud = line.substr(0, first_tab);
        frame_options.output = line.substr(
            first_tab + 1, second_tab == std::string::npos
                ? std::string::npos : second_tab - first_tab - 1);
        const Eigen::VectorXd configuration = second_tab == std::string::npos
            ? default_configuration
            : jointConfiguration(commaList(line.substr(second_tab + 1),
                                           "stream configuration"), robot.dof());
        run_frame(frame_options, configuration, frame_index++);
      }
      if(std::cin.bad()) throw std::runtime_error("Failed reading stream stdin");
      if(!frame_index) throw std::invalid_argument("Stream stdin contains no frames");
    }
    return EXIT_SUCCESS;
  }
  catch(const std::exception& error)
  {
    usage(argv[0]);
    std::cerr << "link_distance_node failed: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
