#include <rokae_demo/link_distance_evaluator.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace rokae_demo
{

LinkDistanceEvaluator::LinkDistanceEvaluator(
    const SphericalRobotModel& robot,
    const halfspace::SmoothTree& field,
    const LinkDistanceOptions& options)
    : robot_(&robot), field_(&field), margin_(options.global_margin),
      dynamic_field_batch_(options.dynamic_field_batch)
{
  if(!std::isfinite(options.alpha) ||
     !std::isfinite(options.max_smoothing_error) ||
     !std::isfinite(options.global_margin) || options.alpha < 0.0 ||
     options.max_smoothing_error < 0.0 || options.global_margin < 0.0 ||
     options.field_workers == 0 ||
     ((options.alpha > 0.0) == (options.max_smoothing_error > 0.0)))
    throw std::invalid_argument(
        "Choose one positive link alpha/error and a nonnegative margin");

  field_workers_ = std::min(options.field_workers, robot.spheres().size());

  std::vector<JointCollisionGroup> groups_by_configuration(robot.dof());
  for(const auto& joint : robot.joints())
    if(joint.configuration_index >= 0)
    {
      auto& group = groups_by_configuration[static_cast<std::size_t>(
          joint.configuration_index)];
      group.link_name = robot.links()[joint.child_link].name;
      group.joint_name = joint.name;
      group.joint_configuration_index = joint.configuration_index;
    }

  for(std::size_t sphere = 0; sphere < robot.spheres().size(); ++sphere)
  {
    const auto& collision = robot.spheres()[sphere];
    std::size_t link = collision.link_index;
    int configuration_index = -1;
    while(link != robot.rootLink())
    {
      const int parent_joint = robot.links()[link].parent_joint;
      if(parent_joint < 0)
        throw std::logic_error("Broken ancestor chain in spherical robot");
      const auto& joint =
          robot.joints()[static_cast<std::size_t>(parent_joint)];
      if(joint.configuration_index >= 0)
      {
        configuration_index = joint.configuration_index;
        break;
      }
      link = joint.parent_link;
    }
    // Root/static-base geometry has no controlling movable joint. It remains
    // available in sphere diagnostics but is not a joint-link distance group.
    if(configuration_index < 0)
      continue;
    auto& group = groups_by_configuration[static_cast<std::size_t>(
        configuration_index)];
    group.sphere_indices.push_back(sphere);
    if(std::find(group.source_links.begin(), group.source_links.end(),
                 collision.link_name) == group.source_links.end())
      group.source_links.push_back(collision.link_name);
  }

  std::size_t maximum_spheres = 0;
  for(auto& group : groups_by_configuration)
    if(!group.sphere_indices.empty())
    {
      maximum_spheres = std::max(maximum_spheres,
                                 group.sphere_indices.size());
      groups_.push_back(std::move(group));
    }
  if(groups_.empty())
    throw std::invalid_argument(
        "Robot model has no collision spheres controlled by movable joints");

  const double weight = maximum_spheres > 1
                            ? std::log(static_cast<double>(maximum_spheres))
                            : 0.0;
  alpha_ = options.alpha > 0.0
               ? options.alpha
               : std::nextafter((weight > 0.0 ? weight : 1.0) /
                                    options.max_smoothing_error,
                                std::numeric_limits<double>::infinity());
  maximum_error_ = weight > 0.0
                       ? std::nextafter(weight / alpha_,
                                        std::numeric_limits<double>::infinity())
                       : 0.0;
  if(!(alpha_ > 0.0) || !std::isfinite(alpha_) ||
     !std::isfinite(maximum_error_))
    throw std::invalid_argument("Unrepresentable link smoothing parameter");
}

LinkDistanceWorkspace LinkDistanceEvaluator::makeWorkspace() const
{ return makeWorkspaceImpl(nullptr); }

LinkDistanceWorkspace LinkDistanceEvaluator::makeWorkspace(
    const halfspace::BatchExecutor& executor) const
{
  if(executor.workerCount()!=field_workers_)
    throw std::invalid_argument("Link distance executor has the wrong worker count");
  return makeWorkspaceImpl(&executor);
}


LinkDistanceWorkspace LinkDistanceEvaluator::makeWorkspace(
    halfspace::GradientBatchWorkspace&& preallocated) const
{
  if(preallocated.workerCount()!=field_workers_)
    throw std::invalid_argument(
        "Link distance preallocated workspace has the wrong worker count");
  field_->rebindGradientBatchWorkspace(preallocated);
  return makeWorkspaceImpl(nullptr,&preallocated);
}

LinkDistanceWorkspace LinkDistanceEvaluator::makeWorkspaceImpl(
    const halfspace::BatchExecutor* executor,
    halfspace::GradientBatchWorkspace* preallocated) const
{
  LinkDistanceWorkspace workspace;
  workspace.robot = robot_->makeWorkspace();
  workspace.field_points.resize(robot_->spheres().size());
  workspace.field_samples.resize(robot_->spheres().size());
  if(preallocated)
    workspace.field_batch_workspace.emplace(std::move(*preallocated));
  else if(executor)
    workspace.field_batch_workspace.emplace(
        field_->makeGradientBatchWorkspace(*executor));
  else
    workspace.field_batch_workspace.emplace(
        field_->makeGradientBatchWorkspace(field_workers_));
  workspace.sphere_clearances.resize(robot_->spheres().size());
  workspace.results.resize(groups_.size());
  for(std::size_t index = 0; index < groups_.size(); ++index)
  {
    const auto& group = groups_[index];
    auto& result = workspace.results[index];
    result.link_name = group.link_name;
    result.joint_name = group.joint_name;
    result.joint_configuration_index = group.joint_configuration_index;
    result.source_links = group.source_links;
    result.sphere_count = group.sphere_indices.size();
    result.gradient_q = Eigen::VectorXd::Zero(
        static_cast<Eigen::Index>(robot_->dof()));
  }
  return workspace;
}

const std::vector<LinkDistanceResult>& LinkDistanceEvaluator::evaluate(
    const Eigen::Ref<const Eigen::VectorXd>& configuration,
    const Eigen::Isometry3d& field_T_base,
    LinkDistanceWorkspace& workspace) const
{
  return evaluateImpl<false>(configuration, field_T_base, workspace, nullptr);
}

const std::vector<LinkDistanceResult>& LinkDistanceEvaluator::evaluateProfiled(
    const Eigen::Ref<const Eigen::VectorXd>& configuration,
    const Eigen::Isometry3d& field_T_base,
    LinkDistanceWorkspace& workspace,
    LinkDistancePhaseTimes& times) const
{
  return evaluateImpl<true>(configuration, field_T_base, workspace, &times);
}

template<bool Profile>
const std::vector<LinkDistanceResult>& LinkDistanceEvaluator::evaluateImpl(
    const Eigen::Ref<const Eigen::VectorXd>& configuration,
    const Eigen::Isometry3d& field_T_base,
    LinkDistanceWorkspace& workspace,
    LinkDistancePhaseTimes* times) const
{
  using Clock = std::chrono::steady_clock;
  Clock::time_point phase_start;
  if constexpr(Profile)
    phase_start = Clock::now();
  if(workspace.field_points.size() != robot_->spheres().size() ||
     workspace.field_samples.size() != robot_->spheres().size() ||
     workspace.sphere_clearances.size() != robot_->spheres().size() ||
     workspace.results.size() != groups_.size() ||
     !workspace.field_batch_workspace ||
     workspace.field_batch_workspace->workerCount() != field_workers_ ||
     workspace.field_batch_workspace->valuesPerWorker() != field_->nodeCount())
    throw std::invalid_argument("Link distance workspace has the wrong size");

  robot_->computeSphereKinematics(configuration, field_T_base, workspace.robot);
  for(std::size_t index = 0; index < robot_->spheres().size(); ++index)
  {
    const Eigen::Vector3d& center = workspace.robot.spheres[index].center_field;
    if(!center.allFinite())
      throw std::overflow_error("Nonfinite sphere center");
    workspace.field_points[index] = {center.x(), center.y(), center.z()};
  }
  if constexpr(Profile)
  {
    const auto now = Clock::now();
    times->kinematics_ms = std::chrono::duration<double, std::milli>(
        now - phase_start).count();
    phase_start = now;
  }
  field_->evaluateWithGradientTrustedBatch(
      workspace.field_points, workspace.field_samples,
      *workspace.field_batch_workspace, dynamic_field_batch_);
  if constexpr(Profile)
  {
    const auto now = Clock::now();
    times->field_ms = std::chrono::duration<double, std::milli>(
        now - phase_start).count();
    phase_start = now;
  }
  for(std::size_t index = 0; index < robot_->spheres().size(); ++index)
  {
    const auto& sample = workspace.field_samples[index];
    if(!std::isfinite(sample.value) || !std::isfinite(sample.gradient.x) ||
       !std::isfinite(sample.gradient.y) || !std::isfinite(sample.gradient.z))
      throw std::overflow_error("Nonfinite smooth field sample");
    const double clearance = sample.value -
                             robot_->spheres()[index].radius - margin_;
    if(!std::isfinite(clearance))
      throw std::overflow_error("Nonfinite sphere clearance proxy");
    workspace.sphere_clearances[index] = clearance;
  }

  for(std::size_t result_index = 0; result_index < groups_.size();
      ++result_index)
  {
    const auto& group = groups_[result_index];
    auto& result = workspace.results[result_index];
    const auto minimum = std::min_element(
        group.sphere_indices.begin(), group.sphere_indices.end(),
        [&](std::size_t lhs, std::size_t rhs) {
          return workspace.sphere_clearances[lhs] <
                 workspace.sphere_clearances[rhs];
        });
    result.nearest_sphere = *minimum;
    result.hard_min_proxy = workspace.sphere_clearances[*minimum];
    result.gradient_q.setZero();

    double weight_sum = 0.0;
    for(const std::size_t sphere : group.sphere_indices)
    {
      const double weight = std::exp(
          -alpha_ * (workspace.sphere_clearances[sphere] -
                     result.hard_min_proxy));
      weight_sum += weight;
      const auto& spatial = workspace.field_samples[sphere].gradient;
      const Eigen::Vector3d spatial_gradient(spatial.x, spatial.y, spatial.z);
      result.gradient_q.noalias() +=
          weight * workspace.robot.spheres[sphere].jacobian.transpose() *
          spatial_gradient;
    }
    if(!(weight_sum >= 1.0) || !std::isfinite(weight_sum))
      throw std::runtime_error("Invalid link smooth-min weight sum");
    result.distance_proxy =
        result.hard_min_proxy - std::log(weight_sum) / alpha_;
    result.gradient_q /= weight_sum;
    result.link_smoothing_bound =
        group.sphere_indices.size() > 1
            ? std::nextafter(
                  std::log(static_cast<double>(group.sphere_indices.size())) /
                      alpha_,
                  std::numeric_limits<double>::infinity())
            : 0.0;
    // The analytic bound applies before the final double subtraction.
    // Cancellation in hard_min - distance can exceed the rounded bound by
    // a few ulps even for an exactly tied group of spheres.
    const double rounding_slack =
        8.0 * std::numeric_limits<double>::epsilon() *
        std::max({1.0, std::abs(result.hard_min_proxy),
                  std::abs(result.distance_proxy)});
    if(!std::isfinite(result.distance_proxy) ||
       !result.gradient_q.allFinite() ||
       result.distance_proxy > result.hard_min_proxy ||
       result.hard_min_proxy - result.distance_proxy >
           result.link_smoothing_bound + rounding_slack)
      throw std::runtime_error("Invalid link smooth-min value or gradient");
  }
  if constexpr(Profile)
    times->link_ms = std::chrono::duration<double, std::milli>(
        Clock::now() - phase_start).count();
  return workspace.results;
}

}  // namespace rokae_demo
