#pragma once

#include <rokae_demo/halfspace_smooth_tree.hpp>
#include <rokae_demo/spherical_robot_model.hpp>

#include <Eigen/Core>

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace rokae_demo
{

struct LinkDistanceOptions
{
  double alpha = 0.0;               // inverse meters; choose exactly one
  double max_smoothing_error = 0.0; // meters
  double global_margin = 0.0;       // meters
  std::size_t field_workers = 1;    // independent sphere-field queries
  bool dynamic_field_batch = false;
};

struct LinkDistanceResult
{
  std::string link_name;
  std::string joint_name;
  int joint_configuration_index = -1;
  std::vector<std::string> source_links;
  double distance_proxy = 0.0;
  Eigen::VectorXd gradient_q;
  std::size_t sphere_count = 0;
  double hard_min_proxy = 0.0;
  std::size_t nearest_sphere = 0; // global sphere index
  double link_smoothing_bound = 0.0;
};

struct LinkDistanceWorkspace
{
  RobotKinematicsWorkspace robot;
  std::vector<halfspace::Vec3> field_points;
  std::vector<halfspace::ValueGradient> field_samples;
  std::optional<halfspace::GradientBatchWorkspace> field_batch_workspace;
  std::vector<double> sphere_clearances;
  std::vector<LinkDistanceResult> results;
};

struct LinkDistancePhaseTimes
{
  double kinematics_ms = 0.0; // includes sphere-center gathering
  double field_ms = 0.0;      // batch LSE value and spatial gradient
  double link_ms = 0.0;       // clearances and link smooth-min reduction
};

class LinkDistanceEvaluator
{
public:
  LinkDistanceEvaluator(const SphericalRobotModel& robot,
                        const halfspace::SmoothTree& field,
                        const LinkDistanceOptions& options);

  LinkDistanceWorkspace makeWorkspace() const;
  // Reuses field-independent worker threads initialized before a dynamic
  // environment arrives. Workspaces sharing it must evaluate sequentially.
  LinkDistanceWorkspace makeWorkspace(
      const halfspace::BatchExecutor& executor) const;
  // Rebinds field-independent worker buffers while retaining their capacity.
  LinkDistanceWorkspace makeWorkspace(
      halfspace::GradientBatchWorkspace&& preallocated) const;

  const std::vector<LinkDistanceResult>& evaluate(
      const Eigen::Ref<const Eigen::VectorXd>& configuration,
      const Eigen::Isometry3d& field_T_base,
      LinkDistanceWorkspace& workspace) const;

  // Diagnostic variant only; per-phase clocks are excluded from the normal
  // evaluate() path and from production query_ms benchmarks.
  const std::vector<LinkDistanceResult>& evaluateProfiled(
      const Eigen::Ref<const Eigen::VectorXd>& configuration,
      const Eigen::Isometry3d& field_T_base,
      LinkDistanceWorkspace& workspace,
      LinkDistancePhaseTimes& times) const;

  double alpha() const { return alpha_; }
  double maximumLinkSmoothingBound() const { return maximum_error_; }
  std::size_t fieldWorkerCount() const { return field_workers_; }

private:
  LinkDistanceWorkspace makeWorkspaceImpl(
      const halfspace::BatchExecutor* executor,
      halfspace::GradientBatchWorkspace* preallocated = nullptr) const;
  template<bool Profile>
  const std::vector<LinkDistanceResult>& evaluateImpl(
      const Eigen::Ref<const Eigen::VectorXd>& configuration,
      const Eigen::Isometry3d& field_T_base,
      LinkDistanceWorkspace& workspace,
      LinkDistancePhaseTimes* times) const;
  struct JointCollisionGroup
  {
    std::string link_name;
    std::string joint_name;
    int joint_configuration_index = -1;
    std::vector<std::string> source_links;
    std::vector<std::size_t> sphere_indices;
  };

  const SphericalRobotModel* robot_ = nullptr;
  const halfspace::SmoothTree* field_ = nullptr;
  double alpha_ = 0.0;
  double maximum_error_ = 0.0;
  double margin_ = 0.0;
  std::size_t field_workers_ = 1;
  bool dynamic_field_batch_ = false;
  std::vector<JointCollisionGroup> groups_;
};

}  // namespace rokae_demo
