#pragma once

#include <Eigen/Geometry>

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace rokae_demo
{

enum class RobotJointType
{
  FIXED,
  REVOLUTE,
  CONTINUOUS,
  PRISMATIC
};

struct RobotCollisionSphere
{
  std::string link_name;
  std::size_t link_index = 0;
  Eigen::Vector3d center_link = Eigen::Vector3d::Zero();
  double radius = 0.0;
};

struct SphericalRobotLink
{
  std::string name;
  int parent_joint = -1;
  std::vector<std::size_t> sphere_indices;
};

struct SphericalRobotJoint
{
  std::string name;
  RobotJointType type = RobotJointType::FIXED;
  std::size_t parent_link = 0;
  std::size_t child_link = 0;
  Eigen::Isometry3d parent_T_joint = Eigen::Isometry3d::Identity();
  Eigen::Vector3d axis_joint = Eigen::Vector3d::UnitX();
  int configuration_index = -1;
};

struct SphereKinematics
{
  Eigen::Vector3d center_field = Eigen::Vector3d::Zero();
  Eigen::Matrix<double, 3, Eigen::Dynamic> jacobian;
};

struct RobotKinematicsWorkspace
{
  std::vector<Eigen::Isometry3d> field_T_link;
  std::vector<Eigen::Vector3d> joint_origin_field;
  std::vector<Eigen::Vector3d> joint_axis_field;
  std::vector<SphereKinematics> spheres;
};

class SphericalRobotModel
{
public:
  static SphericalRobotModel fromUrdf(const std::filesystem::path& path);

  std::size_t dof() const { return joint_names_.size(); }
  std::size_t rootLink() const { return root_link_; }
  const std::string& name() const { return name_; }
  const std::vector<std::string>& jointNames() const { return joint_names_; }
  const std::vector<SphericalRobotLink>& links() const { return links_; }
  const std::vector<SphericalRobotJoint>& joints() const { return joints_; }
  const std::vector<RobotCollisionSphere>& spheres() const { return spheres_; }

  RobotKinematicsWorkspace makeWorkspace() const;

  void computeSphereKinematics(
      const Eigen::Ref<const Eigen::VectorXd>& configuration,
      const Eigen::Isometry3d& field_T_base,
      RobotKinematicsWorkspace& workspace) const;

private:
  std::string name_;
  std::vector<std::string> joint_names_;
  std::vector<SphericalRobotLink> links_;
  std::vector<SphericalRobotJoint> joints_;
  std::vector<RobotCollisionSphere> spheres_;
  std::vector<std::size_t> topological_links_;
  std::size_t root_link_ = 0;
};

}  // namespace rokae_demo
