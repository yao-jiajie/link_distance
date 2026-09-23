#include <rokae_demo/spherical_robot_model.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace rokae_demo
{
namespace
{
namespace PT = boost::property_tree;

double finiteDouble(const std::string& text, const char* description)
{
  std::size_t parsed = 0;
  double value = 0.0;
  try
  {
    value = std::stod(text, &parsed);
  }
  catch(const std::exception&)
  {
    throw std::runtime_error(std::string("Invalid ") + description + ": " + text);
  }
  if(parsed != text.size() || !std::isfinite(value))
    throw std::runtime_error(std::string("Invalid ") + description + ": " + text);
  return value;
}

Eigen::Vector3d vector3(const std::string& text, const char* description)
{
  std::istringstream input(text);
  Eigen::Vector3d result;
  if(!(input >> result.x() >> result.y() >> result.z()))
    throw std::runtime_error(std::string("Invalid ") + description + ": " + text);
  input >> std::ws;
  if(!input.eof() || !result.allFinite())
    throw std::runtime_error(std::string("Invalid ") + description + ": " + text);
  return result;
}

Eigen::Isometry3d origin(const PT::ptree& node)
{
  Eigen::Isometry3d result = Eigen::Isometry3d::Identity();
  const auto child = node.get_child_optional("origin");
  if(!child)
    return result;
  const Eigen::Vector3d xyz = vector3(
      child->get<std::string>("<xmlattr>.xyz", "0 0 0"), "origin xyz");
  const Eigen::Vector3d rpy = vector3(
      child->get<std::string>("<xmlattr>.rpy", "0 0 0"), "origin rpy");
  result.linear() =
      (Eigen::AngleAxisd(rpy.z(), Eigen::Vector3d::UnitZ()) *
       Eigen::AngleAxisd(rpy.y(), Eigen::Vector3d::UnitY()) *
       Eigen::AngleAxisd(rpy.x(), Eigen::Vector3d::UnitX()))
          .toRotationMatrix();
  result.translation() = xyz;
  return result;
}

RobotJointType jointType(const std::string& text)
{
  if(text == "fixed")
    return RobotJointType::FIXED;
  if(text == "revolute")
    return RobotJointType::REVOLUTE;
  if(text == "continuous")
    return RobotJointType::CONTINUOUS;
  if(text == "prismatic")
    return RobotJointType::PRISMATIC;
  throw std::runtime_error("Unsupported URDF joint type: " + text);
}

Eigen::Isometry3d jointMotion(RobotJointType type,
                              const Eigen::Vector3d& axis,
                              double value)
{
  Eigen::Isometry3d result = Eigen::Isometry3d::Identity();
  if(type == RobotJointType::REVOLUTE || type == RobotJointType::CONTINUOUS)
    result.linear() = Eigen::AngleAxisd(value, axis).toRotationMatrix();
  else if(type == RobotJointType::PRISMATIC)
    result.translation() = value * axis;
  return result;
}

}  // namespace

SphericalRobotModel SphericalRobotModel::fromUrdf(
    const std::filesystem::path& path)
{
  PT::ptree document;
  try
  {
    PT::read_xml(path.string(), document,
                 PT::xml_parser::trim_whitespace |
                     PT::xml_parser::no_comments);
  }
  catch(const std::exception& error)
  {
    throw std::runtime_error("Cannot parse URDF " + path.string() + ": " +
                             error.what());
  }

  const auto robot_node = document.get_child_optional("robot");
  if(!robot_node)
    throw std::runtime_error("URDF has no robot root element");

  SphericalRobotModel model;
  model.name_ = robot_node->get<std::string>("<xmlattr>.name", "");
  if(model.name_.empty())
    throw std::runtime_error("URDF robot name is empty");

  std::map<std::string, std::size_t> link_index;
  for(const auto& entry : *robot_node)
  {
    if(entry.first != "link")
      continue;
    const std::string name =
        entry.second.get<std::string>("<xmlattr>.name", "");
    if(name.empty() || link_index.count(name))
      throw std::runtime_error("URDF contains an empty or duplicate link name");
    const std::size_t index = model.links_.size();
    link_index.emplace(name, index);
    model.links_.push_back({name, -1, {}});
  }
  if(model.links_.empty())
    throw std::runtime_error("URDF contains no links");

  for(const auto& entry : *robot_node)
  {
    if(entry.first != "link")
      continue;
    const std::string name =
        entry.second.get<std::string>("<xmlattr>.name");
    const std::size_t index = link_index.at(name);
    for(const auto& link_child : entry.second)
    {
      if(link_child.first != "collision")
        continue;
      const auto geometry = link_child.second.get_child_optional("geometry");
      if(!geometry)
        throw std::runtime_error("Collision on link " + name +
                                 " has no geometry");
      std::size_t geometry_count = 0;
      double radius = 0.0;
      for(const auto& shape : *geometry)
      {
        if(shape.first == "<xmlattr>")
          continue;
        ++geometry_count;
        if(shape.first != "sphere")
          throw std::runtime_error("Non-sphere collision geometry on link " +
                                   name + ": " + shape.first);
        radius = finiteDouble(
            shape.second.get<std::string>("<xmlattr>.radius", ""),
            "sphere radius");
      }
      if(geometry_count != 1 || !(radius > 0.0))
        throw std::runtime_error("Collision on link " + name +
                                 " must contain one positive sphere");

      RobotCollisionSphere sphere;
      sphere.link_name = name;
      sphere.link_index = index;
      sphere.center_link = origin(link_child.second).translation();
      sphere.radius = radius;
      model.links_[index].sphere_indices.push_back(model.spheres_.size());
      model.spheres_.push_back(std::move(sphere));
    }
  }
  if(model.spheres_.empty())
    throw std::runtime_error("URDF contains no collision spheres");

  for(const auto& entry : *robot_node)
  {
    if(entry.first != "joint")
      continue;
    SphericalRobotJoint joint;
    joint.name = entry.second.get<std::string>("<xmlattr>.name", "");
    if(joint.name.empty())
      throw std::runtime_error("URDF contains an unnamed joint");
    if(std::any_of(model.joints_.begin(), model.joints_.end(),
                   [&](const auto& existing) {
                     return existing.name == joint.name;
                   }))
      throw std::runtime_error("Duplicate URDF joint name: " + joint.name);
    joint.type = jointType(
        entry.second.get<std::string>("<xmlattr>.type", ""));
    const auto parent = entry.second.get_child_optional("parent");
    const auto child = entry.second.get_child_optional("child");
    if(!parent || !child)
      throw std::runtime_error("Joint " + joint.name +
                               " has no parent or child");
    const std::string parent_name =
        parent->get<std::string>("<xmlattr>.link", "");
    const std::string child_name =
        child->get<std::string>("<xmlattr>.link", "");
    if(!link_index.count(parent_name) || !link_index.count(child_name))
      throw std::runtime_error("Joint " + joint.name +
                               " references an unknown link");
    joint.parent_link = link_index.at(parent_name);
    joint.child_link = link_index.at(child_name);
    if(model.links_[joint.child_link].parent_joint >= 0)
      throw std::runtime_error("Link " + child_name +
                               " has multiple parent joints");
    joint.parent_T_joint = origin(entry.second);

    const auto axis = entry.second.get_child_optional("axis");
    joint.axis_joint =
        axis ? vector3(axis->get<std::string>("<xmlattr>.xyz", "1 0 0"),
                       "joint axis")
             : Eigen::Vector3d::UnitX();
    if(joint.type != RobotJointType::FIXED)
    {
      const double norm = joint.axis_joint.norm();
      if(!(norm > std::numeric_limits<double>::epsilon()) ||
         !std::isfinite(norm))
        throw std::runtime_error("Movable joint has an invalid axis: " +
                                 joint.name);
      joint.axis_joint /= norm;
      if(entry.second.get_child_optional("mimic"))
        throw std::runtime_error(
            "Movable mimic joints are not supported: " + joint.name);
      joint.configuration_index =
          static_cast<int>(model.joint_names_.size());
      model.joint_names_.push_back(joint.name);
    }

    model.links_[joint.child_link].parent_joint =
        static_cast<int>(model.joints_.size());
    model.joints_.push_back(std::move(joint));
  }

  std::vector<std::size_t> roots;
  for(std::size_t index = 0; index < model.links_.size(); ++index)
    if(model.links_[index].parent_joint < 0)
      roots.push_back(index);
  if(roots.size() != 1)
    throw std::runtime_error("URDF must contain exactly one root link");
  model.root_link_ = roots.front();

  std::vector<std::vector<std::size_t>> children(model.links_.size());
  for(std::size_t joint = 0; joint < model.joints_.size(); ++joint)
    children[model.joints_[joint].parent_link].push_back(
        model.joints_[joint].child_link);
  std::deque<std::size_t> queue{model.root_link_};
  std::set<std::size_t> visited;
  while(!queue.empty())
  {
    const std::size_t link = queue.front();
    queue.pop_front();
    if(!visited.insert(link).second)
      throw std::runtime_error("URDF kinematic graph contains a cycle");
    model.topological_links_.push_back(link);
    for(const std::size_t child : children[link])
      queue.push_back(child);
  }
  if(visited.size() != model.links_.size())
    throw std::runtime_error("URDF kinematic graph is disconnected");
  return model;
}

RobotKinematicsWorkspace SphericalRobotModel::makeWorkspace() const
{
  RobotKinematicsWorkspace workspace;
  workspace.field_T_link.resize(links_.size(), Eigen::Isometry3d::Identity());
  workspace.joint_origin_field.resize(joints_.size(), Eigen::Vector3d::Zero());
  workspace.joint_axis_field.resize(joints_.size(), Eigen::Vector3d::Zero());
  workspace.spheres.resize(spheres_.size());
  for(auto& sphere : workspace.spheres)
    sphere.jacobian = Eigen::MatrixXd::Zero(3, static_cast<Eigen::Index>(dof()));
  return workspace;
}

void SphericalRobotModel::computeSphereKinematics(
    const Eigen::Ref<const Eigen::VectorXd>& configuration,
    const Eigen::Isometry3d& field_T_base,
    RobotKinematicsWorkspace& workspace) const
{
  if(static_cast<std::size_t>(configuration.size()) != dof() ||
     !configuration.allFinite())
    throw std::invalid_argument("Robot configuration has the wrong size or nonfinite values");
  if(!field_T_base.matrix().allFinite())
    throw std::invalid_argument("Robot base transform is nonfinite");
  if(workspace.field_T_link.size() != links_.size() ||
     workspace.joint_origin_field.size() != joints_.size() ||
     workspace.joint_axis_field.size() != joints_.size() ||
     workspace.spheres.size() != spheres_.size())
    throw std::invalid_argument("Robot kinematics workspace has the wrong size");

  workspace.field_T_link[root_link_] = field_T_base;
  for(const std::size_t link_index : topological_links_)
  {
    if(link_index == root_link_)
      continue;
    const int joint_index = links_[link_index].parent_joint;
    if(joint_index < 0)
      throw std::logic_error("Non-root link has no parent joint");
    const auto& joint = joints_[static_cast<std::size_t>(joint_index)];
    const Eigen::Isometry3d field_T_joint =
        workspace.field_T_link[joint.parent_link] * joint.parent_T_joint;
    workspace.joint_origin_field[static_cast<std::size_t>(joint_index)] =
        field_T_joint.translation();
    workspace.joint_axis_field[static_cast<std::size_t>(joint_index)] =
        field_T_joint.linear() * joint.axis_joint;
    const double value = joint.configuration_index >= 0
                             ? configuration[joint.configuration_index]
                             : 0.0;
    workspace.field_T_link[link_index] =
        field_T_joint * jointMotion(joint.type, joint.axis_joint, value);
  }

  for(std::size_t sphere_index = 0; sphere_index < spheres_.size();
      ++sphere_index)
  {
    const auto& sphere = spheres_[sphere_index];
    auto& output = workspace.spheres[sphere_index];
    if(output.jacobian.rows() != 3 ||
       static_cast<std::size_t>(output.jacobian.cols()) != dof())
      throw std::invalid_argument("Sphere Jacobian workspace has the wrong size");
    output.center_field = workspace.field_T_link[sphere.link_index] *
                          sphere.center_link;
    output.jacobian.setZero();

    std::size_t link = sphere.link_index;
    while(link != root_link_)
    {
      const int joint_index = links_[link].parent_joint;
      if(joint_index < 0)
        throw std::logic_error("Broken ancestor chain in robot model");
      const auto& joint = joints_[static_cast<std::size_t>(joint_index)];
      if(joint.configuration_index >= 0)
      {
        Eigen::Vector3d derivative;
        if(joint.type == RobotJointType::PRISMATIC)
          derivative = workspace.joint_axis_field[static_cast<std::size_t>(joint_index)];
        else
          derivative =
              workspace.joint_axis_field[static_cast<std::size_t>(joint_index)]
                  .cross(output.center_field -
                         workspace.joint_origin_field[static_cast<std::size_t>(joint_index)]);
        output.jacobian.col(joint.configuration_index) = derivative;
      }
      link = joint.parent_link;
    }
  }
}

}  // namespace rokae_demo
