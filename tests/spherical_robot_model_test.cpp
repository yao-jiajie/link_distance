#include <rokae_demo/link_distance_evaluator.hpp>
#include <rokae_demo/spherical_robot_model.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>

namespace HS = rokae_demo::halfspace;

namespace
{
void require(bool condition, const char* message)
{
  if(!condition)
    throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance,
          const char* message)
{
  if(std::abs(actual - expected) > tolerance)
  {
    std::cerr << message << ": actual=" << actual
              << " expected=" << expected << '\n';
    throw std::runtime_error(message);
  }
}
}  // namespace

int main(int argc, char** argv)
{
  try
  {
    require(argc == 3, "Expected test and Panda URDF paths");
    const auto robot = rokae_demo::SphericalRobotModel::fromUrdf(argv[1]);
    require(robot.name() == "spherical_test_robot", "Wrong robot name");
    require(robot.dof() == 2, "Wrong DOF");
    require(robot.links().size() == 3, "Wrong link count");
    require(robot.spheres().size() == 4, "Wrong sphere count");
    require(robot.jointNames()[0] == "joint1" &&
                robot.jointNames()[1] == "joint2",
            "Wrong joint order");

    Eigen::Vector2d q;
    q << std::acos(-1.0) / 2.0, 0.25;
    auto workspace = robot.makeWorkspace();
    robot.computeSphereKinematics(q, Eigen::Isometry3d::Identity(), workspace);
    const auto& link1_first = workspace.spheres[1];
    near(link1_first.center_field.x(), 0.0, 1e-12, "link1 sphere x");
    near(link1_first.center_field.y(), 1.0, 1e-12, "link1 sphere y");
    near(link1_first.center_field.z(), 1.0, 1e-12, "link1 sphere z");
    near(link1_first.jacobian(0, 0), -1.0, 1e-12,
         "revolute Jacobian x");
    near(link1_first.jacobian(1, 0), 0.0, 1e-12,
         "revolute Jacobian y");
    const auto& link2 = workspace.spheres[3];
    near(link2.center_field.x(), 0.0, 1e-12, "link2 sphere x");
    near(link2.center_field.y(), 1.25, 1e-12, "link2 sphere y");
    near(link2.center_field.z(), 2.0, 1e-12, "link2 sphere z");
    near(link2.jacobian(0, 0), -1.25, 1e-12,
         "link2 revolute Jacobian");
    near(link2.jacobian(1, 1), 1.0, 1e-12,
         "link2 prismatic Jacobian");

    // A one-leaf SmoothTree with q=-x returns the analytic field D(x)=x.
    std::vector<HS::HalfspaceLeaf> leaves(1);
    leaves[0].id = 0;
    leaves[0].normal = {-1.0, 0.0, 0.0};
    leaves[0].offset = 0.0;
    const auto expression = HS::makeLeaf(0);
    HS::SmoothTree field(expression, leaves, {10.0, 0.0, false, false});
    rokae_demo::LinkDistanceEvaluator evaluator(
        robot, field, {12.0, 0.0, 0.01});
    auto distance_workspace = evaluator.makeWorkspace();
    Eigen::Vector2d query;
    query << 0.3, 0.1;
    const auto& results = evaluator.evaluate(
        query, Eigen::Isometry3d::Identity(), distance_workspace);
    require(results.size() == 2, "Expected one distance group per movable joint");
    require(results[0].link_name == "link1" &&
                results[0].joint_name == "joint1" &&
                results[0].sphere_count == 2 &&
                results[0].source_links == std::vector<std::string>{"link1"},
            "Wrong joint1 collision group");
    require(results[1].link_name == "link2" &&
                results[1].joint_name == "joint2" &&
                results[1].sphere_count == 1 &&
                results[1].source_links == std::vector<std::string>{"link2"},
            "Wrong joint2 collision group");
    for(const auto& result : results)
    {
      require(result.distance_proxy <= result.hard_min_proxy,
              "Smooth minimum exceeds hard minimum");
      require(result.hard_min_proxy - result.distance_proxy <=
                  result.link_smoothing_bound,
              "Smooth minimum violates error bound");
    }

    rokae_demo::LinkDistanceOptions parallel_options{12.0, 0.0, 0.01};
    parallel_options.field_workers = 4;
    rokae_demo::LinkDistanceEvaluator parallel_evaluator(
        robot, field, parallel_options);
    auto parallel_workspace = parallel_evaluator.makeWorkspace();
    const auto& parallel_results = parallel_evaluator.evaluate(
        query, Eigen::Isometry3d::Identity(), parallel_workspace);
    require(parallel_results.size() == results.size(),
            "Parallel link result count changed");
    for(std::size_t link = 0; link < results.size(); ++link)
    {
      require(parallel_results[link].distance_proxy ==
                  results[link].distance_proxy &&
              parallel_results[link].hard_min_proxy ==
                  results[link].hard_min_proxy &&
              parallel_results[link].nearest_sphere ==
                  results[link].nearest_sphere &&
              (parallel_results[link].gradient_q.array() ==
               results[link].gradient_q.array()).all(),
              "Parallel link evaluation changed a value or gradient");
    }
    HS::BatchExecutor preallocated_executor(4);
    auto preallocated_field_workspace=
        preallocated_executor.makeWorkspace<HS::ValueGradient>(
            field.nodeCount()+7);
    auto rebound_workspace=parallel_evaluator.makeWorkspace(
        std::move(preallocated_field_workspace));
    const auto& rebound_results=parallel_evaluator.evaluate(
        query, Eigen::Isometry3d::Identity(), rebound_workspace);
    for(std::size_t link = 0; link < results.size(); ++link)
      require(rebound_results[link].distance_proxy ==
                  parallel_results[link].distance_proxy &&
                  (rebound_results[link].gradient_q.array() ==
                   parallel_results[link].gradient_q.array()).all(),
              "Preallocated link workspace changed a value or gradient");
    parallel_options.dynamic_field_batch = true;
    rokae_demo::LinkDistanceEvaluator dynamic_evaluator(robot, field, parallel_options);
    auto dynamic_workspace = dynamic_evaluator.makeWorkspace();
    const auto& dynamic_results = dynamic_evaluator.evaluate(
        query, Eigen::Isometry3d::Identity(), dynamic_workspace);
    for(std::size_t link = 0; link < results.size(); ++link)
      require(dynamic_results[link].distance_proxy == results[link].distance_proxy &&
                  (dynamic_results[link].gradient_q.array() ==
                   results[link].gradient_q.array()).all(),
              "Dynamic link evaluation changed a value or gradient");
    auto profiled_workspace = parallel_evaluator.makeWorkspace();
    rokae_demo::LinkDistancePhaseTimes phase_times;
    const auto& profiled_results = parallel_evaluator.evaluateProfiled(
        query, Eigen::Isometry3d::Identity(), profiled_workspace,
        phase_times);
    require(phase_times.kinematics_ms >= 0.0 && phase_times.field_ms >= 0.0 &&
                phase_times.link_ms >= 0.0,
            "Invalid link-distance phase times");
    for(std::size_t link = 0; link < results.size(); ++link)
      require(profiled_results[link].distance_proxy == results[link].distance_proxy &&
                  (profiled_results[link].gradient_q.array() ==
                   results[link].gradient_q.array()).all(),
              "Profiled link evaluation changed a value or gradient");

    const double step = 1e-6;
    for(std::size_t link = 0; link < results.size(); ++link)
      for(Eigen::Index joint = 0; joint < query.size(); ++joint)
      {
        Eigen::Vector2d plus = query;
        Eigen::Vector2d minus = query;
        plus[joint] += step;
        minus[joint] -= step;
        auto plus_workspace = evaluator.makeWorkspace();
        auto minus_workspace = evaluator.makeWorkspace();
        const double plus_value = evaluator
                                      .evaluate(plus,
                                                Eigen::Isometry3d::Identity(),
                                                plus_workspace)[link]
                                      .distance_proxy;
        const double minus_value = evaluator
                                       .evaluate(minus,
                                                 Eigen::Isometry3d::Identity(),
                                                 minus_workspace)[link]
                                       .distance_proxy;
        const double numerical = (plus_value - minus_value) / (2.0 * step);
        near(results[link].gradient_q[joint], numerical, 2e-7,
             "End-to-end joint gradient");
      }
    const auto panda = rokae_demo::SphericalRobotModel::fromUrdf(argv[2]);
    require(panda.name() == "panda", "Wrong Panda robot name");
    require(panda.dof() == 7, "Panda must have seven movable joints");
    require(panda.spheres().size() == 59, "Panda must contain 59 spheres");
    const std::map<std::string, std::size_t> expected_spheres{
        {"panda_link0", 1},       {"panda_link1", 4},
        {"panda_link2", 4},       {"panda_link3", 4},
        {"panda_link4", 4},       {"panda_link5", 12},
        {"panda_link6", 3},       {"panda_link7", 5},
        {"panda_link8", 0},       {"panda_hand", 18},
        {"panda_leftfinger", 2},  {"panda_rightfinger", 2},
        {"panda_grasptarget", 0}};
    require(panda.links().size() == expected_spheres.size(),
            "Wrong Panda link count");
    for(const auto& link : panda.links())
    {
      const auto expected = expected_spheres.find(link.name);
      require(expected != expected_spheres.end(), "Unexpected Panda link");
      require(link.sphere_indices.size() == expected->second,
              "Wrong per-link Panda sphere count");
    }

    Eigen::VectorXd panda_q(7);
    panda_q << 0.0, -0.4, 0.0, -2.0, 0.0, 1.6, 0.8;
    auto panda_workspace = panda.makeWorkspace();
    panda.computeSphereKinematics(
        panda_q, Eigen::Isometry3d::Identity(), panda_workspace);
    require(panda_workspace.spheres.size() == 59,
            "Wrong Panda sphere kinematics count");
    std::vector<int> panda_sphere_groups(panda.spheres().size(), -1);
    std::vector<std::size_t> panda_group_sphere_counts(panda.dof(), 0);
    for(std::size_t sphere = 0; sphere < panda.spheres().size(); ++sphere)
    {
      require(panda_workspace.spheres[sphere].center_field.allFinite(),
              "Nonfinite Panda sphere center");
      require(panda_workspace.spheres[sphere].jacobian.allFinite(),
              "Nonfinite Panda sphere Jacobian");
      Eigen::Vector3d center_group_link = panda.spheres()[sphere].center_link;
      std::size_t ancestor_link = panda.spheres()[sphere].link_index;
      std::size_t group_link = ancestor_link;
      while(ancestor_link != panda.rootLink())
      {
        const int parent_joint = panda.links()[ancestor_link].parent_joint;
        require(parent_joint >= 0, "Broken Panda sphere ancestor chain");
        const auto& joint =
            panda.joints()[static_cast<std::size_t>(parent_joint)];
        if(joint.configuration_index >= 0)
        {
          panda_sphere_groups[sphere] = joint.configuration_index;
          group_link = joint.child_link;
          break;
        }
        require(joint.type == rokae_demo::RobotJointType::FIXED,
                "Unexpected joint without a configuration index");
        center_group_link = joint.parent_T_joint * center_group_link;
        ancestor_link = joint.parent_link;
      }
      const int group = panda_sphere_groups[sphere];
      if(group < 0)
      {
        require(panda.spheres()[sphere].link_name == "panda_link0" &&
                    panda_workspace.spheres[sphere].jacobian.norm() < 1e-14,
                "Only the static Panda base sphere may be ungrouped");
      }
      else
      {
        ++panda_group_sphere_counts[static_cast<std::size_t>(group)];
        require(panda.links()[group_link].name ==
                    "panda_link" + std::to_string(group + 1),
                "Sphere resolved to the wrong movable-joint child frame");
        const Eigen::Vector3d reconstructed =
            panda_workspace.field_T_link[group_link] * center_group_link;
        require((reconstructed -
                 panda_workspace.spheres[sphere].center_field)
                        .norm() < 1e-12,
                "Sphere fixed chain does not reconstruct its field center");
        for(Eigen::Index downstream = group + 1;
            downstream < panda_q.size(); ++downstream)
          require(panda_workspace.spheres[sphere]
                          .jacobian.col(downstream)
                          .norm() < 1e-14,
                  "Sphere Jacobian depends on a downstream joint");
      }
      for(Eigen::Index joint = 0; joint < panda_q.size(); ++joint)
      {
        Eigen::VectorXd plus = panda_q;
        Eigen::VectorXd minus = panda_q;
        plus[joint] += step;
        minus[joint] -= step;
        auto plus_workspace = panda.makeWorkspace();
        auto minus_workspace = panda.makeWorkspace();
        panda.computeSphereKinematics(
            plus, Eigen::Isometry3d::Identity(), plus_workspace);
        panda.computeSphereKinematics(
            minus, Eigen::Isometry3d::Identity(), minus_workspace);
        const Eigen::Vector3d numerical =
            (plus_workspace.spheres[sphere].center_field -
             minus_workspace.spheres[sphere].center_field) /
            (2.0 * step);
        require((panda_workspace.spheres[sphere].jacobian.col(joint) -
                 numerical)
                        .norm() < 5e-7,
                "Panda sphere Jacobian finite-difference failure");
      }
    }

    rokae_demo::LinkDistanceEvaluator panda_evaluator(
        panda, field, {20.0, 0.0, 0.005});
    auto panda_distance_workspace = panda_evaluator.makeWorkspace();
    const auto& panda_results = panda_evaluator.evaluate(
        panda_q, Eigen::Isometry3d::Identity(), panda_distance_workspace);
    require(panda_results.size() == 7,
            "Panda must return one group per arm joint");
    const std::vector<std::size_t> expected_group_spheres{
        4, 4, 4, 4, 12, 3, 27};
    require(panda_group_sphere_counts == expected_group_spheres,
            "Per-sphere joint-frame audit produced wrong group counts");
    for(std::size_t group = 0; group < panda_results.size(); ++group)
    {
      require(panda_results[group].link_name ==
                  "panda_link" + std::to_string(group + 1) &&
                  panda_results[group].joint_name ==
                      "panda_joint" + std::to_string(group + 1) &&
                  panda_results[group].joint_configuration_index ==
                      static_cast<int>(group) &&
                  panda_results[group].sphere_count ==
                      expected_group_spheres[group],
              "Wrong Panda joint collision group");
    }
    require(panda_results.back().source_links ==
                std::vector<std::string>{"panda_link7", "panda_hand",
                                         "panda_leftfinger",
                                         "panda_rightfinger"},
            "Fixed Panda tool links were not merged into link7");
    for(std::size_t sphere = 0; sphere < panda.spheres().size(); ++sphere)
      if(panda_sphere_groups[sphere] >= 0)
      {
        const auto& group = panda_results[static_cast<std::size_t>(
            panda_sphere_groups[sphere])];
        require(std::find(group.source_links.begin(), group.source_links.end(),
                          panda.spheres()[sphere].link_name) !=
                    group.source_links.end(),
                "Sphere owner link is absent from its joint group");
      }
    for(std::size_t link = 0; link < panda_results.size(); ++link)
    {
      require(panda_results[link].distance_proxy <=
                  panda_results[link].hard_min_proxy,
              "Panda smooth minimum exceeds hard minimum");
      require(panda_results[link].hard_min_proxy -
                      panda_results[link].distance_proxy <=
                  panda_results[link].link_smoothing_bound,
              "Panda smooth minimum violates error bound");
      for(Eigen::Index joint = 0; joint < panda_q.size(); ++joint)
      {
        Eigen::VectorXd plus = panda_q;
        Eigen::VectorXd minus = panda_q;
        plus[joint] += step;
        minus[joint] -= step;
        auto plus_workspace = panda_evaluator.makeWorkspace();
        auto minus_workspace = panda_evaluator.makeWorkspace();
        const double plus_value = panda_evaluator
                                      .evaluate(plus,
                                                Eigen::Isometry3d::Identity(),
                                                plus_workspace)[link]
                                      .distance_proxy;
        const double minus_value = panda_evaluator
                                       .evaluate(minus,
                                                 Eigen::Isometry3d::Identity(),
                                                 minus_workspace)[link]
                                       .distance_proxy;
        const double numerical = (plus_value - minus_value) / (2.0 * step);
        near(panda_results[link].gradient_q[joint], numerical, 5e-7,
             "Panda link gradient finite-difference failure");
      }
    }
    return EXIT_SUCCESS;
  }
  catch(const std::exception& error)
  {
    std::cerr << "spherical_robot_model_test failed: " << error.what()
              << '\n';
    return EXIT_FAILURE;
  }
}
