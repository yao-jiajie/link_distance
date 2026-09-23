#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace rokae_demo::halfspace
{

struct Vec3
{
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

struct HalfspaceLeaf
{
  int id = -1;
  Vec3 normal;
  Vec3 point;
  double offset = 0.0;
  std::vector<int> convex_part_ids;
  std::vector<int> source_face_ids;
  std::string source_type;
};

enum class TreeOperator
{
  LEAF,
  MIN,
  MAX
};

struct Expression;
using ExpressionPtr = std::shared_ptr<const Expression>;

struct Expression
{
  TreeOperator op = TreeOperator::LEAF;
  int leaf_id = -1;
  std::vector<ExpressionPtr> children;
};

struct SerializedNode
{
  int id = -1;
  TreeOperator op = TreeOperator::LEAF;
  int leaf_id = -1;
  std::vector<int> children;
};

struct SerializedTree
{
  std::vector<SerializedNode> nodes;
  int root_id = -1;
};

inline ExpressionPtr makeLeaf(int leaf_id)
{
  auto expression = std::make_shared<Expression>();
  expression->leaf_id = leaf_id;
  return expression;
}

inline ExpressionPtr makeNode(TreeOperator op,
                              std::vector<ExpressionPtr> children)
{
  if(op == TreeOperator::LEAF)
    throw std::invalid_argument("A LEAF node must be created with makeLeaf");
  if(children.empty())
    throw std::invalid_argument("A MIN/MAX node must have children");
  auto expression = std::make_shared<Expression>();
  expression->op = op;
  expression->children = std::move(children);
  return expression;
}

inline const char* operatorName(TreeOperator op)
{
  switch(op)
  {
    case TreeOperator::LEAF:
      return "leaf";
    case TreeOperator::MIN:
      return "min";
    case TreeOperator::MAX:
      return "max";
  }
  return "unknown";
}

inline const char* booleanAlias(TreeOperator op)
{
  switch(op)
  {
    case TreeOperator::LEAF:
      return "halfspace";
    case TreeOperator::MIN:
      return "OR";
    case TreeOperator::MAX:
      return "AND";
  }
  return "unknown";
}

inline double evaluate(const ExpressionPtr& expression,
                       const std::vector<HalfspaceLeaf>& leaves,
                       const Vec3& point)
{
  if(!expression)
    throw std::invalid_argument("Cannot evaluate a null expression");
  if(expression->op == TreeOperator::LEAF)
  {
    if(expression->leaf_id < 0 ||
       static_cast<std::size_t>(expression->leaf_id) >= leaves.size())
      throw std::out_of_range("Halfspace leaf id is out of range");
    const HalfspaceLeaf& leaf = leaves[expression->leaf_id];
    return leaf.normal.x * point.x + leaf.normal.y * point.y +
           leaf.normal.z * point.z - leaf.offset;
  }

  double result = expression->op == TreeOperator::MIN
                      ? std::numeric_limits<double>::infinity()
                      : -std::numeric_limits<double>::infinity();
  for(const ExpressionPtr& child : expression->children)
  {
    const double value = evaluate(child, leaves, point);
    result = expression->op == TreeOperator::MIN ? std::min(result, value)
                                                 : std::max(result, value);
  }
  return result;
}

inline bool evaluateInside(const ExpressionPtr& expression,
                           const std::vector<HalfspaceLeaf>& leaves,
                           const Vec3& point,
                           double tolerance = 0.0)
{
  return evaluate(expression, leaves, point) <= tolerance;
}

inline std::string structuralKey(const ExpressionPtr& expression)
{
  if(expression->op == TreeOperator::LEAF)
    return "L" + std::to_string(expression->leaf_id);
  std::vector<std::string> child_keys;
  child_keys.reserve(expression->children.size());
  for(const ExpressionPtr& child : expression->children)
    child_keys.push_back(structuralKey(child));
  std::sort(child_keys.begin(), child_keys.end());
  std::ostringstream out;
  out << (expression->op == TreeOperator::MIN ? "m(" : "M(");
  for(const std::string& key : child_keys)
    out << key.size() << ':' << key;
  out << ')';
  return out.str();
}

inline ExpressionPtr simplifyOnce(const ExpressionPtr& expression)
{
  if(!expression || expression->op == TreeOperator::LEAF)
    return expression;

  std::vector<ExpressionPtr> children;
  for(const ExpressionPtr& original_child : expression->children)
  {
    const ExpressionPtr child = simplifyOnce(original_child);
    if(child->op == expression->op)
      children.insert(children.end(), child->children.begin(),
                      child->children.end());
    else
      children.push_back(child);
  }

  std::vector<std::pair<std::string, ExpressionPtr>> keyed;
  keyed.reserve(children.size());
  for(const ExpressionPtr& child : children)
    keyed.emplace_back(structuralKey(child), child);
  std::sort(keyed.begin(), keyed.end(),
            [](const auto& lhs, const auto& rhs) {
              return lhs.first < rhs.first;
            });
  children.clear();
  std::string previous;
  bool have_previous = false;
  for(auto& entry : keyed)
  {
    if(!have_previous || entry.first != previous)
    {
      previous = entry.first;
      have_previous = true;
      children.push_back(std::move(entry.second));
    }
  }

  if(children.size() == 1)
    return children.front();
  return makeNode(expression->op, std::move(children));
}

inline ExpressionPtr simplifyToFixedPoint(const ExpressionPtr& expression)
{
  if(!expression)
    throw std::invalid_argument("Cannot simplify a null expression");
  ExpressionPtr current = expression;
  for(;;)
  {
    ExpressionPtr next = simplifyOnce(current);
    if(structuralKey(next) == structuralKey(current))
      return next;
    current = std::move(next);
  }
}

inline std::size_t treeDepth(const ExpressionPtr& expression)
{
  if(!expression)
    return 0;
  std::size_t depth = 0;
  for(const ExpressionPtr& child : expression->children)
    depth = std::max(depth, treeDepth(child));
  return depth + 1;
}

inline std::size_t treeNodeReferenceCount(const ExpressionPtr& expression)
{
  if(!expression)
    return 0;
  std::size_t count = 1;
  for(const ExpressionPtr& child : expression->children)
    count += treeNodeReferenceCount(child);
  return count;
}

inline SerializedTree serializeTree(const ExpressionPtr& expression)
{
  if(!expression)
    throw std::invalid_argument("Cannot serialize a null expression");
  SerializedTree result;
  std::unordered_map<const Expression*, int> ids;
  const auto visit = [&](const auto& self, const ExpressionPtr& node) -> int {
    const auto found = ids.find(node.get());
    if(found != ids.end())
      return found->second;
    SerializedNode serialized;
    serialized.op = node->op;
    serialized.leaf_id = node->leaf_id;
    serialized.children.reserve(node->children.size());
    for(const ExpressionPtr& child : node->children)
      serialized.children.push_back(self(self, child));
    serialized.id = static_cast<int>(result.nodes.size());
    result.nodes.push_back(std::move(serialized));
    ids.emplace(node.get(), result.nodes.back().id);
    return result.nodes.back().id;
  };
  result.root_id = visit(visit, expression);
  return result;
}

inline std::string indentation(int spaces)
{
  return std::string(static_cast<std::size_t>(std::max(0, spaces)), ' ');
}

inline std::string booleanExpression(const ExpressionPtr& expression,
                                     int indent)
{
  if(expression->op == TreeOperator::LEAF)
    return indentation(indent) + "(psi_" +
           std::to_string(expression->leaf_id) + "(x) <= 0)";
  const std::string joiner = expression->op == TreeOperator::MIN
                                 ? " ||\n"
                                 : " &&\n";
  std::ostringstream out;
  out << indentation(indent) << "(\n";
  for(std::size_t i = 0; i < expression->children.size(); ++i)
  {
    if(i != 0)
      out << joiner;
    out << booleanExpression(expression->children[i], indent + 2);
  }
  out << "\n" << indentation(indent) << ')';
  return out.str();
}

inline std::string scalarExpression(const ExpressionPtr& expression,
                                    int indent)
{
  if(expression->op == TreeOperator::LEAF)
    return indentation(indent) + "psi_" +
           std::to_string(expression->leaf_id) + "(x)";
  std::ostringstream out;
  out << indentation(indent) << operatorName(expression->op) << "(\n";
  for(std::size_t i = 0; i < expression->children.size(); ++i)
  {
    if(i != 0)
      out << ",\n";
    out << scalarExpression(expression->children[i], indent + 2);
  }
  out << "\n" << indentation(indent) << ')';
  return out.str();
}

}  // namespace rokae_demo::halfspace
