// Read-only benchmark adapter: reconstruct the exported tree, then call the
// EXISTING evaluator. Parsing, reference comparisons and output are not timed.
#include <rokae_demo/halfspace_logic_tree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>

namespace HS = rokae_demo::halfspace;
namespace PT = boost::property_tree;
using Clock = std::chrono::steady_clock;

void require(bool value, const char* message)
{
  if(!value) throw std::runtime_error(message);
}

struct Tree
{
  std::vector<HS::HalfspaceLeaf> leaves;
  HS::ExpressionPtr root;
};

HS::Vec3 vector3(const PT::ptree& data)
{
  require(data.size() == 3, "Expected a three-component vector");
  double v[3]{};
  unsigned i = 0;
  for(const auto& item : data)
  {
    v[i] = item.second.get_value<double>();
    require(std::isfinite(v[i]), "Nonfinite vector component");
    ++i;
  }
  return {v[0], v[1], v[2]};
}

Tree readTree(const std::string& path)
{
  PT::ptree json;
  PT::read_json(path, json);
  require(json.get<std::string>("length_unit") == "m" &&
          json.get<std::string>("halfspace_offset_unit") == "m", "Tree units must be meters");
  Tree result;
  result.leaves.resize(json.get_child("leaves").size());
  std::set<int> leaf_ids;
  for(const auto& item : json.get_child("leaves"))
  {
    const auto& leaf = item.second;
    const int id = leaf.get<int>("id");
    require(id >= 0 && static_cast<std::size_t>(id) < result.leaves.size() && leaf_ids.insert(id).second,
            "Leaf IDs must be unique and contiguous");
    auto& p = result.leaves[id];
    p.id = id;
    p.normal = vector3(leaf.get_child("normal"));
    p.point = vector3(leaf.get_child("point"));
    p.offset = leaf.get<double>("offset");
    require(std::isfinite(p.offset), "Nonfinite plane offset");
  }
  std::map<int, const PT::ptree*> nodes;
  for(const auto& item : json.get_child("nodes"))
  {
    const int id = item.second.get<int>("id");
    require(id >= 0 && nodes.emplace(id, &item.second).second, "Invalid or duplicate node ID");
  }
  std::set<int> pending;
  std::map<int, HS::ExpressionPtr> built;
  const auto visit = [&](const auto& self, int id, unsigned depth) -> HS::ExpressionPtr {
    require(depth <= 128, "Tree depth exceeds benchmark safety limit");
    if(built.count(id)) return built.at(id);
    require(nodes.count(id) && pending.insert(id).second, "Missing reference or cyclic tree");
    const auto& n = *nodes.at(id);
    const auto type = n.get<std::string>("type");
    HS::ExpressionPtr expression;
    if(type == "leaf")
    {
      const int leaf = n.get<int>("leaf_id");
      require(leaf_ids.count(leaf), "Missing leaf reference");
      require(n.get_child("children").empty(), "Leaf cannot have children");
      expression = HS::makeLeaf(leaf);
    }
    else
    {
      require(type == "min" || type == "max", "Unknown tree operator");
      std::vector<HS::ExpressionPtr> children;
      for(const auto& child : n.get_child("children")) children.push_back(self(self, child.second.template get_value<int>(), depth+1));
      expression = HS::makeNode(type == "min" ? HS::TreeOperator::MIN : HS::TreeOperator::MAX, std::move(children));
    }
    pending.erase(id);
    built.emplace(id, expression);
    return expression;
  };
  result.root = visit(visit, json.get<int>("root_id"), 0);
  require(built.size() == nodes.size(), "Unreachable nodes in exported tree");
  return result;
}

std::vector<HS::Vec3> readPoints(const std::string& path)
{
  std::ifstream in(path);
  require(bool(in), "Cannot read query XYZ");
  std::vector<HS::Vec3> result;
  std::string line;
  while(std::getline(in, line))
  {
    const auto comment = line.find('#');
    if(comment != std::string::npos) line.resize(comment);
    std::istringstream row(line);
    row >> std::ws;
    if(row.eof()) continue;
    HS::Vec3 p;
    require(bool(row >> p.x >> p.y >> p.z), "Malformed query XYZ");
    row >> std::ws;
    require(row.eof() && std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z), "Nonfinite or extra query fields");
    result.push_back(p);
  }
  require(!in.bad() && !result.empty(), "Empty or unreadable query file");
  return result;
}

std::size_t positive(const std::string& value)
{
  require(!value.empty() && value.find_first_not_of("0123456789") == std::string::npos, "Expected a positive integer");
  const auto n = std::stoull(value);
  require(n > 0 && n <= std::numeric_limits<std::size_t>::max(), "Integer is zero or too large");
  return static_cast<std::size_t>(n);
}

int main(int argc, char** argv)
{
  try
  {
    require(argc >= 4, "Usage: halfspace_tree_benchmark tree.json queries.xyz result.json [--reference raw.json] [--timed-count 256] [--rounds 5] [--repeats 7]");
    std::size_t count = 256, rounds = 5, repeats = 7;
    std::string reference;
    std::set<std::string> seen;
    for(int i = 4; i < argc; i += 2)
    {
      const std::string option = argv[i];
      require(i+1 < argc && seen.insert(option).second, "Missing or repeated benchmark option");
      if(option == "--reference") reference = argv[i+1];
      else if(option == "--timed-count") count = positive(argv[i+1]);
      else if(option == "--rounds") rounds = positive(argv[i+1]);
      else if(option == "--repeats") repeats = positive(argv[i+1]);
      else throw std::runtime_error("Unsupported benchmark option");
    }
    const auto tree = readTree(argv[1]);
    const auto queries = readPoints(argv[2]);
    count = std::min(count, queries.size());
    require(rounds <= std::numeric_limits<std::size_t>::max()/count, "Query count overflow");
    std::vector<double> values;
    for(const auto& p : queries)
    {
      const auto value = HS::evaluate(tree.root, tree.leaves, p);
      require(std::isfinite(value), "Nonfinite tree evaluation");
      values.push_back(value);
    }
    std::size_t decision_mismatches = 0, scalar_mismatches = 0;
    if(!reference.empty())
    {
      const auto old = readTree(reference);
      for(std::size_t i = 0; i < queries.size(); ++i)
      {
        const auto value = HS::evaluate(old.root, old.leaves, queries[i]);
        require(std::isfinite(value), "Nonfinite reference evaluation");
        decision_mismatches += (value <= 0) != (values[i] <= 0);
        scalar_mismatches += value != values[i];
      }
    }
    // Full evaluation above also warms code/data. No short circuit, batching,
    // new tree representation or CBF controller is being benchmarked here.
    volatile double sink = 0;
    std::vector<double> timings;
    for(std::size_t repeat = 0; repeat < repeats; ++repeat)
    {
      double checksum = 0;
      const auto start = Clock::now();
      for(std::size_t r = 0; r < rounds; ++r)
        for(std::size_t i = 0; i < count; ++i) checksum += HS::evaluate(tree.root, tree.leaves, queries[i]);
      sink = checksum;
      timings.push_back(std::chrono::duration<double, std::micro>(Clock::now()-start).count() / (count*rounds));
      require(std::isfinite(checksum), "Benchmark checksum overflow");
    }
    auto sorted = timings;
    std::sort(sorted.begin(), sorted.end());
    const auto median = sorted.size()%2 ? sorted[sorted.size()/2] : (sorted[sorted.size()/2-1]+sorted[sorted.size()/2])/2;
    std::ofstream out(argv[3]);
    require(bool(out), "Cannot write query benchmark result");
    out.exceptions(std::ios::badbit | std::ios::failbit);
    out << std::setprecision(17) << "{\"evaluator\":\"existing_recursive_halfspace_evaluate\",\"length_unit\":\"m\","
        << "\"query_count\":" << queries.size() << ",\"timed_count\":" << count << ",\"rounds\":" << rounds
        << ",\"repeats\":" << repeats << ",\"reference_checked\":" << (reference.empty() ? "false" : "true")
        << ",\"reference_decision_mismatches\":" << decision_mismatches
        << ",\"reference_scalar_mismatches\":" << scalar_mismatches << ",\"checksum\":" << sink
        << ",\"median_us_per_query\":" << median << ",\"min_us_per_query\":" << sorted.front()
        << ",\"max_us_per_query\":" << sorted.back() << ",\"values\":[";
    for(std::size_t i = 0; i < values.size(); ++i) out << (i ? "," : "") << values[i];
    out << "],\"timings_us_per_query\":[";
    for(std::size_t i = 0; i < timings.size(); ++i) out << (i ? "," : "") << timings[i];
    out << "]}\n";
    out.close();
    return decision_mismatches || scalar_mismatches ? 2 : 0;
  }
  catch(const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
