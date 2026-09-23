#pragma once
#include <rokae_demo/halfspace_batch.hpp>
#include <rokae_demo/halfspace_logic_tree.hpp>

namespace rokae_demo::halfspace
{
struct PreinternedProgramTag { explicit PreinternedProgramTag()=default; };
inline constexpr PreinternedProgramTag preinterned_program{};

// Immutable execution snapshot of the EXISTING serialized DAG. No Boolean
// rewriting, reassociation, new expression operators or geometry changes.
// Compile once per geometry frame; rebuild if tree/coefficients change.
// Share the program between threads, but give EACH worker its own workspace.
class CompiledTree
{
  friend class SmoothTree; // execution-only LSE adapter; never rewrites the DAG
  struct Instruction { TreeOperator op; int leaf,begin,count; };
  struct Coefficients { Vec3 normal; double offset; };
  std::vector<Instruction> instructions_;
  std::vector<int> children_;
  std::vector<Coefficients> planes_;
  int root_=-1;
  bool structural_sharing_=false;
  std::size_t source_nodes_=0;
  void shareIdenticalSubexpressions()
  {
    // Validate the ENTIRE input before interning. Hashes select buckets only;
    // equality checks exact operators, leaf IDs and ordered child occurrences.
    // Compact in place: writes end at/before the current source range.
    std::unordered_multimap<std::size_t,int> buckets; buckets.reserve(source_nodes_);
    std::vector<int> remap(source_nodes_),mapped;
    std::size_t written_nodes=0,written_edges=0;
    for(std::size_t i=0;i<source_nodes_;++i)
    {
      const auto n=instructions_[i]; mapped.clear();
      for(int k=n.begin;k<n.begin+n.count;++k) mapped.push_back(remap[children_[k]]);
      std::size_t hash=static_cast<std::size_t>(n.op);
      const auto mix=[&](std::size_t v) { hash^=v+std::size_t{0x9e3779b9}+(hash<<6)+(hash>>2); };
      if(n.op==TreeOperator::LEAF) mix(static_cast<std::size_t>(n.leaf));
      mix(mapped.size()); for(const auto child:mapped) mix(static_cast<std::size_t>(child));
      int id=-1; const auto range=buckets.equal_range(hash);
      for(auto it=range.first;it!=range.second;++it)
      {
        const auto& old=instructions_[it->second];
        if(old.op==n.op && old.count==n.count && (n.op!=TreeOperator::LEAF || old.leaf==n.leaf) &&
           std::equal(mapped.begin(),mapped.end(),children_.begin()+old.begin))
        { id=it->second; break; }
      }
      if(id<0)
      {
        id=static_cast<int>(written_nodes++);
        instructions_[id]={n.op,n.leaf,static_cast<int>(written_edges),n.count};
        std::copy(mapped.begin(),mapped.end(),children_.begin()+written_edges);
        written_edges+=mapped.size(); buckets.emplace(hash,id);
      }
      remap[i]=id;
    }
    root_=remap[root_]; instructions_.resize(written_nodes); children_.resize(written_edges);
    if(static_cast<std::size_t>(root_)!=written_nodes-1)
      throw std::logic_error("Interned root must remain last in a reachable postorder DAG");
    // Retain capacities rather than paying for extra shrinking copies.
  }
public:
  explicit CompiledTree(const SerializedTree& tree,const std::vector<HalfspaceLeaf>& leaves,bool share_subexpressions=false)
    : structural_sharing_(share_subexpressions)
  {
    const auto limit=static_cast<std::size_t>(std::numeric_limits<int>::max());
    if(tree.nodes.empty() || tree.nodes.size()>limit || tree.root_id<0 ||
       static_cast<std::size_t>(tree.root_id)!=tree.nodes.size()-1 || leaves.size()>limit)
      throw std::invalid_argument("Compiled tree requires a nonempty postorder DAG with last root");
    source_nodes_=tree.nodes.size();
    std::size_t source_edges=0;
    for(const auto& n:tree.nodes)
    {
      if(n.children.size()>limit-source_edges)
        throw std::invalid_argument("Invalid compiled edge count");
      source_edges+=n.children.size();
    }
    root_=tree.root_id; instructions_.reserve(tree.nodes.size());
    children_.reserve(source_edges); planes_.reserve(leaves.size());
    for(const auto& p:leaves)
    {
      if(!std::isfinite(p.normal.x) || !std::isfinite(p.normal.y) || !std::isfinite(p.normal.z) || !std::isfinite(p.offset))
        throw std::invalid_argument("Compiled tree requires finite plane coefficients");
      planes_.push_back({p.normal,p.offset});
    }
    for(std::size_t i=0;i<tree.nodes.size();++i)
    {
      const auto& n=tree.nodes[i];
      if(n.id!=static_cast<int>(i) || n.children.size()>limit-children_.size())
        throw std::invalid_argument("Invalid compiled node ID/edge count");
      if(n.op==TreeOperator::LEAF)
      {
        if(!n.children.empty() || n.leaf_id<0 || static_cast<std::size_t>(n.leaf_id)>=planes_.size())
          throw std::invalid_argument("Invalid compiled leaf");
      }
      else if((n.op!=TreeOperator::MIN && n.op!=TreeOperator::MAX) || n.children.empty())
        throw std::invalid_argument("Invalid compiled MIN/MAX");
      for(const auto child:n.children)
        if(child<0 || child>=n.id) throw std::invalid_argument("Compiled DAG has cycle/forward/invalid child");
      instructions_.push_back({n.op,n.leaf_id,static_cast<int>(children_.size()),static_cast<int>(n.children.size())});
      children_.insert(children_.end(),n.children.begin(),n.children.end());
    }
    std::vector<unsigned char> reachable(instructions_.size(),0); reachable.back()=1;
    for(std::size_t i=reachable.size();i-->0;)
      if(reachable[i]) for(const auto child:tree.nodes[i].children) reachable[child]=1;
    if(std::find(reachable.begin(),reachable.end(),0)!=reachable.end())
      throw std::invalid_argument("Compiled DAG contains unreachable nodes");
    if(share_subexpressions) shareIdenticalSubexpressions();
  }
  explicit CompiledTree(const ExpressionPtr& tree,const std::vector<HalfspaceLeaf>& leaves,bool share_subexpressions=false)
    : CompiledTree(serializeTree(tree),leaves,share_subexpressions) {}
  // The caller guarantees that exact ordered subexpressions were already
  // interned while constructing this postorder DAG. Structural sharing is
  // recorded without paying for a second hash/equality pass.
  explicit CompiledTree(const SerializedTree& tree,const std::vector<HalfspaceLeaf>& leaves,
                        PreinternedProgramTag)
    : CompiledTree(tree,leaves,false)
  { structural_sharing_=true; }

  std::size_t nodeCount() const { return instructions_.size(); }
  std::size_t sourceNodeCount() const { return source_nodes_; }
  std::size_t sharedNodeCount() const { return source_nodes_-nodeCount(); }
  bool structuralSharingEnabled() const { return structural_sharing_; }
  std::size_t edgeCount() const { return children_.size(); }
  std::size_t workPerQuery() const { return nodeCount()+edgeCount(); }
  std::size_t programBytes() const
  { return sizeof(*this)+instructions_.capacity()*sizeof(Instruction)+children_.capacity()*sizeof(int)+planes_.capacity()*sizeof(Coefficients); }
  std::vector<double> makeWorkspace() const { return std::vector<double>(nodeCount()); }
  ValueBatchWorkspace makeBatchWorkspace(std::size_t worker_count) const
  { return ValueBatchWorkspace(worker_count,nodeCount()); }

  // No allocation, hashing, recursive calls or cross-query cache in this loop.
  // Child order, +/- infinity initialization and plane arithmetic deliberately
  // match evaluate(), including signed zero and finite-input overflow behavior.
  double evaluate(const Vec3& p,std::vector<double>& workspace) const
  {
    if(workspace.size()!=nodeCount()) throw std::invalid_argument("Wrong compiled tree workspace size");
    if(!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
      throw std::invalid_argument("Compiled tree query must be finite");
    for(std::size_t i=0;i<instructions_.size();++i)
    {
      const auto& n=instructions_[i];
      if(n.op==TreeOperator::LEAF)
      {
        const auto& leaf=planes_[n.leaf];
        workspace[i]=leaf.normal.x*p.x+leaf.normal.y*p.y+leaf.normal.z*p.z-leaf.offset;
      }
      else
      {
        double value=n.op==TreeOperator::MIN ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity();
        for(int k=n.begin;k<n.begin+n.count;++k)
          value=n.op==TreeOperator::MIN ? std::min(value,workspace[children_[k]]) : std::max(value,workspace[children_[k]]);
        workspace[i]=value;
      }
    }
    return workspace[root_];
  }

  // Output order matches points. Each worker runs the scalar evaluator with a
  // private reusable workspace, so arithmetic and signed-zero behavior within
  // every point are unchanged. A failed point may leave other outputs written.
  void evaluateBatch(const std::vector<Vec3>& points,std::vector<double>& values,
                     ValueBatchWorkspace& workspace) const
  {
    if(values.size()!=points.size()) throw std::invalid_argument("Wrong compiled batch output size");
    if(workspace.workers_.empty() || !workspace.pool_ ||
       workspace.values_per_worker_!=nodeCount())
      throw std::invalid_argument("Wrong compiled batch workspace size");
    detail::evaluatePointBatch(points.size(),*workspace.pool_,[&](std::size_t worker,std::size_t begin,std::size_t end) {
      auto& work=workspace.workers_[worker];
      for(auto i=begin;i<end;++i) values[i]=evaluate(points[i],work);
    });
  }
};
} // namespace rokae_demo::halfspace
