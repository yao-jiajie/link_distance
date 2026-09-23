#include <rokae_demo/boundary_direct_tree.hpp>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <optional>
#include <tuple>
#include <unordered_map>
#ifdef ROKAE_USE_ABSEIL_DIRECT_HASH
#include <absl/container/flat_hash_map.h>
#endif

namespace rokae_demo
{
namespace
{
using Clock=std::chrono::steady_clock;
using Index=std::array<std::size_t,3>;
double ms(Clock::time_point t) { return std::chrono::duration<double,std::milli>(Clock::now()-t).count(); }
void require(bool ok,const char* why) { if(!ok) throw std::runtime_error(why); }
std::size_t index(const Index& p,const Index& d) { return (p[0]*d[1]+p[1])*d[2]+p[2]; }
std::size_t product(const Index& d,std::size_t cap)
{
  std::size_t n=1;
  for(const auto w:d) { require(w && w<=cap/n,"Direct canonical/prefix cell limit"); n*=w; }
  return n;
}
std::vector<unsigned char> canonical(const SparseVoxelOccupancy& o,const BoundaryPlaneRegistry& registry,
                                     Index& dimensions,std::size_t cap)
{
  for(unsigned a=0;a<3;++a) dimensions[a]=registry.coordinates[a].size()+1;
  const auto n=product(dimensions,cap);
  struct DenseAxisLookup
  {
    std::int64_t lower=0;
    std::vector<std::size_t> positions;
  };
  std::array<DenseAxisLookup,3> lookups;
  const auto occupancy_scaled=o.occupied_voxels.size()<=std::numeric_limits<std::size_t>::max()/4
    ? o.occupied_voxels.size()*4 : std::numeric_limits<std::size_t>::max();
  const auto dense_limit=std::max(std::size_t{4096},occupancy_scaled);
  for(unsigned a=0;a<3;++a)
  {
    const auto& coordinates=registry.coordinates[a];
    require(coordinates.size()>=2,"Incomplete direct canonical coordinates");
    const auto span=static_cast<std::uint64_t>(coordinates.back())-
                    static_cast<std::uint64_t>(coordinates.front());
    if(span>dense_limit) continue;
    auto& lookup=lookups[a]; lookup.lower=coordinates.front();
    lookup.positions.resize(static_cast<std::size_t>(span));
    std::size_t position=0;
    for(std::size_t offset=0;offset<lookup.positions.size();++offset)
    {
      const auto coordinate=lookup.lower+static_cast<std::int64_t>(offset);
      while(position<coordinates.size() && coordinates[position]<=coordinate) ++position;
      lookup.positions[offset]=position;
    }
  }
  std::vector<std::size_t> counts(n,0); std::vector<unsigned char> occupied(n,0);
  std::vector<std::size_t> touched; touched.reserve(o.occupied_voxels.size());
  for(const auto& key:o.occupied_voxels)
  {
    Index p{};
    for(unsigned a=0;a<3;++a)
    {
      const auto& c=registry.coordinates[a];
      const auto& lookup=lookups[a];
      if(!lookup.positions.empty())
      {
        const auto offset=static_cast<std::uint64_t>(key[a])-static_cast<std::uint64_t>(lookup.lower);
        require(offset<lookup.positions.size(),"Voxel outside dense direct canonical lookup");
        p[a]=lookup.positions[static_cast<std::size_t>(offset)];
      }
      else p[a]=static_cast<std::size_t>(std::upper_bound(c.begin(),c.end(),key[a])-c.begin());
      require(p[a]>0 && p[a]<c.size(),"Voxel outside direct canonical interior");
    }
    const auto id=index(p,dimensions);
    if(!counts[id]) touched.push_back(id);
    ++counts[id];
  }
  const auto stride_x=dimensions[1]*dimensions[2],stride_y=dimensions[2];
  for(const auto id:touched)
  {
    const auto x=id/stride_x,rest=id%stride_x,y=rest/stride_y,z=rest%stride_y;
    const Index p{x,y,z}; const auto count=counts[id];
    std::size_t volume=1;
    for(unsigned a=0;a<3;++a)
    {
      const auto& c=registry.coordinates[a];
      const auto width=static_cast<std::uint64_t>(c[p[a]])-static_cast<std::uint64_t>(c[p[a]-1]);
      require(width && width<=count/volume,"Mixed compressed canonical cell"); volume*=static_cast<std::size_t>(width);
    }
    require(volume==count,"Canonical occupied volume mismatch"); occupied[id]=1;
  }
  return occupied;
}
std::array<std::size_t,3> supportOffsets(const BoundaryPlaneRegistry& registry)
{
  std::array<std::size_t,3> offsets{}; std::size_t id=0;
  for(unsigned axis=0;axis<3;++axis)
  {
    offsets[axis]=id;
    for(const auto coordinate:registry.coordinates[axis])
    {
      require(id<registry.supports.size(),"Incomplete direct support registry");
      const auto& support=registry.supports[id++];
      require(support.axis==axis && support.coordinate==coordinate,"Unordered direct support registry");
    }
  }
  require(id==registry.supports.size(),"Oversized direct support registry");
  return offsets;
}
struct Prefix
{
  Index dimensions{}; std::size_t stride_x=0,stride_y=0;
  std::size_t maximum_count=0;
  std::vector<std::uint32_t> compact_values;
  std::vector<std::uint64_t> wide_values;
  template<class Value> void fill(std::vector<Value>& values,const BoundaryDirectResult& r,
                                  std::size_t padded_cells)
  {
    values.resize(padded_cells,0);
    const auto source_stride_y=r.dimensions[2];
    const auto source_stride_x=r.dimensions[1]*source_stride_y;
    for(std::size_t x=1;x<dimensions[0];++x) for(std::size_t y=1;y<dimensions[1];++y)
      {
        auto id=x*stride_x+y*stride_y+1;
        auto source=(x-1)*source_stride_x+(y-1)*source_stride_y;
        Value row=0;
        for(std::size_t z=1;z<dimensions[2];++z)
        { row+=static_cast<Value>(r.occupied[source++]); values[id++]=row; }
      }
    for(std::size_t x=1;x<dimensions[0];++x)
    {
      for(std::size_t y=1;y<dimensions[1];++y)
        for(std::size_t z=1;z<dimensions[2];++z)
        { const auto id=x*stride_x+y*stride_y+z; values[id]+=values[id-stride_y]; }
      for(std::size_t y=1;y<dimensions[1];++y)
        for(std::size_t z=1;z<dimensions[2];++z)
        { const auto id=x*stride_x+y*stride_y+z; values[id]+=values[id-stride_x]; }
    }
  }
  template<class Value> std::size_t sumValues(const std::vector<Value>& values,
                                               const Index& lo,const Index& hi) const
  {
    const auto x0=lo[0]*stride_x,x1=hi[0]*stride_x;
    const auto y0=lo[1]*stride_y,y1=hi[1]*stride_y;
    Value count=values[x1+y1+hi[2]];
    count-=values[x0+y1+hi[2]];
    count-=values[x1+y0+hi[2]];
    count+=values[x0+y0+hi[2]];
    count-=values[x1+y1+lo[2]];
    count+=values[x0+y1+lo[2]];
    count+=values[x1+y0+lo[2]];
    count-=values[x0+y0+lo[2]];
    require(count<=maximum_count,"Invalid canonical prefix count");
    return static_cast<std::size_t>(count);
  }
  Prefix(const BoundaryDirectResult& r)
  {
    for(unsigned a=0;a<3;++a) dimensions[a]=r.dimensions[a]+1;
    stride_y=dimensions[2]; stride_x=dimensions[1]*stride_y;
    // With each dimension >=3, padded product < 3 * canonical product.
    maximum_count=r.occupied.size();
    const auto padded_cells=product(dimensions,r.occupied.size()*3);
    if(maximum_count<=std::numeric_limits<std::uint32_t>::max())
      fill(compact_values,r,padded_cells);
    else fill(wide_values,r,padded_cells);
  }
  std::size_t sum(const Index& lo,const Index& hi) const
  {
    return compact_values.empty() ? sumValues(wide_values,lo,hi) : sumValues(compact_values,lo,hi);
  }
};
std::vector<HS::ExpressionPtr> materializeExpressions(const HS::SerializedTree& tree)
{
  require(!tree.nodes.empty() && tree.root_id>=0 &&
          static_cast<std::size_t>(tree.root_id)<tree.nodes.size(),"Invalid serialized direct tree");
  std::vector<HS::ExpressionPtr> expressions(tree.nodes.size());
  for(std::size_t i=0;i<tree.nodes.size();++i)
  {
    const auto& node=tree.nodes[i];
    require(node.id==static_cast<int>(i),"Invalid serialized direct node order");
    if(node.op==HS::TreeOperator::LEAF)
      expressions[i]=HS::makeLeaf(node.leaf_id);
    else
    {
      std::vector<HS::ExpressionPtr> children; children.reserve(node.children.size());
      for(const auto child:node.children)
      {
        require(child>=0 && child<node.id,"Invalid serialized direct child");
        children.push_back(expressions[static_cast<std::size_t>(child)]);
      }
      expressions[i]=HS::makeNode(node.op,std::move(children));
    }
  }
  return expressions;
}
HS::ExpressionPtr materializeExpression(const HS::SerializedTree& tree)
{ auto expressions=materializeExpressions(tree); return expressions[static_cast<std::size_t>(tree.root_id)]; }
// Constants are construction tags, never final leaves or numerical infinities.
struct Term
{
  enum Kind { FREE, OCCUPIED, EXPRESSION } kind=FREE;
  HS::ExpressionPtr expression;
  std::size_t nodes=0, planes=0;
  unsigned depth=0;
  std::size_t identity=0; // construction identity, independent of structural interning
  int serialized_id=-1;
};
struct InternKey
{
  HS::TreeOperator op=HS::TreeOperator::MIN;
  unsigned char count=0;
  std::array<const void*,3> children{};
  bool operator==(const InternKey& other) const
  { return op==other.op && count==other.count && children==other.children; }
};
struct InternKeyHash
{
  std::size_t operator()(const InternKey& key) const
  {
    std::size_t hash=static_cast<std::size_t>(key.op);
    const auto mix=[&](std::size_t value) {
      hash^=value+std::size_t{0x9e3779b9}+(hash<<6)+(hash>>2);
    };
    mix(key.count);
    for(unsigned char i=0;i<key.count;++i)
      mix(reinterpret_cast<std::size_t>(key.children[i]));
    return hash;
  }
};
struct InternedNode
{
  HS::ExpressionPtr expression;
  int serialized_id=-1;
};
struct SerializedInternKey
{
  HS::TreeOperator op;
  unsigned char count;
  std::array<int,3> children;
  bool operator==(const SerializedInternKey& other) const
  {
    return op==other.op && count==other.count &&
      std::equal(children.begin(),children.begin()+count,other.children.begin());
  }
};
struct SerializedInternKeyHash
{
  std::size_t operator()(const SerializedInternKey& key) const
  {
    std::size_t hash=static_cast<std::size_t>(key.op);
    const auto mix=[&](std::size_t value) {
      hash^=value+std::size_t{0x9e3779b9}+(hash<<6)+(hash>>2);
    };
    mix(key.count);
    for(unsigned char i=0;i<key.count;++i) mix(static_cast<std::size_t>(key.children[i]));
    return hash;
  }
};
class SerializedInterner
{
  struct Slot
  {
    SerializedInternKey key;
    int value;
    Slot() noexcept:value(-1) {}
  };
  std::vector<Slot> slots;
  std::size_t count=0;

  static std::size_t capacityFor(std::size_t requested)
  {
    std::size_t capacity=16;
    while(capacity<requested)
    {
      require(capacity<=std::numeric_limits<std::size_t>::max()/2,
              "Direct serialized interner capacity overflow");
      capacity*=2;
    }
    return capacity;
  }
  void rehash(std::size_t capacity)
  {
    auto old=std::move(slots);
    slots.clear();
    slots.resize(capacityFor(capacity));
    count=0;
    for(const auto& slot:old) if(slot.value>=0) insert(slot.key,slot.value);
  }
public:
  void reserve(std::size_t expected)
  {
    const auto requested=expected<=std::numeric_limits<std::size_t>::max()/2
      ? expected*2 : std::numeric_limits<std::size_t>::max();
    if(slots.size()<requested) rehash(requested);
  }
  template<class Create> int findOrCreate(const SerializedInternKey& key,Create&& create)
  {
    if(slots.empty()) rehash(16);
    if((count+1)*10>=slots.size()*7) rehash(slots.size()*2);
    const auto mask=slots.size()-1;
    auto position=SerializedInternKeyHash{}(key)&mask;
    while(slots[position].value>=0)
    {
      if(slots[position].key==key) return slots[position].value;
      position=(position+1)&mask;
    }
    const auto value=create();
    require(value>=0,"Invalid serialized intern value");
    slots[position].key=key; slots[position].value=value; ++count;
    return value;
  }
  void insert(const SerializedInternKey& key,int value)
  {
    require(value>=0,"Invalid serialized intern value");
    if(slots.empty()) rehash(16);
    if((count+1)*10>=slots.size()*7) rehash(slots.size()*2);
    const auto mask=slots.size()-1;
    auto position=SerializedInternKeyHash{}(key)&mask;
    while(slots[position].value>=0)
    {
      if(slots[position].key==key)
      {
        require(slots[position].value==value,"Conflicting serialized intern value");
        return;
      }
      position=(position+1)&mask;
    }
    slots[position].key=key; slots[position].value=value; ++count;
  }
};
class Builder
{
  BoundaryDirectResult& r; const BoundaryDirectOptions& limits; const SparseVoxelOccupancy& occupancy;
  Prefix prefix;
  std::array<std::size_t,3> support_offsets{};
  std::vector<std::array<std::optional<Term>,2>> leaves;
#ifdef ROKAE_USE_ABSEIL_DIRECT_HASH
  absl::flat_hash_map<InternKey,InternedNode,InternKeyHash> interned_nodes;
#else
  std::unordered_map<InternKey,InternedNode,InternKeyHash> interned_nodes;
#endif
  SerializedInterner serialized_interned_nodes;
  std::size_t next_identity=1;
  void allocation() { require(r.allocated_nodes<limits.max_nodes,"Direct expression allocation limit"); ++r.allocated_nodes; }
  Term leaf(std::size_t support,int side)
  {
    require(support<leaves.size() && (side==-1 || side==1),"Invalid direct leaf key");
    auto& cached=leaves[support][side>0]; if(cached) return *cached;
    const auto& source=r.registry.supports[support]; const double coordinate=static_cast<double>(source.coordinate)*occupancy.voxel_size;
    Plane plane; double normal[3]{},point[3]{}; normal[source.axis]=side; point[source.axis]=coordinate;
    plane.normal={normal[0],normal[1],normal[2]}; plane.point={point[0],point[1],point[2]}; plane.offset=side*coordinate;
    const auto id=findOrAddLeaf(plane,-1,0,r.tree.leaves,"boundary_direct_support",true);
    require(static_cast<std::size_t>(id)==r.leaf_sources.size(),"Unexpected direct leaf alias");
    r.leaf_sources.emplace_back(support,side); allocation();
    HS::SerializedNode serialized;
    serialized.id=static_cast<int>(r.serialized.nodes.size());
    serialized.op=HS::TreeOperator::LEAF; serialized.leaf_id=id;
    r.serialized.nodes.push_back(std::move(serialized));
    Term result{Term::EXPRESSION,limits.materialize_expression ? HS::makeLeaf(id) : HS::ExpressionPtr{},1,1,1,next_identity++,
                static_cast<int>(r.serialized.nodes.size()-1)};
    cached=result; return result;
  }
  Term node(HS::TreeOperator op,std::initializer_list<Term> terms)
  {
    Term out{Term::EXPRESSION,{},1,0,1,next_identity++,-1};
    std::vector<HS::ExpressionPtr> children;
    if(limits.materialize_expression) children.reserve(terms.size());
    std::array<int,3> serialized_children{}; unsigned char child_count=0;
    for(const auto& t:terms)
    {
      require(t.kind==Term::EXPRESSION,"Constant escaped Phi folding");
      require(t.nodes<=limits.max_expanded_nodes-out.nodes,"Direct expanded expression limit");
      require(t.serialized_id>=0,"Direct child missing serialized ID");
      require(child_count<serialized_children.size(),"Direct serialized child arity limit");
      out.nodes+=t.nodes; out.planes+=t.planes; out.depth=std::max(out.depth,t.depth+1);
      if(limits.materialize_expression) children.push_back(t.expression);
      serialized_children[child_count++]=t.serialized_id;
    }
    require(out.depth<=limits.max_depth,"Direct expression depth limit");
    allocation();
    InternKey key; key.op=op;
    SerializedInternKey serialized_key{}; serialized_key.op=op;
    if(limits.intern_subexpressions && limits.materialize_expression)
    {
        require(children.size()<=key.children.size(),"Direct intern key arity limit");
        key.count=static_cast<unsigned char>(children.size());
        for(std::size_t i=0;i<children.size();++i) key.children[i]=children[i].get();
        const auto found=interned_nodes.find(key);
        if(found!=interned_nodes.end())
        {
          out.expression=found->second.expression;
          out.serialized_id=found->second.serialized_id;
          return out;
        }
    }
    else if(limits.intern_subexpressions)
    {
      serialized_key.count=child_count;
      for(unsigned char i=0;i<child_count;++i) serialized_key.children[i]=serialized_children[i];
      out.serialized_id=serialized_interned_nodes.findOrCreate(serialized_key,[&] {
        HS::SerializedNode serialized;
        serialized.id=static_cast<int>(r.serialized.nodes.size());
        serialized.op=op;
        serialized.children.assign(serialized_children.begin(),serialized_children.begin()+child_count);
        r.serialized.nodes.push_back(std::move(serialized));
        return static_cast<int>(r.serialized.nodes.size()-1);
      });
      return out;
    }
    if(limits.materialize_expression) out.expression=HS::makeNode(op,std::move(children));
    HS::SerializedNode serialized;
    serialized.id=static_cast<int>(r.serialized.nodes.size());
    serialized.op=op;
    serialized.children.assign(serialized_children.begin(),serialized_children.begin()+child_count);
    r.serialized.nodes.push_back(std::move(serialized));
    out.serialized_id=static_cast<int>(r.serialized.nodes.size()-1);
    if(limits.intern_subexpressions && limits.materialize_expression)
      interned_nodes.emplace(key,InternedNode{out.expression,out.serialized_id});
    return out;
  }
  Term phi(std::size_t support,const Term& a,const Term& b)
  {
    if(a.kind==b.kind && (a.kind!=Term::EXPRESSION || a.identity==b.identity)) { ++r.phi_folded; return a; }
    if(a.kind!=Term::EXPRESSION || b.kind!=Term::EXPRESSION)
    {
      ++r.phi_folded;
      if(a.kind==Term::FREE && b.kind==Term::OCCUPIED) return leaf(support,1);
      if(a.kind==Term::OCCUPIED && b.kind==Term::FREE) return leaf(support,-1);
      if(a.kind==Term::FREE) return node(HS::TreeOperator::MIN,{leaf(support,1),b});
      if(a.kind==Term::OCCUPIED) return node(HS::TreeOperator::MAX,{b,leaf(support,-1)});
      if(b.kind==Term::FREE) return node(HS::TreeOperator::MIN,{a,leaf(support,-1)});
      return node(HS::TreeOperator::MAX,{a,leaf(support,1)});
    }
    ++r.phi_full;
    const auto s=leaf(support,1),neg_s=leaf(support,-1);
    const auto ab=node(HS::TreeOperator::MIN,{a,b}), an=node(HS::TreeOperator::MIN,{a,neg_s}), sb=node(HS::TreeOperator::MIN,{s,b});
    return node(HS::TreeOperator::MAX,{ab,an,sb});
  }
public:
  Builder(BoundaryDirectResult& result,const BoundaryDirectOptions& opts,const SparseVoxelOccupancy& o)
    : r(result),limits(opts),occupancy(o),prefix(result),
      support_offsets(supportOffsets(result.registry)),leaves(result.registry.supports.size())
  {
    if(limits.intern_subexpressions)
    {
      if(limits.materialize_expression)
        interned_nodes.reserve(std::min(limits.max_nodes,std::size_t{16384}));
      else serialized_interned_nodes.reserve(std::min(limits.max_nodes,std::size_t{16384}));
    }
    r.serialized.nodes.reserve(std::min(limits.max_nodes,std::size_t{16384}));
    if(limits.retain_partition_diagnostics)
      r.splits.reserve(std::min(limits.max_nodes,std::size_t{16384}));
  }
  void finishSerialized(int root)
  {
    // The builder can create leaves/subexpressions that Phi later folds away.
    // Compact the emitted IDs in the same child-first order as serializeTree,
    // without revisiting ExpressionPtr or hashing its pointers.
    auto emitted=std::move(r.serialized.nodes);
    std::vector<int> remap(emitted.size(),-1);
    r.serialized.nodes.reserve(emitted.size());
    const auto visit=[&](const auto& self,int source)->int {
      require(source>=0 && static_cast<std::size_t>(source)<emitted.size(),
              "Invalid emitted direct node");
      auto& mapped=remap[source];
      if(mapped>=0) return mapped;
      auto& node=emitted[source];
      for(auto& child:node.children) child=self(self,child);
      mapped=static_cast<int>(r.serialized.nodes.size());
      node.id=mapped;
      r.serialized.nodes.push_back(std::move(node));
      return mapped;
    };
    r.serialized.root_id=visit(visit,root);
  }
  Term build(const Index& lo,const Index& hi,unsigned depth)
  {
    require(depth<=limits.max_depth,"Direct partition depth limit"); r.recursion_depth=std::max(r.recursion_depth,depth);
    const auto n=prefix.sum(lo,hi); std::size_t volume=1; for(unsigned a=0;a<3;++a) volume*=hi[a]-lo[a];
    if(!n)
    {
      ++r.free_terminals; return Term{};
    }
    if(n==volume)
    {
      ++r.occupied_terminals; return Term{Term::OCCUPIED,{},0,0,0,0,-1};
    }
    using Score=std::tuple<unsigned,std::size_t,unsigned,std::size_t>;
    Score best{3,0,0,0}; unsigned axis=0; std::size_t cut=0,left_count=0;
    for(unsigned a=0;a<3;++a)
    {
      const auto slab=volume/(hi[a]-lo[a]);
      for(std::size_t k=lo[a]+1;k<hi[a];++k)
      {
        require(r.split_checks<limits.max_split_checks,"Direct split candidate limit"); ++r.split_checks;
        auto left=hi; left[a]=k; const auto lc=prefix.sum(lo,left);
        const auto lv=slab*(k-lo[a]),rv=volume-lv,rc=n-lc;
        const unsigned mixed=unsigned(lc>0 && lc<lv)+unsigned(rc>0 && rc<rv);
        const Score score{mixed,std::max(lv,rv),a,k};
        if(score<best) { best=score; axis=a; cut=k; left_count=lc; }
      }
    }
    require(cut>lo[axis] && cut<hi[axis],"Mixed indivisible direct region");
    const auto support=support_offsets[axis]+cut-1;
    if(limits.retain_partition_diagnostics)
      r.splits.push_back({lo,hi,support,cut,left_count,n-left_count});
    auto left_hi=hi,right_lo=lo; left_hi[axis]=right_lo[axis]=cut;
    // Explicit evaluation order keeps leaf IDs deterministic across compilers.
    const auto a=build(lo,left_hi,depth+1);
    const auto b=build(right_lo,hi,depth+1);
    return phi(support,a,b);
  }
};

} // namespace

BoundaryDirectResult buildBoundaryDirectTree(const SparseVoxelOccupancy& o,const OctreeBoundary& boundary,const BoundaryDirectOptions& options)
{
  require(options.max_canonical_cells>0 && options.max_canonical_cells<=std::numeric_limits<std::size_t>::max()/8 &&
          options.max_canonical_cells<=static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()/8) &&
          options.max_nodes>0 && options.max_nodes<=static_cast<std::size_t>(std::numeric_limits<int>::max()) &&
          options.max_expanded_nodes>0 && options.max_split_checks>0 && options.max_depth>0 && options.max_depth<=256,
          "Invalid direct work limits");
  const auto core=Clock::now(); auto start=Clock::now(); BoundaryDirectResult r;
  r.registry=buildBoundaryPlaneRegistry(boundary); r.registry_ms=ms(start); start=Clock::now();
  r.occupied=canonical(o,r.registry,r.dimensions,options.max_canonical_cells); r.canonical_ms=ms(start);
  const auto tree_start=Clock::now(); start=tree_start;
  Builder builder(r,options,o); r.prefix_ms=ms(start); start=Clock::now();
  const auto root=builder.build({0,0,0},r.dimensions,0); r.recursion_ms=ms(start);
  require(root.kind==Term::EXPRESSION,"Unbounded/empty geometry cannot export a constant-only direct tree");
  if(options.materialize_expression)
  { r.tree.raw=root.expression; r.tree.simplified=root.expression; }
  // Otherwise serialized is the primary executable representation and the
  // legacy ExpressionPtr graph is reconstructed only by diagnostics.
  r.tree.leaf_references=r.expanded_plane_refs=root.planes; r.expanded_nodes=root.nodes; r.expression_depth=root.depth;
  start=Clock::now(); builder.finishSerialized(root.serialized_id); r.serialization_ms=ms(start);
  r.construction_interned=options.intern_subexpressions;
  r.constructed_nodes=r.serialized.nodes.size(); r.tree_ms=ms(tree_start); r.core_ms=ms(core); return r;
}

BoundaryDirectResult buildBoundaryDirectTree(const SparseVoxelOccupancy& o,
                                             const BoundaryDirectOptions& options)
{
  require(options.max_canonical_cells>0 && options.max_canonical_cells<=std::numeric_limits<std::size_t>::max()/8 &&
          options.max_canonical_cells<=static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()/8) &&
          options.max_nodes>0 && options.max_nodes<=static_cast<std::size_t>(std::numeric_limits<int>::max()) &&
          options.max_expanded_nodes>0 && options.max_split_checks>0 && options.max_depth>0 && options.max_depth<=256,
          "Invalid direct work limits");
  const auto core=Clock::now(); auto start=Clock::now(); BoundaryDirectResult r;
  r.registry=buildBoundaryPlaneRegistry(o); r.registry_ms=ms(start); start=Clock::now();
  r.occupied=canonical(o,r.registry,r.dimensions,options.max_canonical_cells); r.canonical_ms=ms(start);
  const auto tree_start=Clock::now(); start=tree_start;
  Builder builder(r,options,o); r.prefix_ms=ms(start); start=Clock::now();
  const auto root=builder.build({0,0,0},r.dimensions,0); r.recursion_ms=ms(start);
  require(root.kind==Term::EXPRESSION,"Unbounded/empty geometry cannot export a constant-only direct tree");
  if(options.materialize_expression)
  { r.tree.raw=root.expression; r.tree.simplified=root.expression; }
  // Otherwise serialized is the primary executable representation and the
  // legacy ExpressionPtr graph is reconstructed only by diagnostics.
  r.tree.leaf_references=r.expanded_plane_refs=root.planes; r.expanded_nodes=root.nodes; r.expression_depth=root.depth;
  start=Clock::now(); builder.finishSerialized(root.serialized_id); r.serialization_ms=ms(start);
  r.construction_interned=options.intern_subexpressions;
  r.constructed_nodes=r.serialized.nodes.size(); r.tree_ms=ms(tree_start); r.core_ms=ms(core); return r;
}


BoundarySignDiagnostics validateBoundaryDirectTree(const SparseVoxelOccupancy& o,const OctreeBoundary& b,
                                                   const BoundaryDirectResult& r,const BoundaryDirectOptions& options)
{
  validateVoxelBoundary(o.occupied_voxels,b); const auto registry=buildBoundaryPlaneRegistry(b);
  require(r.registry.coordinates==registry.coordinates && r.registry.supports.size()==registry.supports.size(),"Direct registry mismatch");
  for(std::size_t i=0;i<registry.supports.size();++i)
  {
    const auto& a=r.registry.supports[i]; const auto& expected=registry.supports[i];
    require(a.axis==expected.axis && a.coordinate==expected.coordinate &&
            a.outward_signs==expected.outward_signs,"Direct support geometry mismatch");
    if(r.registry.patch_provenance_complete)
      require(a.patches==expected.patches,"Direct support provenance mismatch");
    else
      require(a.patches.empty(),"Partial direct support provenance");
  }
  Index dimensions; const auto cells=canonical(o,registry,dimensions,options.max_canonical_cells);
  require(dimensions==r.dimensions && cells==r.occupied,"Direct canonical occupancy mismatch");
  auto tree=r.tree;
  if(!tree.raw)
  { tree.raw=materializeExpression(r.serialized); tree.simplified=tree.raw; }
  return diagnoseBoundaryPlaneSigns(o.occupied_voxels,o.voxel_size,registry,r.leaf_sources,tree,options.signs);
}

int boundaryDirectExpectedSign(const BoundaryDirectResult& r,double h,const HS::Vec3& point)
{
  require(std::isfinite(point.x) && std::isfinite(point.y) && std::isfinite(point.z),"Nonfinite direct point");
  const double p[3]{point.x,point.y,point.z}; Index position{}; unsigned on[3]{};
  for(unsigned a=0;a<3;++a)
  {
    const auto& c=r.registry.coordinates[a];
    const auto it=std::lower_bound(c.begin(),c.end(),p[a],[&](auto k,double v) { return static_cast<double>(k)*h<v; });
    position[a]=static_cast<std::size_t>(it-c.begin()); on[a]=it!=c.end() && static_cast<double>(*it)*h==p[a];
  }
  bool any=false,all=true;
  for(unsigned x=0;x<=on[0];++x) for(unsigned y=0;y<=on[1];++y) for(unsigned z=0;z<=on[2];++z)
  {
    const bool filled=r.occupied.at(index({position[0]+x,position[1]+y,position[2]+z},r.dimensions)); any|=filled; all&=filled;
  }
  return all ? 1 : any ? 0 : -1;
}

void writeBoundaryDirectDiagnostics(const std::filesystem::path& path,const BoundaryDirectResult& r,const BoundarySignDiagnostics& d)
{
  std::ofstream out(path); require(bool(out),"Cannot write direct diagnostics"); out.exceptions(std::ios::badbit|std::ios::failbit);
  out << std::boolalpha << std::setprecision(17);
  const auto array=[&](const auto& values) { out << '['; bool first=true; for(const auto v:values) { if(!first) out << ','; first=false; out << v; } out << ']'; };
  out << "{\"method\":\"integer_plane_arrangement_three_sign_propagation\",\"applied\":" << d.applied
      << ",\"strict_sign_certified\":" << d.strict_sign_certified << ",\"closure_set_certified\":" << d.closure_set_certified
      << ",\"reason\":\"" << d.reason << "\",\"strata_count\":" << d.strata_count << ",\"unsafe_free_count\":" << d.unsafe_free_count
      << ",\"missed_free_count\":" << d.missed_free_count << ",\"true_boundary_sign_errors\":" << d.boundary_sign_errors
      << ",\"raw_simplified_sign_errors\":" << d.raw_simplified_sign_errors << ",\"free_zeros_by_dimension\":";
  array(d.free_zeros); out << ",\"occupied_zeros_by_dimension\":"; array(d.occupied_zeros);
  out << ",\"dimension_order\":[\"vertex\",\"edge\",\"face\",\"volume\"],\"coordinates_grid\":[";
  for(unsigned a=0;a<3;++a) { if(a) out << ','; array(r.registry.coordinates[a]); }
  out << "],\"stratum_index_rule\":\"2*k+1 is plane k; even indices are open intervals including unbounded\",\"witnesses\":[";
  for(std::size_t i=0;i<d.witnesses.size();++i)
  { const auto& w=d.witnesses[i]; out << (i ? "," : "") << "{\"strata\":"; array(w.strata); out << ",\"dimension\":" << w.dimension << ",\"expected\":" << w.expected << ",\"actual\":" << w.actual << '}'; }
  out << "]}\n"; out.close();
}
} // namespace rokae_demo
