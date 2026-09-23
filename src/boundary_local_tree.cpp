#include <rokae_demo/boundary_local_tree.hpp>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <optional>
#include <tuple>

namespace rokae_demo
{
namespace
{
using Clock = std::chrono::steady_clock;
double ms(Clock::time_point t) { return std::chrono::duration<double,std::milli>(Clock::now()-t).count(); }
void require(bool ok, const char* message) { if(!ok) throw std::runtime_error(message); }
std::size_t indexOf(const std::vector<std::int64_t>& c,std::int64_t p)
{
  const auto i = std::lower_bound(c.begin(),c.end(),p);
  require(i != c.end() && *i == p,"Missing registered coordinate");
  return static_cast<std::size_t>(i-c.begin());
}
int symbolicValue(const HS::SerializedTree& tree,const std::vector<int>& leaves,std::vector<int>& values)
{
  for(const auto& n : tree.nodes)
  {
    if(n.op == HS::TreeOperator::LEAF) values[n.id] = leaves.at(n.leaf_id);
    else
    {
      int value = n.op == HS::TreeOperator::MIN ? 1 : -1;
      for(const auto id : n.children) value = n.op == HS::TreeOperator::MIN ? std::min(value,values[id]) : std::max(value,values[id]);
      values[n.id] = value;
    }
  }
  return values.at(tree.root_id);
}
struct SymbolicBits
{
  std::uint64_t negative=0,positive=0;
};
SymbolicBits symbolicValue(const HS::SerializedTree& tree,const std::vector<SymbolicBits>& leaves,
                           std::vector<SymbolicBits>& values,std::uint64_t valid)
{
  for(const auto& n : tree.nodes)
  {
    if(n.op == HS::TreeOperator::LEAF) values[n.id] = leaves.at(n.leaf_id);
    else
    {
      SymbolicBits value;
      if(n.op == HS::TreeOperator::MIN) value.positive = valid;
      else value.negative = valid;
      for(const auto id : n.children)
      {
        if(n.op == HS::TreeOperator::MIN)
        {
          value.negative |= values[id].negative;
          value.positive &= values[id].positive;
        }
        else
        {
          value.negative &= values[id].negative;
          value.positive |= values[id].positive;
        }
      }
      values[n.id] = value;
    }
  }
  return values.at(tree.root_id);
}
std::size_t symbolicCost(const HS::SerializedTree& t)
{
  std::size_t n = t.nodes.size(); for(const auto& v : t.nodes) n += v.children.size(); return n;
}
} // namespace

BoundarySignDiagnostics diagnoseBoundaryPlaneSigns(const std::vector<VoxelIndex>& occupied,double voxel_size,
  const BoundaryPlaneRegistry& r,const std::vector<std::pair<std::size_t,int>>& leaf_sources,
  const TreeRepresentation& tree,const BoundaryLocalOptions& limits)
{
  const auto start = Clock::now(); BoundarySignDiagnostics d;
  const auto skip = [&](const char* reason) { d.reason = reason; d.time_ms = ms(start); return d; };
  require(std::isfinite(voxel_size) && voxel_size > 0,"Invalid sign voxel size");
  require(leaf_sources.size() == tree.leaves.size(),"Invalid sign leaf dictionary");
  for(std::size_t i = 0; i < tree.leaves.size(); ++i)
  {
    const auto& source = leaf_sources[i]; const auto& plane = r.supports.at(source.first);
    require(plane.axis < 3 && (source.second == 1 || source.second == -1),"Invalid sign plane reference");
    const double normal[3]{tree.leaves[i].normal.x,tree.leaves[i].normal.y,tree.leaves[i].normal.z};
    for(unsigned a=0;a<3;++a) require(normal[a] == (a == plane.axis ? source.second : 0),"Changed sign leaf normal");
    require(tree.leaves[i].offset == source.second*(static_cast<double>(plane.coordinate)*voxel_size),"Changed sign leaf offset");
  }
  std::array<std::size_t,3> m{}, width{}; std::size_t product = 1;
  for(unsigned a = 0; a < 3; ++a)
  {
    const auto& coordinates=r.coordinates[a];
    require(std::is_sorted(coordinates.begin(),coordinates.end()) &&
            std::adjacent_find(coordinates.begin(),coordinates.end())==coordinates.end(),"Invalid sign coordinates");
    m[a] = coordinates.size(); require(m[a] >= 2,"Empty sign coordinate registry");
    require(m[a] <= (std::numeric_limits<std::size_t>::max()-1)/2,"Sign grid overflow");
    width[a] = 2*m[a]+1;
    if(width[a] > limits.max_strata/product) return skip("sign_strata_limit");
    product *= width[a];
  }
  // Direct deliberately aliases raw/simplified. Identical pointers have exactly
  // the same values on every stratum: avoid evaluating that same DAG twice.
  const bool same_tree=tree.raw==tree.simplified;
  const auto raw = HS::serializeTree(tree.raw);
  const auto simplified = same_tree ? HS::SerializedTree{} : HS::serializeTree(tree.simplified);
  const auto cost = symbolicCost(raw)+(same_tree ? 0 : symbolicCost(simplified))+tree.leaves.size();
  if(cost > limits.max_sign_operations/product) return skip("sign_operation_limit");
  const auto nx = m[0]-1, ny = m[1]-1, nz = m[2]-1;
  std::vector<std::size_t> count(nx*ny*nz,0); // bounded by the larger all-strata product above
  for(const auto& voxel : occupied)
  {
    std::size_t ix[3];
    for(unsigned a = 0; a < 3; ++a)
    {
      const auto& c = r.coordinates[a]; const auto p = std::upper_bound(c.begin(),c.end(),voxel[a]);
      require(p != c.begin() && p != c.end(),"Voxel outside sign grid"); ix[a] = static_cast<std::size_t>(p-c.begin()-1);
    }
    ++count[(ix[0]*ny+ix[1])*nz+ix[2]];
  }
  for(std::size_t x = 0; x < nx; ++x) for(std::size_t y = 0; y < ny; ++y) for(std::size_t z = 0; z < nz; ++z)
  {
    const auto n = count[(x*ny+y)*nz+z]; if(!n) continue;
    const std::size_t ix[3]{x,y,z}; std::size_t capacity = 1;
    for(unsigned a = 0; a < 3; ++a)
    {
      const auto w = static_cast<std::uint64_t>(r.coordinates[a][ix[a]+1]-r.coordinates[a][ix[a]]);
      require(w && w <= n/capacity,"Nonhomogeneous occupied sign cell"); capacity *= static_cast<std::size_t>(w);
    }
    require(capacity == n,"Occupied sign cell volume mismatch");
  }
  const auto filled = [&](std::size_t x,std::size_t y,std::size_t z) {
    if(x == 0 || y == 0 || z == 0 || x == m[0] || y == m[1] || z == m[2]) return false;
    return count[((x-1)*ny+y-1)*nz+z-1] != 0;
  };
  std::vector<std::size_t> leaf_plane; std::vector<unsigned> leaf_axis; std::vector<int> leaf_sign;
  for(const auto& link : leaf_sources)
  {
    const auto& p = r.supports.at(link.first);
    require(link.second == 1 || link.second == -1,"Invalid leaf side");
    leaf_plane.push_back(2*indexOf(r.coordinates[p.axis],p.coordinate)+1); leaf_axis.push_back(p.axis); leaf_sign.push_back(link.second);
  }
  const auto record = [&](const std::array<std::size_t,3>& q,int actual,int after) {
    const auto x=q[0],y=q[1],z=q[2];
    const unsigned dimension = 3-unsigned(x%2)-unsigned(y%2)-unsigned(z%2);
    bool any = false, all = true;
    for(unsigned dx = 0; dx <= x%2; ++dx) for(unsigned dy = 0; dy <= y%2; ++dy) for(unsigned dz = 0; dz <= z%2; ++dz)
    { const bool occupied = filled(x/2+dx,y/2+dy,z/2+dz); any |= occupied; all &= occupied; }
    const int expected = all ? 1 : any ? 0 : -1;
    ++d.strata_count; d.raw_simplified_sign_errors += actual != after;
    d.closure_set_errors += (actual <= 0) != (expected <= 0);
    d.unsafe_free_count += actual < 0 && expected >= 0;
    d.missed_free_count += actual > 0 && expected < 0;
    d.boundary_sign_errors += expected == 0 && actual != 0;
    if(actual == 0 && expected < 0) ++d.free_zeros[dimension];
    if(actual == 0 && expected > 0) ++d.occupied_zeros[dimension];
    if(actual != expected && d.witnesses.size() < 24) d.witnesses.push_back({q,dimension,expected,actual});
  };
  if(!limits.packed_sign_propagation)
  {
    std::vector<int> signs(tree.leaves.size()),raw_values(raw.nodes.size()),simple_values(simplified.nodes.size());
    for(std::size_t x = 0; x < width[0]; ++x) for(std::size_t y = 0; y < width[1]; ++y) for(std::size_t z = 0; z < width[2]; ++z)
    {
      const std::array<std::size_t,3> q{x,y,z};
      for(std::size_t i = 0; i < signs.size(); ++i)
        signs[i] = leaf_sign[i]*(q[leaf_axis[i]] < leaf_plane[i] ? -1 : q[leaf_axis[i]] > leaf_plane[i] ? 1 : 0);
      const int actual = symbolicValue(raw,signs,raw_values);
      record(q,actual,same_tree ? actual : symbolicValue(simplified,signs,simple_values));
    }
  }
  else
  {
    constexpr std::size_t lanes=64;
    std::array<std::array<std::size_t,3>,lanes> queries{};
    std::vector<SymbolicBits> signs(tree.leaves.size()),raw_values(raw.nodes.size()),simple_values(simplified.nodes.size());
    for(std::size_t begin=0;begin<product;begin+=lanes)
    {
      const auto lane_count=std::min(lanes,product-begin);
      const std::uint64_t valid=lane_count==lanes ? ~std::uint64_t{0} : (std::uint64_t{1}<<lane_count)-1;
      for(std::size_t lane=0;lane<lane_count;++lane)
      {
        auto index=begin+lane;
        queries[lane][2]=index%width[2]; index/=width[2];
        queries[lane][1]=index%width[1]; queries[lane][0]=index/width[1];
      }
      std::fill(signs.begin(),signs.end(),SymbolicBits{});
      for(std::size_t i=0;i<signs.size();++i)
      {
        for(std::size_t lane=0;lane<lane_count;++lane)
        {
          const auto coordinate=queries[lane][leaf_axis[i]];
          if(coordinate==leaf_plane[i]) continue;
          const bool positive=leaf_sign[i]*(coordinate<leaf_plane[i] ? -1 : 1)>0;
          (positive ? signs[i].positive : signs[i].negative)|=std::uint64_t{1}<<lane;
        }
      }
      const auto actual=symbolicValue(raw,signs,raw_values,valid);
      const auto after=same_tree ? actual : symbolicValue(simplified,signs,simple_values,valid);
      require(!(actual.negative&actual.positive) && !(after.negative&after.positive),"Invalid packed symbolic sign");
      for(std::size_t lane=0;lane<lane_count;++lane)
      {
        const auto bit=std::uint64_t{1}<<lane;
        const int a=actual.positive&bit ? 1 : actual.negative&bit ? -1 : 0;
        const int b=after.positive&bit ? 1 : after.negative&bit ? -1 : 0;
        record(queries[lane],a,b);
      }
    }
  }
  d.applied = true; d.closure_set_certified = d.closure_set_errors == 0 && d.raw_simplified_sign_errors == 0;
  d.strict_sign_certified = d.closure_set_certified && d.unsafe_free_count == 0 && d.missed_free_count == 0 &&
    d.boundary_sign_errors == 0 && std::accumulate(d.free_zeros.begin(),d.free_zeros.end(),std::size_t(0)) == 0 &&
    std::accumulate(d.occupied_zeros.begin(),d.occupied_zeros.end(),std::size_t(0)) == 0;
  d.reason = d.strict_sign_certified ? "strict_sign_certified" : "remaining_sign_errors"; d.time_ms = ms(start); return d;
}
static BoundarySignDiagnostics diagnoseSigns(const std::vector<VoxelIndex>& occupied,double h,const BoundaryTreeResult& r,
                                              const TreeRepresentation& tree,const BoundaryLocalOptions& limits)
{
  if(r.fallback) { BoundarySignDiagnostics d; d.reason="obstacle_fallback_not_closure_tree"; return d; }
  require(tree.leaves.size()==r.tree.leaves.size() && r.leaf_sources.size()==tree.leaves.size(),"Changed local leaf dictionary");
  for(std::size_t i=0;i<tree.leaves.size();++i)
  {
    const auto& a=tree.leaves[i]; const auto& b=r.tree.leaves[i];
    require(a.normal.x==b.normal.x && a.normal.y==b.normal.y && a.normal.z==b.normal.z && a.offset==b.offset,
            "Local construction changed plane coefficients");
  }
  return diagnoseBoundaryPlaneSigns(occupied,h,{r.coordinates,r.supports},r.leaf_sources,tree,limits);
}
BoundarySignDiagnostics diagnoseBoundarySigns(const OctreeOccupancy& o,const BoundaryTreeResult& r,
                                              const TreeRepresentation& t,const BoundaryLocalOptions& limits)
{ return diagnoseSigns(o.occupied_voxels,o.voxel_size,r,t,limits); }
BoundarySignDiagnostics diagnoseBoundarySigns(const SparseVoxelOccupancy& o,const BoundaryTreeResult& r,
                                              const TreeRepresentation& t,const BoundaryLocalOptions& limits)
{ return diagnoseSigns(o.occupied_voxels,o.voxel_size,r,t,limits); }
namespace
{
using Expr = HS::ExpressionPtr;
struct Bounds { std::array<std::optional<std::int64_t>,3> lo{},hi{}; };
struct Region { Expr expression; Bounds bounds; std::vector<std::size_t> sources; bool active = true; };
struct Face { std::size_t region; int sign; std::int64_t u0,u1,v0,v1; bool exterior; };
// nullptr is a transient TRUE (-infinity) operand, never exported as a node or plane.
Expr combine(HS::TreeOperator op,const std::vector<Expr>& input)
{
  std::vector<Expr> out; std::set<int> leaf_ids; std::set<const HS::Expression*> other_ids;
  const auto insert = [&](const Expr& e) {
    if(e->op == HS::TreeOperator::LEAF ? leaf_ids.insert(e->leaf_id).second : other_ids.insert(e.get()).second) out.push_back(e);
  };
  for(const auto& e : input)
  {
    if(!e) { if(op == HS::TreeOperator::MIN) return {}; else continue; }
    if(e->op == op) for(const auto& child : e->children) insert(child);
    else insert(e);
  }
  if(out.empty()) return {};
  return out.size() == 1 ? out.front() : HS::makeNode(op,std::move(out));
}
std::vector<Expr> factors(const Expr& e)
{ return e->op == HS::TreeOperator::MAX ? e->children : std::vector<Expr>{e}; }
Bounds envelope(const Bounds& a,const Bounds& b)
{
  Bounds out;
  for(unsigned k = 0; k < 3; ++k)
  {
    if(a.lo[k] && b.lo[k]) out.lo[k] = std::min(*a.lo[k],*b.lo[k]);
    if(a.hi[k] && b.hi[k]) out.hi[k] = std::max(*a.hi[k],*b.hi[k]);
  }
  return out;
}
struct Counts { std::size_t nodes = 0, leaves = 0, depth = 0; };
Counts countReferences(const Expr& e,std::size_t limit)
{
  std::map<const HS::Expression*,Counts> cache;
  const auto visit = [&](const auto& self,const Expr& p)->Counts {
    if(!p) return {};
    const auto found = cache.find(p.get()); if(found != cache.end()) return found->second;
    Counts c{1,p->op == HS::TreeOperator::LEAF ? 1u : 0u,1};
    for(const auto& child : p->children)
    {
      const auto v = self(self,child);
      c.nodes += std::min(v.nodes,limit+1-std::min(c.nodes,limit+1));
      c.leaves += std::min(v.leaves,limit+1-std::min(c.leaves,limit+1)); c.depth = std::max(c.depth,v.depth+1);
    }
    cache.emplace(p.get(),c); return c;
  };
  return visit(visit,e);
}
} // namespace

BoundaryLocalResult buildBoundaryLocalTree(const BoundaryTreeResult& reference,const BoundaryLocalOptions& options)
{
  const auto core = Clock::now(); BoundaryLocalResult r; r.tree = reference.tree; r.obstacle_fallback = reference.fallback;
  require(options.max_attempts && options.max_adjacency_checks && options.max_expanded_nodes &&
          options.max_expanded_nodes < std::numeric_limits<std::size_t>::max(),"Invalid local construction limits");
  if(reference.fallback) { r.stop_reason = "flat_work_cap_fallback"; r.core_ms = ms(core); return r; }
  const auto graph_start = Clock::now();
  std::map<std::pair<unsigned,std::int64_t>,std::size_t> support_ids;
  for(std::size_t i = 0; i < reference.supports.size(); ++i)
    support_ids.emplace(std::make_pair(reference.supports[i].axis,reference.supports[i].coordinate),i);
  std::map<std::pair<std::size_t,int>,int> leaf_ids;
  for(std::size_t i = 0; i < reference.leaf_sources.size(); ++i) leaf_ids.emplace(reference.leaf_sources[i],static_cast<int>(i));
  const auto leaf = [&](unsigned a,std::int64_t c,int sign) {
    return HS::makeLeaf(leaf_ids.at({support_ids.at({a,c}),sign}));
  };
  const auto boundsTree = [&](const Bounds& b) {
    std::vector<Expr> f;
    for(unsigned a = 0; a < 3; ++a)
    { if(b.lo[a]) f.push_back(leaf(a,*b.lo[a],-1)); if(b.hi[a]) f.push_back(leaf(a,*b.hi[a],1)); }
    return combine(HS::TreeOperator::MAX,f);
  };
  std::vector<Region> regions;
  std::map<std::size_t,std::array<std::vector<Face>,2>> buckets;
  for(const auto& p : reference.prisms)
  {
    Bounds b; for(unsigned a = 0; a < 3; ++a) { b.lo[a] = p.lower[a]; b.hi[a] = p.upper[a]; }
    const auto id = regions.size(); regions.push_back({boundsTree(b),b,{id},true});
    for(unsigned a = 0; a < 3; ++a) for(unsigned high = 0; high < 2; ++high)
    {
      const unsigned u = (a+1)%3,v = (a+2)%3; const auto c = high ? p.upper[a] : p.lower[a];
      buckets[support_ids.at({a,c})][high].push_back({id,high ? 1 : -1,p.lower[u],p.upper[u],p.lower[v],p.upper[v],false});
    }
  }
  for(unsigned a = 0; a < 3; ++a) for(unsigned high = 0; high < 2; ++high)
  {
    const auto c = high ? reference.coordinates[a].back() : reference.coordinates[a].front();
    Bounds b; if(high) b.lo[a] = c; else b.hi[a] = c;
    const auto id = regions.size(); regions.push_back({boundsTree(b),b,{id},true});
    const unsigned u = (a+1)%3,v = (a+2)%3;
    // The exterior is unbounded. This rectangle clips its adjacency search ONLY,
    // not its geometry; all finite prisms lie within these transverse extents.
    buckets[support_ids.at({a,c})][high ? 0 : 1].push_back({id,high ? -1 : 1,
      reference.coordinates[u].front(),reference.coordinates[u].back(),reference.coordinates[v].front(),reference.coordinates[v].back(),true});
  }
  bool graph_complete = true;
  for(const auto& bucket : buckets)
  {
    for(const auto& a : bucket.second[0])
    {
      for(const auto& b : bucket.second[1])
      {
        if(r.adjacency_checks == options.max_adjacency_checks) { graph_complete = false; break; }
        ++r.adjacency_checks;
        const auto u0 = std::max(a.u0,b.u0),u1 = std::min(a.u1,b.u1),v0 = std::max(a.v0,b.v0),v1 = std::min(a.v1,b.v1);
        if(u0 < u1 && v0 < v1) r.adjacency.push_back({a.region,b.region,bucket.first,a.sign,u0,u1,v0,v1,a.exterior || b.exterior});
      }
      if(!graph_complete) break;
    }
    if(!graph_complete) break;
  }
  std::sort(r.adjacency.begin(),r.adjacency.end(),[](const auto& a,const auto& b) {
    const auto area = [](const auto& f) { return static_cast<long double>(f.u1-f.u0)*static_cast<long double>(f.v1-f.v0); };
    if(area(a) != area(b)) return area(a) > area(b);
    return std::tie(a.support,a.a,a.b,a.u0,a.v0) < std::tie(b.support,b.a,b.b,b.u0,b.v0);
  });
  r.adjacency_ms = ms(graph_start); const auto construction_start = Clock::now();
  std::vector<std::size_t> owner(regions.size()); std::iota(owner.begin(),owner.end(),0);
  std::set<std::tuple<std::size_t,std::size_t,std::size_t>> attempted;
  std::size_t total_refs = 1;
  for(const auto& region : regions) total_refs += countReferences(region.expression,options.max_expanded_nodes).nodes;
  if(!graph_complete) r.stop_reason = "adjacency_check_limit";
  else if(total_refs > options.max_expanded_nodes) r.stop_reason = "initial_expanded_node_limit";
  bool progress = r.stop_reason.empty();
  while(progress)
  {
    progress = false;
    for(const auto& face : r.adjacency)
    {
      const auto a_id = owner[face.a], b_id = owner[face.b]; if(a_id == b_id) continue;
      if(!attempted.emplace(a_id,b_id,face.support).second) continue;
      if(r.attempts == options.max_attempts) { r.stop_reason = "local_attempt_limit"; break; }
      ++r.attempts;
      const auto a = factors(regions[a_id].expression),b = factors(regions[b_id].expression);
      const int a_gate = leaf_ids.at({face.support,face.a_sign}),b_gate = leaf_ids.at({face.support,-face.a_sign});
      const auto has = [](const auto& f,int id) { return std::any_of(f.begin(),f.end(),[&](const auto& p) { return p->op == HS::TreeOperator::LEAF && p->leaf_id == id; }); };
      if(!has(a,a_gate) || !has(b,b_gate)) { ++r.factoring_failed; continue; }
      const auto common = envelope(regions[a_id].bounds,regions[b_id].bounds);
      const auto residual = [&](const std::vector<Expr>& f,int gate) {
        std::vector<Expr> out;
        for(const auto& p : f)
        {
          if(p->op == HS::TreeOperator::LEAF)
          {
            if(p->leaf_id == gate) continue;
            const auto& link = reference.leaf_sources.at(p->leaf_id); const auto& s = reference.supports.at(link.first);
            if(link.second < 0 && common.lo[s.axis] && *common.lo[s.axis] >= s.coordinate) continue;
            if(link.second > 0 && common.hi[s.axis] && *common.hi[s.axis] <= s.coordinate) continue;
          }
          out.push_back(p);
        }
        return combine(HS::TreeOperator::MAX,out);
      };
      const auto ra = residual(a,a_gate),rb = residual(b,b_gate),s = HS::makeLeaf(a_gate),neg_s = HS::makeLeaf(b_gate);
      // min(max(ra,s),max(rb,-s)) = max(min(ra,rb),min(ra,-s),min(s,rb),min(s,-s)).
      // Dropping the last MAX child preserves <=0 since min(s,-s)<=0 everywhere.
      const auto expression = combine(HS::TreeOperator::MAX,{boundsTree(common),
        combine(HS::TreeOperator::MIN,{ra,rb}),combine(HS::TreeOperator::MIN,{ra,neg_s}),combine(HS::TreeOperator::MIN,{s,rb})});
      if(!expression) { ++r.factoring_failed; continue; } // no constant node or fabricated plane
      const auto next = countReferences(expression,options.max_expanded_nodes);
      const auto old_a = countReferences(regions[a_id].expression,options.max_expanded_nodes);
      const auto old_b = countReferences(regions[b_id].expression,options.max_expanded_nodes);
      if(next.depth > 128 || next.nodes > options.max_expanded_nodes ||
         total_refs-old_a.nodes-old_b.nodes > options.max_expanded_nodes-next.nodes)
      { ++r.factoring_failed; continue; }
      total_refs = total_refs-old_a.nodes-old_b.nodes+next.nodes;
      auto sources = regions[a_id].sources; sources.insert(sources.end(),regions[b_id].sources.begin(),regions[b_id].sources.end());
      regions[a_id].active = regions[b_id].active = false; const auto id = regions.size();
      for(const auto source : sources) owner[source] = id;
      regions.push_back({expression,common,std::move(sources),true}); ++r.logical_merges; ++r.cancelled_interface_pairs; progress = true;
    }
    if(!r.stop_reason.empty()) break;
  }
  if(r.stop_reason.empty()) r.stop_reason = "no_more_matching_root_gate_templates";
  // Keep the original expression byte-compatible if no rewrite was accepted.
  if(r.logical_merges)
  {
    std::vector<Expr> roots; for(const auto& region : regions) if(region.active) roots.push_back(region.expression);
    r.tree.raw = combine(HS::TreeOperator::MIN,roots);
  }
  r.constructed_nodes = HS::serializeTree(r.tree.raw).nodes.size();
  const auto counts = countReferences(r.tree.raw,std::numeric_limits<std::size_t>::max()-1);
  r.expanded_nodes = counts.nodes; r.tree.leaf_references = counts.leaves;
  r.construction_ms = ms(construction_start); const auto simplify_start = Clock::now();
  if(r.logical_merges) r.tree.simplified = HS::simplifyToFixedPoint(r.tree.raw);
  r.simplify_ms = ms(simplify_start); r.core_ms = ms(core); return r;
}

bool boundaryLocalStrictFree(const BoundaryLocalResult& r,const OctreeOccupancy& occupancy,const HS::Vec3& p,bool* zero)
{
  require(std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z),"Nonfinite local query");
  if(zero) *zero = false;
  const double q = HS::evaluate(r.tree.simplified,r.tree.leaves,p);
  if(r.obstacle_fallback) return q > 0;
  if(r.diagnostics.strict_sign_certified) return q < 0;
  if(q != 0) return q < 0;
  if(zero) *zero = true;
  return !occupancy.contains(p);
}

void writeBoundaryLocalReport(const std::filesystem::path& path,const BoundaryTreeResult& reference,
                             const BoundarySignDiagnostics& flat,const BoundaryLocalResult* local)
{
  std::ofstream out(path); require(bool(out),"Cannot write boundary sign diagnostics");
  out.exceptions(std::ios::failbit | std::ios::badbit); out << std::setprecision(17) << std::boolalpha;
  const auto array = [&](const auto& values) { out << '['; bool first = true; for(const auto& v : values) { if(!first) out << ','; first = false; out << v; } out << ']'; };
  const auto report = [&](const BoundarySignDiagnostics& d) {
    out << "{\"applied\":" << d.applied << ",\"strict_sign_certified\":" << d.strict_sign_certified
      << ",\"closure_set_certified\":" << d.closure_set_certified << ",\"reason\":\"" << d.reason << "\",\"strata_count\":" << d.strata_count
      << ",\"unsafe_free_count\":" << d.unsafe_free_count << ",\"missed_free_count\":" << d.missed_free_count
      << ",\"boundary_sign_errors\":" << d.boundary_sign_errors << ",\"closure_set_errors\":" << d.closure_set_errors
      << ",\"raw_simplified_sign_errors\":" << d.raw_simplified_sign_errors << ",\"free_zeros_by_dimension\":";
    array(d.free_zeros); out << ",\"occupied_zeros_by_dimension\":"; array(d.occupied_zeros);
    out << ",\"time_ms\":" << d.time_ms << ",\"witnesses\":[";
    for(std::size_t i = 0; i < d.witnesses.size(); ++i)
    {
      const auto& w = d.witnesses[i]; out << (i ? "," : "") << "{\"strata\":"; array(w.strata);
      out << ",\"dimension\":" << w.dimension << ",\"expected_sign\":" << w.expected << ",\"actual_sign\":" << w.actual << '}';
    }
    out << "]}";
  };
  out << "{\"schema\":1,\"method\":\"integer_plane_arrangement_three_sign_propagation\",\"dimension_order\":[\"vertex\",\"edge\",\"face\",\"volume\"],"
      << "\"stratum_index_rule\":\"2*k+1 is plane k; even indices are open intervals, including unbounded intervals\",\"coordinates_grid\":[";
  for(unsigned a = 0; a < 3; ++a) { if(a) out << ','; array(reference.coordinates[a]); }
  out << "],\"flat\":"; report(flat); out << ",\"local\":";
  if(local) report(local->diagnostics); else out << "null";
  if(local)
  {
    out << ",\"stop_reason\":\"" << local->stop_reason << "\",\"adjacency_faces\":[";
    for(std::size_t i = 0; i < local->adjacency.size(); ++i)
    {
      const auto& f = local->adjacency[i];
      out << (i ? "," : "") << "{\"a\":" << f.a << ",\"b\":" << f.b << ",\"support_id\":" << f.support
        << ",\"a_sign\":" << f.a_sign << ",\"uv_bounds\":[" << f.u0 << ',' << f.u1 << ',' << f.v0 << ',' << f.v1
        << "],\"exterior\":" << f.exterior << '}';
    }
    out << ']';
  }
  out << "}\n"; out.close();
}
} // namespace rokae_demo
