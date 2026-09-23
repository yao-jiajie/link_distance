#include <rokae_demo/halfspace_compiled_tree.hpp>
#include <rokae_demo/halfspace_validation_cache.hpp>
#include <atomic>
#include <cstring>
#include <iostream>
#include <random>
#include <thread>

namespace HS=rokae_demo::halfspace;
void check(bool v,const char* why) { if(!v) throw std::runtime_error(why); }
template<class F> void rejects(F f)
{ bool rejected=false; try { f(); } catch(const std::exception&) { rejected=true; } check(rejected,"Invalid compiled input accepted"); }
bool same(double a,double b) { return std::memcmp(&a,&b,sizeof(double))==0 || (std::isnan(a) && std::isnan(b)); }
int main()
{
  try
  {
    std::vector<HS::HalfspaceLeaf> leaves(8);
    for(int i=0;i<8;++i) { leaves[i].normal={.25*(i-3),.125*(i%3-1),i%2 ? 1. : -1.}; leaves[i].offset=.04*(i-2); }
    leaves[0].normal={1,0,0}; leaves[0].offset=0;
    std::vector<HS::ExpressionPtr> nodes; for(int i=0;i<8;++i) nodes.push_back(HS::makeLeaf(i));
    std::mt19937_64 rng(19821);
    for(unsigned i=0;i<22;++i)
      nodes.push_back(HS::makeNode(i%2 ? HS::TreeOperator::MIN : HS::TreeOperator::MAX,
        {nodes.back(),nodes[rng()%nodes.size()],nodes[rng()%nodes.size()]}));
    const auto root=nodes.back(); const auto serialized=HS::serializeTree(root);
    const HS::CompiledTree program(root,leaves), imported(serialized,leaves),interned(serialized,leaves,true);
    auto scratch=program.makeWorkspace(),other=imported.makeWorkspace(),interned_ws=interned.makeWorkspace();
    const auto* storage=scratch.data();
    const auto verify=[&](const HS::Vec3& p) {
      const auto expected=HS::evaluate(root,leaves,p);
      check(same(expected,program.evaluate(p,scratch)),"DAG/recursive scalar bits differ");
      check(same(expected,imported.evaluate(p,other)),"Imported DAG differs");
      check(same(expected,interned.evaluate(p,interned_ws)),"Interned DAG scalar bits differ");
      check(scratch.data()==storage,"Workspace reallocated");
    };
    for(int i=0;i<400;++i) verify({.01*static_cast<int>(rng()%501)-2.5,.03*static_cast<int>(rng()%121)-1.8,.04*static_cast<int>(rng()%91)-1.8});
    for(double v:{0.,-0.,std::numeric_limits<double>::denorm_min(),-std::numeric_limits<double>::denorm_min(),
                  std::numeric_limits<double>::max(),-std::numeric_limits<double>::max()}) verify({v,v,v});
    std::vector<HS::Vec3> batch_points;
    for(unsigned i=0;i<257;++i)
      batch_points.push_back({.01*static_cast<int>(i%31)-.15,.02*static_cast<int>(i%17)-.16,.03*static_cast<int>(i%13)-.18});
    batch_points.push_back({-0.,0.,std::numeric_limits<double>::denorm_min()});
    std::vector<double> batch_expected(batch_points.size()),batch_values(batch_points.size());
    for(std::size_t i=0;i<batch_points.size();++i) batch_expected[i]=program.evaluate(batch_points[i],scratch);
    for(std::size_t workers_count:{1u,4u})
    {
      auto batch=program.makeBatchWorkspace(workers_count); const auto bytes=batch.storageBytes();
      check(batch.workerCount()==workers_count && batch.valuesPerWorker()==program.nodeCount(),"Compiled batch workspace shape mismatch");
      program.evaluateBatch(batch_points,batch_values,batch);
      for(std::size_t i=0;i<batch_values.size();++i)
        check(same(batch_values[i],batch_expected[i]),"Compiled batch changed scalar bits/order");
      std::fill(batch_values.begin(),batch_values.end(),NAN); program.evaluateBatch(batch_points,batch_values,batch);
      check(batch.storageBytes()==bytes,"Compiled batch workspace storage changed during reuse");
    }
    auto empty_batch=program.makeBatchWorkspace(2); std::vector<HS::Vec3> no_points; std::vector<double> no_values;
    program.evaluateBatch(no_points,no_values,empty_batch);
    rejects([&] { program.makeBatchWorkspace(0); });
    auto moved_batch=std::move(empty_batch); rejects([&] { program.evaluateBatch(batch_points,batch_values,empty_batch); });
    program.evaluateBatch(batch_points,batch_values,moved_batch);
    // Repeated/alternating queries cannot reuse stale values. Same immutable
    // program is safe across threads with separate caller-owned workspaces.
    std::atomic<bool> failed{false}; std::vector<std::thread> workers;
    for(unsigned worker=0;worker<4;++worker) workers.emplace_back([&,worker] {
      auto ws=interned.makeWorkspace();
      for(unsigned i=0;i<100;++i)
      {
        const HS::Vec3 p{.04*i,.02*worker,-.01*i};
        if(!same(interned.evaluate(p,ws),HS::evaluate(root,leaves,p))) failed=true;
      }
    });
    for(auto& worker:workers) worker.join();
    check(!failed,"Concurrent compiled scalar mismatch");
    auto wrong=scratch; wrong.pop_back(); rejects([&] { program.evaluate({},wrong); });
    rejects([&] { program.evaluate({INFINITY,0,0},scratch); });
    rejects([&] { program.evaluate({NAN,0,0},scratch); });
    // Captured coefficients are immutable snapshots, not dangling references.
    const HS::CompiledTree single(HS::makeLeaf(0),leaves); auto single_ws=single.makeWorkspace();
    leaves[0].offset=123; check(single.evaluate({1,0,0},single_ws)==1,"Coefficient snapshot changed");
    leaves[0].offset=0;
    for(unsigned mode=0;mode<9;++mode)
    {
      auto bad=serialized;
      if(mode==0) bad.root_id=-1;
      if(mode==1) bad.nodes[0].id=8;
      if(mode==2) bad.nodes[0].leaf_id=-1;
      if(mode==3) bad.nodes[0].children={0};
      if(mode==4) bad.nodes.back().children.clear();
      if(mode==5) bad.nodes.back().children={bad.root_id};
      if(mode==6) bad.nodes.back().children={-1};
      if(mode==7) bad.nodes.back().op=static_cast<HS::TreeOperator>(90);
      if(mode==8) bad.nodes.back().children={0}; // now other nodes are unreachable
      rejects([&] { HS::CompiledTree rejected(bad,leaves); });
      rejects([&] { HS::CompiledTree rejected(bad,leaves,true); });
    }
    auto bad_leaves=leaves; bad_leaves[0].offset=INFINITY;
    rejects([&] { HS::CompiledTree rejected(serialized,bad_leaves); });
    rejects([&] { HS::CompiledTree rejected(HS::ExpressionPtr{},leaves); });
    // A deeply SHARED expression would require exponential recursive work.
    // Compilation and evaluation visit its linear-size DAG without expansion.
    auto shared=HS::makeLeaf(0);
    for(unsigned i=0;i<48;++i) shared=HS::makeNode(HS::TreeOperator::MIN,{shared,shared});
    const HS::CompiledTree compact(shared,leaves); auto ws=compact.makeWorkspace();
    check(compact.nodeCount()==49 && compact.edgeCount()==96 && compact.evaluate({.5,0,0},ws)==.5,"Shared DAG expanded");
    std::vector<double> wrong_output(batch_points.size()-1);
    rejects([&] { program.evaluateBatch(batch_points,wrong_output,empty_batch); });
    auto compact_batch=compact.makeBatchWorkspace(2);
    rejects([&] { program.evaluateBatch(batch_points,batch_values,compact_batch); });
    auto invalid_points=batch_points; invalid_points.back().x=INFINITY;
    rejects([&] { program.evaluateBatch(invalid_points,batch_values,moved_batch); });
    program.evaluateBatch(batch_points,batch_values,moved_batch);
    for(std::size_t i=0;i<batch_points.size();++i)
      check(same(batch_values[i],batch_expected[i]),"Compiled batch failed after worker exception");
    auto excess_workers=program.makeBatchWorkspace(8);
    for(std::size_t count:{1u,3u})
    {
      std::vector<HS::Vec3> small(batch_points.begin(),batch_points.begin()+count); std::vector<double> values(count);
      program.evaluateBatch(small,values,excess_workers);
      for(std::size_t i=0;i<count;++i) check(same(values[i],batch_expected[i]),"Compiled small batch lost points");
    }
    // Different pointers with identical ordered expressions can share execution;
    // permutations, arities and operators must NOT be treated as identical.
    const auto pair=[&](HS::TreeOperator op,bool reverse=false,bool repeat=false) {
      std::vector<HS::ExpressionPtr> c{HS::makeLeaf(reverse ? 1 : 0),HS::makeLeaf(reverse ? 0 : 1)};
      if(repeat) c.push_back(HS::makeLeaf(1));
      return HS::makeNode(op,c);
    };
    const auto clones=HS::makeNode(HS::TreeOperator::MAX,{pair(HS::TreeOperator::MIN),pair(HS::TreeOperator::MIN),
      pair(HS::TreeOperator::MIN,true),pair(HS::TreeOperator::MIN,false,true),pair(HS::TreeOperator::MAX)});
    const auto original=HS::serializeTree(clones);
    const HS::CompiledTree plain(original,leaves),cse(original,leaves,true);
    check(!plain.structuralSharingEnabled() && cse.structuralSharingEnabled(),"Compiled sharing mode was not retained");
    check(cse.nodeCount()==7 && cse.edgeCount()==14 && cse.sourceNodeCount()==plain.nodeCount(),"Incorrect ordered CSE grouping");
    check(cse.sharedNodeCount()==plain.nodeCount()-7,"Incorrect shared node count");
    auto pw=plain.makeWorkspace(),cw=cse.makeWorkspace(); const auto storage_cse=cw.data();
    for(double x:{-0.,0.,-.01,.02,1e300}) for(double y:{-0.,0.,.03,-.04,1e300})
      check(same(plain.evaluate({x,y,0},pw),cse.evaluate({x,y,0},cw)),"CSE changed scalar/signed zero");
    check(cw.data()==storage_cse && original.nodes.size()==HS::serializeTree(clones).nodes.size(),"CSE mutated source or reallocated query workspace");
    auto broken=original; broken.nodes.back().children={0};
    rejects([&] { HS::CompiledTree rejected(broken,leaves,true); }); // validate BEFORE CSE
    auto cache_program=std::make_shared<const HS::CompiledTree>(original,leaves,true);
    std::weak_ptr<const HS::CompiledTree> retained=cache_program;
    HS::ValidationValueCache cache(cache_program,2); auto cache_work=cache_program->makeWorkspace();
    const HS::Vec3 positive_zero{0.,0.,0},negative_zero{-0.,0.,0};
    check(same(cache.evaluate(positive_zero,cache_work),cache.evaluate(positive_zero,cache_work)) &&
          cache.hits()==1 && cache.evaluations()==1 && cache.size()==1,"Validation cache missed an exact coordinate");
    cache.evaluate(negative_zero,cache_work);
    check(cache.evaluations()==2 && cache.size()==2,"Validation cache merged signed-zero coordinates");
    const HS::Vec3 adjacent{std::nextafter(0.,1.),0.,0.};
    cache.evaluate(adjacent,cache_work); cache.evaluate(adjacent,cache_work);
    check(cache.evaluations()==4 && cache.hits()==1 && cache.size()==2,"Validation cache limit was not enforced");
    auto short_work=cache_work; short_work.pop_back();
    rejects([&] { cache.evaluate({},short_work); });
    rejects([&] { cache.evaluate({INFINITY,0,0},cache_work); });
    cache_program.reset(); check(!retained.expired() && &cache.program()==retained.lock().get(),"Validation cache lost its immutable snapshot");
    std::cout<<"Compiled tree: scalar/batch bits, validation cache, malformed DAGs, workspaces, threads and sharing passed\n";
  }
  catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
