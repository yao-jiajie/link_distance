#include <rokae_demo/halfspace_smooth_tree.hpp>
#include <iostream>
#include <cstring>
#include <random>

namespace HS=rokae_demo::halfspace;
void check(bool b,const char* why) { if(!b) throw std::runtime_error(why); }
bool same(double a,double b) { return std::memcmp(&a,&b,sizeof(double))==0; }
template<class F> void rejects(F f)
{ bool caught=false; try { f(); } catch(const std::exception&) { caught=true; } check(caught,"Invalid smooth input accepted"); }
int main()
{
  try
  {
    std::vector<HS::HalfspaceLeaf> leaves(3);
    leaves[0].normal={1,0,0}; leaves[1].normal={0,1,0}; leaves[2].normal={0,0,1};
    const auto a=HS::makeLeaf(0),b=HS::makeLeaf(1),c=HS::makeLeaf(2);
    for(auto op:{HS::TreeOperator::MIN,HS::TreeOperator::MAX})
      for(double beta:{.001,10.,1000.,1000000.})
      {
        const auto tree=HS::makeNode(op,{a,b,c}); HS::SmoothTree smooth(tree,leaves,{beta,0});
        auto sw=smooth.makeValueWorkspace(); auto gw=smooth.makeGradientWorkspace();
        for(const auto p:{HS::Vec3{0,0,0},HS::Vec3{.01,.01,.01},HS::Vec3{-.04,.03,.01},HS::Vec3{1e4,-1e4,0}})
        {
          const double d=smooth.evaluate(p,sw),q=HS::evaluate(tree,leaves,p); const auto s=smooth.evaluateWithGradient(p,gw);
          check(d==s.value && d<=-q,"Smooth scalar/value-gradient or upper ordering mismatch");
          check(-q-d<=smooth.errorBound()+1e-9,"Smooth node error bound failed");
          check(std::abs(s.gradient.x+s.gradient.y+s.gradient.z+1)<1e-14,"Weights are not a convex combination");
        }
        const auto at_tie=smooth.evaluateWithGradient({0,0,0},gw);
        check(std::abs(at_tie.gradient.x+1./3)<1e-14,"Tie gradient incorrect");
        check(std::abs(at_tie.value-(op==HS::TreeOperator::MIN ? 0 : -std::log(3.)/beta))<1e-12,"Tie normalization incorrect");
      }
    // Preserve duplicate occurrences: MAX(a,a) has log(2)/beta, MIN(a,a) does not.
    auto shared=HS::makeNode(HS::TreeOperator::MAX,{a,a}); HS::SmoothTree duplicate(shared,leaves,{100,0});
    auto dw=duplicate.makeValueWorkspace(); check(duplicate.nodeCount()==2,"Shared leaf was expanded");
    check(std::abs(duplicate.evaluate({},dw)+std::log(2.)/100)<1e-16,"Duplicate MAX term discarded");
    auto nested=HS::makeNode(HS::TreeOperator::MIN,{shared,HS::makeNode(HS::TreeOperator::MAX,{b,c})});
    HS::SmoothTree smooth(nested,leaves,{0,.001}); auto sw=smooth.makeValueWorkspace(); auto gw=smooth.makeGradientWorkspace();
    check(smooth.errorBound()<=.00100000000001,"Automatic error budget violated");
    for(int i=-4;i<=4;++i) for(int j=-4;j<=4;++j)
    {
      const HS::Vec3 p{.01*i,.01*j,.003}; const auto s=smooth.evaluateWithGradient(p,gw);
      const double hard=-HS::evaluate(nested,leaves,p);
      check(s.value<=hard && hard-s.value<=smooth.errorBound()+1e-14,"Nested upper/error bound failed");
      const double coords[3]{p.x,p.y,p.z},g[3]{s.gradient.x,s.gradient.y,s.gradient.z},h=1e-7;
      for(unsigned axis=0;axis<3;++axis)
      {
        double lo[3]{coords[0],coords[1],coords[2]},hi[3]{coords[0],coords[1],coords[2]}; lo[axis]-=h; hi[axis]+=h;
        const double fd=(smooth.evaluate({hi[0],hi[1],hi[2]},sw)-smooth.evaluate({lo[0],lo[1],lo[2]},sw))/(2*h);
        check(std::abs(fd-g[axis])<2e-6,"Analytic nested gradient mismatch");
      }
    }
    auto bad=sw; bad.pop_back(); rejects([&] { smooth.evaluate({},bad); });
    rejects([&] { smooth.evaluate({NAN,0,0},sw); });
    for(const auto options:{HS::SmoothTreeOptions{},HS::SmoothTreeOptions{1,.1},HS::SmoothTreeOptions{-1,0},HS::SmoothTreeOptions{INFINITY,0}})
      rejects([&] { HS::SmoothTree s(nested,leaves,options); });
    HS::SmoothTree extreme(nested,leaves,{1e300,0}); auto ew=extreme.makeValueWorkspace();
    check(std::isfinite(extreme.evaluate({1e150,-1e150,0},ew)),"Stable shift overflowed");
    const auto cloned_pair=[&] { return HS::makeNode(HS::TreeOperator::MIN,{HS::makeLeaf(0),HS::makeLeaf(1)}); };
    const auto clones=HS::makeNode(HS::TreeOperator::MAX,{cloned_pair(),cloned_pair(),
      HS::makeNode(HS::TreeOperator::MIN,{HS::makeLeaf(1),HS::makeLeaf(0)})});
    for(const auto parameters:{HS::SmoothTreeOptions{0,.001},HS::SmoothTreeOptions{10,0},HS::SmoothTreeOptions{1e6,0}})
    {
      auto sharing=parameters; sharing.share_subexpressions=true;
      HS::SmoothTree before(clones,leaves,parameters),after(clones,leaves,sharing);
      check(after.nodeCount()==5 && after.sourceNodeCount()==10 && after.sharedNodeCount()==5,"Smooth CSE counts incorrect");
      check(after.expCallsPerQuery()==4 && after.logCallsPerQuery()==3,"Smooth transcendental counts incorrect");
      check(same(before.beta(),after.beta()) && same(before.errorBound(),after.errorBound()) &&
        same(before.errorWeight(),after.errorWeight()) && same(before.lipschitzBound(),after.lipschitzBound()),"CSE changed smoothing parameters/bounds");
      auto bw=before.makeValueWorkspace(),aw=after.makeValueWorkspace(); auto bg=before.makeGradientWorkspace(),ag=after.makeGradientWorkspace();
      const auto storage=ag.data();
      for(double x:{-0.,0.,-.01,.03,1e150}) for(double y:{-0.,0.,.02,-.04,-1e150})
      {
        const HS::Vec3 p{x,y,0}; const auto b=before.evaluateWithGradient(p,bg),a=after.evaluateWithGradient(p,ag);
        check(same(before.evaluate(p,bw),after.evaluate(p,aw)) && same(b.value,a.value) &&
          same(b.gradient.x,a.gradient.x) && same(b.gradient.y,a.gradient.y) && same(b.gradient.z,a.gradient.z),"CSE changed smooth value/gradient bits");
        check(ag.data()==storage,"CSE gradient workspace reallocated");
      }
      const auto tie=after.evaluateWithGradient({},ag);
      check(std::abs(tie.value+std::log(3.)/after.beta())<1e-14,"CSE removed duplicate occurrences in LSE parent");
    }
    auto compiled=std::make_shared<const HS::CompiledTree>(clones,leaves,true);
    HS::SmoothTree reused(compiled,{0,.001,true,true,true});
    check(&reused.compiledProgram()==compiled.get() && reused.nodeCount()==compiled->nodeCount(),
          "Smooth tree copied a shared compiled program");
    std::vector<HS::Vec3> batch_points;
    for(unsigned i=0;i<259;++i)
      batch_points.push_back({.01*static_cast<int>(i%23)-.11,.02*static_cast<int>(i%17)-.16,.03*static_cast<int>(i%11)-.15});
    batch_points.push_back({-0.,0.,std::numeric_limits<double>::denorm_min()});
    batch_points.push_back({1e4,-1e4,0}); // trusted exact-zero exp shortcut
    std::vector<double> expected_values(batch_points.size()),batch_values(batch_points.size());
    std::vector<HS::ValueGradient> expected_gradients(batch_points.size()),batch_gradients(batch_points.size()),trusted_gradients(batch_points.size());
    auto scalar_value_work=reused.makeValueWorkspace(); auto scalar_gradient_work=reused.makeGradientWorkspace();
    for(std::size_t i=0;i<batch_points.size();++i)
    {
      expected_values[i]=reused.evaluate(batch_points[i],scalar_value_work);
      expected_gradients[i]=reused.evaluateWithGradient(batch_points[i],scalar_gradient_work);
    }
    for(std::size_t workers_count:{1u,4u,8u})
    {
      auto value_batch=reused.makeValueBatchWorkspace(workers_count);
      auto gradient_batch=reused.makeGradientBatchWorkspace(workers_count);
      const auto value_bytes=value_batch.storageBytes(),gradient_bytes=gradient_batch.storageBytes();
      reused.evaluateBatch(batch_points,batch_values,value_batch);
      reused.evaluateWithGradientBatch(batch_points,batch_gradients,gradient_batch);
      reused.evaluateWithGradientTrustedBatch(batch_points,trusted_gradients,gradient_batch);
      for(std::size_t i=0;i<batch_points.size();++i)
      {
        check(same(batch_values[i],expected_values[i]) && same(batch_gradients[i].value,expected_gradients[i].value) &&
          same(batch_gradients[i].gradient.x,expected_gradients[i].gradient.x) &&
          same(batch_gradients[i].gradient.y,expected_gradients[i].gradient.y) &&
          same(batch_gradients[i].gradient.z,expected_gradients[i].gradient.z),"Smooth batch changed values/gradients/order");
        check(same(trusted_gradients[i].value,expected_gradients[i].value) &&
          same(trusted_gradients[i].gradient.x,expected_gradients[i].gradient.x) &&
          same(trusted_gradients[i].gradient.y,expected_gradients[i].gradient.y) &&
          same(trusted_gradients[i].gradient.z,expected_gradients[i].gradient.z),"Trusted smooth batch changed values/gradients/order");
      }
      reused.evaluateWithGradientTrustedBatch(batch_points,trusted_gradients,gradient_batch,true);
      for(std::size_t i=0;i<batch_points.size();++i)
        check(same(trusted_gradients[i].value,expected_gradients[i].value) &&
          same(trusted_gradients[i].gradient.x,expected_gradients[i].gradient.x) &&
          same(trusted_gradients[i].gradient.y,expected_gradients[i].gradient.y) &&
          same(trusted_gradients[i].gradient.z,expected_gradients[i].gradient.z),
          "Dynamic trusted batch changed values/gradients/order");
      reused.evaluateBatch(batch_points,batch_values,value_batch);
      reused.evaluateWithGradientBatch(batch_points,batch_gradients,gradient_batch);
      check(value_batch.storageBytes()==value_bytes && gradient_batch.storageBytes()==gradient_bytes,
            "Smooth batch workspace storage changed during reuse");
    }
    HS::SmoothTree paired_reused(compiled,{0,.001,true,true,true,true});
    for(std::size_t workers_count:{1u,4u,8u})
    {
      auto pair_work=paired_reused.makeGradientBatchWorkspace(workers_count);
      const auto bytes=pair_work.storageBytes();
      for(bool dynamic:{false,true})
      {
        paired_reused.evaluateWithGradientTrustedBatch(
            batch_points,trusted_gradients,pair_work,dynamic);
        for(std::size_t i=0;i<batch_points.size();++i)
          check(same(trusted_gradients[i].value,expected_gradients[i].value) &&
            same(trusted_gradients[i].gradient.x,expected_gradients[i].gradient.x) &&
            same(trusted_gradients[i].gradient.y,expected_gradients[i].gradient.y) &&
            same(trusted_gradients[i].gradient.z,expected_gradients[i].gradient.z),
            "Paired smooth batch changed values/gradients/order");
      }
      check(pair_work.storageBytes()==bytes,
            "Paired smooth batch allocated during reuse");
    }
    HS::BatchExecutor prewarmed_executor(4);
    check(prewarmed_executor.workerCount()==4,"Prewarmed batch executor worker count changed");
    auto prewarmed_work=prewarmed_executor.makeWorkspace<HS::ValueGradient>(
        paired_reused.nodeCount()+17,true);
    const auto prewarmed_bytes=prewarmed_work.storageBytes();
    for(unsigned rebound=0;rebound<2;++rebound)
    {
      paired_reused.rebindGradientBatchWorkspace(prewarmed_work);
      check(prewarmed_work.valuesPerWorker()==paired_reused.nodeCount(),
            "Rebound smooth batch workspace has wrong size");
      paired_reused.evaluateWithGradientTrustedBatch(
          batch_points,trusted_gradients,prewarmed_work,true);
      for(std::size_t i=0;i<batch_points.size();++i)
        check(same(trusted_gradients[i].value,expected_gradients[i].value) &&
          same(trusted_gradients[i].gradient.x,expected_gradients[i].gradient.x) &&
          same(trusted_gradients[i].gradient.y,expected_gradients[i].gradient.y) &&
          same(trusted_gradients[i].gradient.z,expected_gradients[i].gradient.z),
          "Prewarmed executor changed smooth batch results");
    }
    check(prewarmed_work.storageBytes()==prewarmed_bytes,
          "Rebinding a smaller smooth batch workspace allocated storage");
    auto empty_values=reused.makeValueBatchWorkspace(2); auto empty_gradients=reused.makeGradientBatchWorkspace(2);
    std::vector<HS::Vec3> no_points; std::vector<double> no_values; std::vector<HS::ValueGradient> no_gradients;
    reused.evaluateBatch(no_points,no_values,empty_values); reused.evaluateWithGradientBatch(no_points,no_gradients,empty_gradients);
    rejects([&] { reused.makeValueBatchWorkspace(0); }); rejects([&] { reused.makeGradientBatchWorkspace(0); });
    rejects([&] { HS::BatchExecutor invalid_executor(0); });
    auto moved_values=std::move(empty_values); auto moved_gradients=std::move(empty_gradients);
    rejects([&] { reused.evaluateBatch(batch_points,batch_values,empty_values); });
    rejects([&] { reused.evaluateWithGradientBatch(batch_points,batch_gradients,empty_gradients); });
    reused.evaluateBatch(batch_points,batch_values,moved_values);
    reused.evaluateWithGradientBatch(batch_points,batch_gradients,moved_gradients);
    std::vector<double> short_values(batch_points.size()-1); std::vector<HS::ValueGradient> short_gradients(batch_points.size()-1);
    rejects([&] { reused.evaluateBatch(batch_points,short_values,empty_values); });
    rejects([&] { reused.evaluateWithGradientBatch(batch_points,short_gradients,empty_gradients); });
    rejects([&] { reused.evaluateWithGradientTrustedBatch(batch_points,short_gradients,moved_gradients); });
    auto wrong_value_batch=duplicate.makeValueBatchWorkspace(2); auto wrong_gradient_batch=duplicate.makeGradientBatchWorkspace(2);
    rejects([&] { reused.evaluateBatch(batch_points,batch_values,wrong_value_batch); });
    rejects([&] { reused.evaluateWithGradientBatch(batch_points,batch_gradients,wrong_gradient_batch); });
    rejects([&] { reused.evaluateWithGradientTrustedBatch(batch_points,batch_gradients,wrong_gradient_batch); });
    auto invalid_batch_points=batch_points; invalid_batch_points.back().z=NAN;
    rejects([&] { reused.evaluateBatch(invalid_batch_points,batch_values,moved_values); });
    rejects([&] { reused.evaluateWithGradientBatch(invalid_batch_points,batch_gradients,moved_gradients); });
    reused.evaluateBatch(batch_points,batch_values,moved_values);
    reused.evaluateWithGradientBatch(batch_points,batch_gradients,moved_gradients);
    for(std::size_t i=0;i<batch_points.size();++i)
      check(same(batch_values[i],expected_values[i]) && same(batch_gradients[i].value,expected_gradients[i].value) &&
        same(batch_gradients[i].gradient.x,expected_gradients[i].gradient.x) &&
        same(batch_gradients[i].gradient.y,expected_gradients[i].gradient.y) &&
        same(batch_gradients[i].gradient.z,expected_gradients[i].gradient.z),"Smooth batch failed after worker exception");
    for(bool kernel:{false,true})
    {
      HS::SmoothTree small_field(clones,leaves,{0,.001,false,kernel,kernel});
      auto small_value_work=small_field.makeValueBatchWorkspace(8);
      auto small_gradient_work=small_field.makeGradientBatchWorkspace(8);
      auto vw=small_field.makeValueWorkspace(); auto gw=small_field.makeGradientWorkspace();
      for(std::size_t count:{1u,3u})
      {
        std::vector<HS::Vec3> small(batch_points.begin(),batch_points.begin()+count);
        std::vector<double> values(count); std::vector<HS::ValueGradient> gradients(count);
        small_field.evaluateBatch(small,values,small_value_work);
        small_field.evaluateWithGradientBatch(small,gradients,small_gradient_work);
        for(std::size_t i=0;i<count;++i)
        {
          const auto expected=small_field.evaluateWithGradient(small[i],gw);
          check(same(values[i],small_field.evaluate(small[i],vw)) && same(gradients[i].value,expected.value) &&
            same(gradients[i].gradient.x,expected.gradient.x) && same(gradients[i].gradient.y,expected.gradient.y) &&
            same(gradients[i].gradient.z,expected.gradient.z),"Smooth small batch lost points or changed kernel");
        }
      }
    }
    rejects([&] { HS::SmoothTree rejected(std::shared_ptr<const HS::CompiledTree>{},{0,.001,true}); });
    rejects([&] { HS::SmoothTree rejected(compiled,{0,.001,false}); });
    const auto retained=[&] {
      auto source_leaves=leaves;
      auto source=std::make_shared<const HS::CompiledTree>(clones,source_leaves,true);
      HS::SmoothTree field(source,{0,.001,true,true});
      source_leaves[0].offset=123;
      return field;
    }();
    auto retained_work=retained.makeGradientWorkspace(),reuse_work=reused.makeGradientWorkspace();
    auto copied=retained;
    check(&copied.compiledProgram()==&retained.compiledProgram(),"Smooth copy duplicated program storage");
    for(const HS::Vec3 p: {HS::Vec3{},HS::Vec3{-.01,.03,-.02},HS::Vec3{.03,-.02,.01}})
    {
      const auto actual=copied.evaluateWithGradient(p,retained_work),expected=reused.evaluateWithGradient(p,reuse_work);
      check(same(actual.value,expected.value) && same(actual.gradient.x,expected.gradient.x) &&
        same(actual.gradient.y,expected.gradient.y) && same(actual.gradient.z,expected.gradient.z),
        "Shared smooth program did not retain an independent coefficient snapshot");
    }
    // Specialization must match the generic evaluator BIT FOR BIT, including
    // gradients and ties. Exercise both winners, duplicates, CSE on/off and
    // fallback arities (unary/ternary/wide), not just the binary happy path.
    const auto compare_kernels=[&](const HS::ExpressionPtr& tree,HS::SmoothTreeOptions options) {
      options.binary_kernel=false; HS::SmoothTree generic(tree,leaves,options);
      options.binary_kernel=true; HS::SmoothTree binary(tree,leaves,options);
      options.predecode_kernel=true; HS::SmoothTree predecoded(tree,leaves,options);
      check(predecoded.predecodedKernelEnabled() && !binary.predecodedKernelEnabled(),
            "Predecoded kernel option was not applied");
      check(generic.nodeCount()==binary.nodeCount() && generic.workPerQuery()==binary.workPerQuery() &&
        generic.expCallsPerQuery()==binary.expCallsPerQuery() && generic.logCallsPerQuery()==binary.logCallsPerQuery(),"Kernel changed program/work counts");
      check(same(generic.beta(),binary.beta()) && same(generic.errorBound(),binary.errorBound()) &&
        same(generic.errorWeight(),binary.errorWeight()),"Kernel changed parameters/bounds");
      auto gv=generic.makeValueWorkspace(),bv=binary.makeValueWorkspace();
      auto gg=generic.makeGradientWorkspace(),bg=binary.makeGradientWorkspace();
      const auto storage_v=bv.data(); const auto storage_g=bg.data();
      const double tiny=std::numeric_limits<double>::denorm_min();
      std::vector<HS::Vec3> trusted_points;
      for(double x:{-0.,0.,tiny,-tiny,-.01,.01,1e150,-1e150})
        for(double y:{-0.,0.,tiny,-tiny,.01,-.01,1e150,-1e150})
        {
          const HS::Vec3 p{x,y,.03}; const auto g=generic.evaluateWithGradient(p,gg),b=binary.evaluateWithGradient(p,bg);
          trusted_points.push_back(p);
          check(same(generic.evaluate(p,gv),binary.evaluate(p,bv)) && same(g.value,b.value) &&
            same(g.gradient.x,b.gradient.x) && same(g.gradient.y,b.gradient.y) && same(g.gradient.z,b.gradient.z),"Binary kernel changed value/gradient bits");
        }
      auto trusted_work=predecoded.makeGradientBatchWorkspace(4);
      std::vector<HS::ValueGradient> trusted_values(trusted_points.size());
      predecoded.evaluateWithGradientTrustedBatch(trusted_points,trusted_values,trusted_work);
      options.paired_batch=true;
      HS::SmoothTree paired(tree,leaves,options);
      auto paired_work=paired.makeGradientBatchWorkspace(1);
      std::vector<HS::ValueGradient> paired_values(trusted_points.size());
      paired.evaluateWithGradientTrustedBatch(trusted_points,paired_values,paired_work);
      for(std::size_t i=0;i<trusted_points.size();++i)
      {
        const auto expected=binary.evaluateWithGradient(trusted_points[i],bg);
        const auto& actual=trusted_values[i];
        check(same(actual.value,expected.value) && same(actual.gradient.x,expected.gradient.x) &&
          same(actual.gradient.y,expected.gradient.y) && same(actual.gradient.z,expected.gradient.z),
          "Predecoded trusted kernel changed value/gradient bits");
        const auto& both=paired_values[i];
        check(same(both.value,actual.value) && same(both.gradient.x,actual.gradient.x) &&
          same(both.gradient.y,actual.gradient.y) && same(both.gradient.z,actual.gradient.z),
          "Paired trusted kernel changed value/gradient bits");
      }
      check(bv.data()==storage_v && bg.data()==storage_g,"Binary kernel allocated query workspace");
      auto wrong=bv; wrong.pop_back(); rejects([&] { binary.evaluate({},wrong); });
      rejects([&] { binary.evaluate({NAN,0,0},bv); });
      rejects([&] { binary.evaluateWithGradient({INFINITY,0,0},bg); });
    };
    for(auto op:{HS::TreeOperator::MIN,HS::TreeOperator::MAX})
      for(unsigned arity:{1u,2u,3u,7u}) for(bool sharing:{false,true})
      {
        std::vector<HS::ExpressionPtr> children; for(unsigned i=0;i<arity;++i) children.push_back(i%2 ? b : a);
        const auto tree=HS::makeNode(op,children);
        for(double beta:{.001,1000.,1e6,1e300}) compare_kernels(tree,{beta,0,sharing});
        compare_kernels(tree,{0,.001,sharing});
      }
    for(auto op:{HS::TreeOperator::MIN,HS::TreeOperator::MAX})
      compare_kernels(HS::makeNode(op,{a,a}),{1000,0,true});
    std::mt19937 rng(849321); std::vector<HS::ExpressionPtr> random_nodes{a,b,c};
    for(unsigned i=0;i<80;++i)
    {
      std::vector<HS::ExpressionPtr> children{random_nodes.back()};
      for(unsigned j=0,n=1+rng()%4;j<n;++j) children.push_back(random_nodes[rng()%random_nodes.size()]);
      random_nodes.push_back(HS::makeNode(i%2 ? HS::TreeOperator::MIN : HS::TreeOperator::MAX,children));
    }
    compare_kernels(random_nodes.back(),{0,.001,false}); compare_kernels(random_nodes.back(),{0,.001,true});
    // Non-axis gradient combinations exercise component cancellation too.
    leaves[0].normal={.25,-.125,.5}; leaves[1].normal={-.75,.125,-.25}; leaves[2].normal={.375,.625,.125};
    compare_kernels(random_nodes.back(),{1000,0,true});
    leaves[0].normal={2,0,0};
    rejects([&] { HS::SmoothTree invalid(random_nodes.back(),leaves,{1000,0,true,false,true}); });
    rejects([&] { HS::SmoothTree invalid(random_nodes.back(),leaves,{1000,0,true,true,false,true}); });
    for(bool kernel:{false,true})
    {
      HS::SmoothTree overflowing(HS::makeNode(HS::TreeOperator::MIN,{a,b}),leaves,{1000,0,true,kernel});
      auto work=overflowing.makeValueWorkspace();
      rejects([&] { overflowing.evaluate({std::numeric_limits<double>::max(),0,0},work); });
    }
    std::cout<<"Smooth tree: scalar/batch LSE and gradients, ties, DAG multiplicity, bounds and numeric limits passed\n";
  }
  catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
