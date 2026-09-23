#include <rokae_demo/boundary_direct_frame_cache.hpp>
#include <iostream>

using namespace rokae_demo;
void check(bool ok,const char* why) { if(!ok) throw std::runtime_error(why); }
template<class F> void rejects(F f,const char* why)
{ bool rejected=false; try { f(); } catch(const std::exception&) { rejected=true; } check(rejected,why); }

std::vector<VoxelPoint> shape(int extra=0)
{
  std::vector<VoxelPoint> points{{.25,.25,.25},{1.25,.25,.25},{.25,1.25,.25}};
  if(extra) points.emplace_back(2.25,.25,.25);
  return points;
}

BoundaryDirectFrameOptions options(double error=.01)
{
  BoundaryDirectFrameOptions result; result.smooth=true; result.smoothing={0,error,true,true};
  result.cache_validation_values=true; return result;
}

int main()
{
  try
  {
    BoundaryDirectFrameCache cache(2); auto first=cache.prepare(shape(),1,options());
    check(!first.cache_hit && cache.size()==1 && cache.misses()==1,"First frame was not a cache miss");
    check(first.snapshot->diagnostics().strict_sign_certified && first.snapshot->compiledProgram() &&
          first.snapshot->smoothTree() && first.snapshot->distanceReference() && first.snapshot->validationCache(),
          "Complete reusable snapshot was not prepared");
    const auto program=first.snapshot->compiledProgram(); const auto generation=first.snapshot->generation();
    auto* query_data=first.snapshot->queryWorkspace().data();
    const HS::Vec3 probe{-.5,-.5,-.5};
    const auto hard=program->evaluate(probe,first.snapshot->queryWorkspace());
    const auto smooth=first.snapshot->smoothTree()->evaluate(probe,first.snapshot->smoothWorkspace());
    (void)first.snapshot->validationCache()->evaluate(probe,first.snapshot->queryWorkspace());

    auto moved=shape(); for(auto& p:moved) p=VoxelPoint(p.x()+.5,p.y()+.5,p.z()+.5);
    std::reverse(moved.begin(),moved.end()); moved.push_back(moved.back());
    auto second=cache.prepare(moved,1,options());
    check(second.cache_hit && second.snapshot==first.snapshot && second.snapshot->generation()==generation,
          "Same occupied frame did not reuse its snapshot");
    check(second.snapshot->compiledProgram()==program && second.snapshot->queryWorkspace().data()==query_data,
          "Program or workspace was rebuilt on a hit");
    check(program->evaluate(probe,second.snapshot->queryWorkspace())==hard &&
          second.snapshot->smoothTree()->evaluate(probe,second.snapshot->smoothWorkspace())==smooth,
          "Reused program changed frame values");
    const auto validation_hits=second.snapshot->validationCache()->hits();
    (void)second.snapshot->validationCache()->evaluate(probe,second.snapshot->queryWorkspace());
    check(second.snapshot->validationCache()->hits()==validation_hits+1,
          "Snapshot-bound validation values were not retained");

    auto changed=cache.prepare(shape(1),1,options());
    check(!changed.cache_hit && changed.snapshot!=first.snapshot && cache.size()==2,
          "Changed occupancy reused stale geometry");
    auto first_again=cache.prepare(shape(),1,options());
    check(first_again.cache_hit && first_again.snapshot==first.snapshot,
          "LRU hit did not retain the requested snapshot");
    auto changed_config=options(); --changed_config.direct.max_depth;
    auto configured=cache.prepare(shape(),1,changed_config);
    check(!configured.cache_hit && cache.evictions()==1,"Changed build configuration did not miss/evict");
    check(cache.prepare(shape(),1,options()).cache_hit,
          "Recently used snapshot was evicted instead of the least-recent entry");
    auto changed_smooth=options(.02); auto smoothing_miss=cache.prepare(shape(),1,changed_smooth);
    check(!smoothing_miss.cache_hit,"Changed LSE configuration reused a program");
    auto evicted_geometry=cache.prepare(shape(1),1,options());
    check(!evicted_geometry.cache_hit,"Capacity limit did not evict least-recent geometry");

    const auto adjacent=std::nextafter(1.,2.); auto adjacent_h=cache.prepare(shape(),adjacent,options());
    check(!adjacent_h.cache_hit,"Adjacent voxel-size bit pattern reused a snapshot");
    cache.clear(); check(cache.size()==0,"Frame cache clear failed");
    check(program->evaluate(probe,first.snapshot->queryWorkspace())==hard,
          "Retained snapshot did not outlive cache eviction");
    rejects([] { BoundaryDirectFrameCache invalid(0); },"Zero cache capacity accepted");
    rejects([&] { cache.prepare({},1,options()); },"Empty frame accepted");
    auto mismatch=options(); mismatch.smoothing.share_subexpressions=false;
    rejects([&] { cache.prepare(shape(),1,mismatch); },"Mismatched smoothing/program configuration accepted");
    check(cache.hits()==3 && cache.misses()==6 && cache.evictions()==4,
          "Frame cache counters changed unexpectedly");

    BoundaryDirectFrameCache uncertified;
    auto capped=options(); capped.direct.signs.max_strata=1;
    const auto failed=uncertified.prepare(shape(),1,capped),failed_again=uncertified.prepare(shape(),1,capped);
    check(!failed.cache_hit && !failed.snapshot->diagnostics().applied && !failed.snapshot->compiledProgram() &&
          !failed_again.cache_hit && uncertified.size()==0 && uncertified.misses()==2,
          "Uncertified frame entered the reusable cache");
    std::cout<<"Boundary direct frame cache: geometry/config identity, LRU, workspaces, values and lifetime passed\n";
  }
  catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
