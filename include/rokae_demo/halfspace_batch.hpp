#pragma once

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <condition_variable>
#include <exception>
#include <functional>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

namespace rokae_demo::halfspace
{
class CompiledTree;
class SmoothTree;

namespace detail
{
// Reusable caller-owned workers. Thread creation belongs to workspace setup,
// not to every batch query. One batch at a time may use a workspace.
class BatchWorkerPool
{
  std::size_t worker_count_=0;
  std::vector<std::thread> threads_;
  std::vector<std::exception_ptr> errors_;
  std::mutex mutex_;
  std::condition_variable start_;
  std::condition_variable done_;
  std::function<void(std::size_t)> task_;
  std::size_t generation_=0,completed_=0;
  bool stopping_=false;

  void worker(std::size_t id)
  {
    std::size_t observed=0;
    for(;;)
    {
      std::unique_lock<std::mutex> lock(mutex_);
      start_.wait(lock,[&] { return stopping_ || generation_!=observed; });
      if(stopping_) return;
      observed=generation_;
      const auto* task=&task_;
      lock.unlock();
      std::exception_ptr error;
      try { (*task)(id); }
      catch(...) { error=std::current_exception(); }
      lock.lock();
      errors_[id]=error;
      ++completed_;
      if(completed_==threads_.size()) done_.notify_one();
    }
  }

public:
  explicit BatchWorkerPool(std::size_t worker_count)
    : worker_count_(worker_count),errors_(worker_count)
  {
    if(!worker_count_) throw std::invalid_argument("Batch worker pool requires at least one worker");
    threads_.reserve(worker_count_-1);
    try
    {
      for(std::size_t id=1;id<worker_count_;++id)
        threads_.emplace_back([this,id] { worker(id); });
    }
    catch(...)
    {
      {
        std::lock_guard<std::mutex> lock(mutex_);
        stopping_=true;
      }
      start_.notify_all();
      for(auto& thread:threads_) if(thread.joinable()) thread.join();
      throw;
    }
  }
  BatchWorkerPool(const BatchWorkerPool&)=delete;
  BatchWorkerPool& operator=(const BatchWorkerPool&)=delete;
  ~BatchWorkerPool()
  {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      stopping_=true;
    }
    start_.notify_all();
    for(auto& thread:threads_) thread.join();
  }
  std::size_t workerCount() const { return worker_count_; }
  std::size_t storageBytes() const
  { return sizeof(*this)+threads_.capacity()*sizeof(std::thread)+errors_.capacity()*sizeof(std::exception_ptr); }

  template<class Function> void run(Function&& function)
  {
    if(worker_count_==1)
    {
      function(0);
      return;
    }
    {
      std::lock_guard<std::mutex> lock(mutex_);
      std::fill(errors_.begin(),errors_.end(),std::exception_ptr{});
      task_=std::forward<Function>(function);
      completed_=0;
      ++generation_;
    }
    start_.notify_all();
    std::exception_ptr caller_error;
    try { task_(0); }
    catch(...) { caller_error=std::current_exception(); }
    std::exception_ptr first_error=caller_error;
    {
      std::unique_lock<std::mutex> lock(mutex_);
      errors_[0]=caller_error;
      done_.wait(lock,[&] { return completed_==threads_.size(); });
      for(const auto& error:errors_) if(!first_error && error) first_error=error;
      task_={};
    }
    if(first_error) std::rethrow_exception(first_error);
  }
};
} // namespace detail

template<class Sample> class EvaluationBatchWorkspace;

// Field-independent worker threads can be created before a changing geometry
// is available, then shared sequentially by workspaces rebound to each field.
// Workspaces using the same executor must never evaluate concurrently.
class BatchExecutor
{
  template<class Sample> friend class EvaluationBatchWorkspace;
  std::shared_ptr<detail::BatchWorkerPool> pool_;
public:
  explicit BatchExecutor(std::size_t worker_count)
    : pool_(std::make_shared<detail::BatchWorkerPool>(worker_count)) {}
  std::size_t workerCount() const { return pool_->workerCount(); }
  template<class Sample>
  EvaluationBatchWorkspace<Sample> makeWorkspace(
      std::size_t values_per_worker,bool pair_workspace=false) const;
};

// One scalar evaluator workspace per worker. Construct through the matching
// tree factory so a batch cannot silently use the wrong node count.
// Reuse sequentially; overlapping batches need separate workspaces and outputs.
// Buffers contain no cached results: another tree of the same size may reuse them.
template<class Sample> class EvaluationBatchWorkspace
{
  friend class BatchExecutor;
  friend class CompiledTree;
  friend class SmoothTree;
  std::vector<std::vector<Sample>> workers_;
  std::vector<std::vector<Sample>> paired_workers_;
  std::shared_ptr<detail::BatchWorkerPool> pool_;
  std::size_t values_per_worker_=0;

  EvaluationBatchWorkspace(std::size_t worker_count,std::size_t values_per_worker,
                           bool pair_workspace=false)
    : workers_(worker_count,std::vector<Sample>(values_per_worker)),
      paired_workers_(pair_workspace ? worker_count : 0,
                      std::vector<Sample>(pair_workspace ? values_per_worker : 0)),
      pool_(std::make_shared<detail::BatchWorkerPool>(worker_count)),
      values_per_worker_(values_per_worker)
  {
    if(!worker_count) throw std::invalid_argument("Batch workspace requires at least one worker");
  }
  EvaluationBatchWorkspace(const BatchExecutor& executor,std::size_t values_per_worker,
                           bool pair_workspace=false)
    : workers_(executor.workerCount(),std::vector<Sample>(values_per_worker)),
      paired_workers_(pair_workspace ? executor.workerCount() : 0,
                      std::vector<Sample>(pair_workspace ? values_per_worker : 0)),
      pool_(executor.pool_),values_per_worker_(values_per_worker)
  {}

public:
  EvaluationBatchWorkspace(const EvaluationBatchWorkspace&)=delete;
  EvaluationBatchWorkspace& operator=(const EvaluationBatchWorkspace&)=delete;
  EvaluationBatchWorkspace(EvaluationBatchWorkspace&&) noexcept=default;
  EvaluationBatchWorkspace& operator=(EvaluationBatchWorkspace&&) noexcept=default;
  std::size_t workerCount() const { return workers_.size(); }
  std::size_t valuesPerWorker() const { return values_per_worker_; }
  std::size_t storageBytes() const
  {
    std::size_t result=sizeof(*this)+
      (workers_.capacity()+paired_workers_.capacity())*sizeof(std::vector<Sample>);
    for(const auto& worker:workers_) result+=worker.capacity()*sizeof(Sample);
    for(const auto& worker:paired_workers_) result+=worker.capacity()*sizeof(Sample);
    if(pool_) result+=pool_->storageBytes();
    return result;
  }
};

template<class Sample>
EvaluationBatchWorkspace<Sample> BatchExecutor::makeWorkspace(
    std::size_t values_per_worker,bool pair_workspace) const
{
  return EvaluationBatchWorkspace<Sample>(
      *this,values_per_worker,pair_workspace);
}

using ValueBatchWorkspace=EvaluationBatchWorkspace<double>;

namespace detail
{
template<class Function>
void evaluatePointBatch(std::size_t point_count,BatchWorkerPool& pool,Function&& function)
{
  if(!point_count) return;
  const auto worker_count=std::min(point_count,pool.workerCount());
  pool.run([&](std::size_t worker) {
    if(worker>=worker_count) return;
    const auto base=point_count/worker_count,remainder=point_count%worker_count;
    const auto begin=worker*base+std::min(worker,remainder);
    const auto end=begin+base+(worker<remainder);
    function(worker,begin,end);
  });
}

// Distribute small groups of independent points to reduce tail latency when
// their libm work differs. Each worker keeps its private evaluation buffer.
template<class Function>
void evaluatePointBatchDynamic(std::size_t point_count,BatchWorkerPool& pool,Function&& function)
{
  if(!point_count) return;
  std::atomic<std::size_t> next{0};
  constexpr std::size_t chunk=2;
  pool.run([&](std::size_t worker) {
    for(;;)
    {
      const auto begin=next.fetch_add(chunk,std::memory_order_relaxed);
      if(begin>=point_count) return;
      function(worker,begin,std::min(begin+chunk,point_count));
    }
  });
}
} // namespace detail
} // namespace rokae_demo::halfspace
