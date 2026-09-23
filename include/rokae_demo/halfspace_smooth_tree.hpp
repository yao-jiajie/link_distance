#pragma once
#include <rokae_demo/halfspace_compiled_tree.hpp>
#include <cstdint>
#include <type_traits>

namespace rokae_demo::halfspace
{
struct SmoothTreeOptions
{
  double beta=0;       // inverse meters; choose EXACTLY ONE parameter
  double max_error=0;  // meters of scalar smoothing error, NOT Euclidean SDF error
  bool share_subexpressions=false; // exact ordered CSE, not Boolean simplification
  bool binary_kernel=false; // specialize arity 2; retain generic arithmetic/order
  bool predecode_kernel=false; // optional trusted-batch instruction dispatch
  bool paired_batch=false; // optional two-point interleaving in trusted batches
};
struct ValueGradient { double value=0; Vec3 gradient; };
using GradientBatchWorkspace=EvaluationBatchWorkspace<ValueGradient>;

// For a certified boundary-direct q (negative free, positive occupied), returns
// d=-q_upper and its derivative. This is a distance PROXY, not an exact SDF.
// Upper MAX = log(sum exp(beta*a))/beta.
// Upper MIN = -log(mean exp(-beta*a))/beta.
// Composition gives q_upper>=q. Preserve EVERY child occurrence and operation;
// hard-tree absorption/idempotence/flattening is not valid after smoothing.
class SmoothTree
{
  std::shared_ptr<const CompiledTree> program_;
  double beta_=0,weight_=0,error_=0,lipschitz_=0,binary_min_zero_correction_=0;
  bool binary_kernel_=false,predecode_kernel_=false,paired_batch_=false;
  struct FastNode
  {
    enum class Kind : std::uint8_t { Leaf, BinaryMin, BinaryMax, TernaryMax, OtherMin, OtherMax };
    Kind kind;
    int first=0,second=0,third=0,begin=0,count=0;
  };
  std::vector<FastNode> fast_nodes_;
  // exp(-750) is below half the least positive binary64 subnormal, so a
  // correctly rounded binary64 exp is exactly +0. The checked API still calls
  // libm unconditionally; the trusted hot path uses this conservative cutoff.
  static constexpr double exact_exp_zero_cutoff_=750.0;
  template<bool Gradient> using Sample=std::conditional_t<Gradient,ValueGradient,double>;
  template<bool Gradient> static double scalar(const Sample<Gradient>& v)
  { if constexpr(Gradient) return v.value; else return v; }

  template<bool Gradient,bool BinaryKernel,bool Checked>
  double upperValue(const Vec3& p,std::vector<Sample<Gradient>>& work) const
  {
    if constexpr(Checked)
    {
      if(work.size()!=nodeCount()) throw std::invalid_argument("Wrong smooth workspace size");
      if(!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
        throw std::invalid_argument("Smooth query must be finite");
    }
    for(std::size_t id=0;id<nodeCount();++id)
    {
      const auto& n=program_->instructions_[id]; double result; Vec3 gradient;
      if(n.op==TreeOperator::LEAF)
      {
        const auto& leaf=program_->planes_[n.leaf];
        result=leaf.normal.x*p.x+leaf.normal.y*p.y+leaf.normal.z*p.z-leaf.offset;
        if constexpr(Gradient) gradient=leaf.normal;
      }
      else if(BinaryKernel && n.count==2)
      {
        // Same tie rule as the ordered generic scan: first occurrence wins.
        // No reciprocal, reassociation, abs-gap or softplus rewrite. The
        // trusted path only skips exp once its binary64 result is exactly zero.
        const bool is_max=n.op==TreeOperator::MAX;
        const int first=program_->children_[n.begin],second=program_->children_[n.begin+1];
        const double a=scalar<Gradient>(work[first]),b=scalar<Gradient>(work[second]);
        const bool choose_second=is_max ? b>a : b<a;
        const double anchor=choose_second ? b : a,other=choose_second ? a : b;
        const double gap=is_max ? anchor-other : other-anchor;
        const double scaled_gap=beta_*gap;
        double e;
        if constexpr(Checked) e=std::exp(-scaled_gap);
        else e=scaled_gap>exact_exp_zero_cutoff_ ? 0.0 : std::exp(-scaled_gap);
        double sum=0; sum+=e;
        if constexpr(Checked)
          if(!(sum>=0 && sum<=1)) throw std::runtime_error("Invalid shifted LSE sum");
        double correction;
        if constexpr(Checked)
          correction=(is_max ? std::log1p(sum) : std::log1p((1-sum)/(1+sum)))/beta_;
        else if(sum==0)
          correction=is_max ? 0.0 : binary_min_zero_correction_;
        else
          correction=(is_max ? std::log1p(sum) : std::log1p((1-sum)/(1+sum)))/beta_;
        result=anchor+correction;
        if constexpr(Checked)
          if(!(correction>=0 && result>=anchor)) throw std::runtime_error("Invalid LSE upper correction");
        if constexpr(Gradient)
        {
          gradient=work[choose_second ? second : first].gradient;
          const auto& g=work[choose_second ? first : second].gradient;
          gradient.x+=e*g.x; gradient.y+=e*g.y; gradient.z+=e*g.z;
          gradient.x/=1+sum; gradient.y/=1+sum; gradient.z/=1+sum;
        }
      }
      else
      {
        const bool is_max=n.op==TreeOperator::MAX;
        int chosen=n.begin;
        double anchor=scalar<Gradient>(work[program_->children_[chosen]]);
        for(int k=n.begin+1;k<n.begin+n.count;++k)
        {
          const double value=scalar<Gradient>(work[program_->children_[k]]);
          if(is_max ? value>anchor : value<anchor) { anchor=value; chosen=k; }
        }
        if constexpr(Gradient) gradient=work[program_->children_[chosen]].gradient;
        double sum=0;
        for(int k=n.begin;k<n.begin+n.count;++k) if(k!=chosen)
        {
          const auto child=program_->children_[k]; const double a=scalar<Gradient>(work[child]);
          const double gap=is_max ? anchor-a : a-anchor;
          // A positive overflow in beta*gap deliberately gives exp(-inf)=0.
          const double scaled_gap=beta_*gap;
          double e;
          if constexpr(Checked) e=std::exp(-scaled_gap);
          else e=scaled_gap>exact_exp_zero_cutoff_ ? 0.0 : std::exp(-scaled_gap);
          sum+=e;
          if constexpr(Gradient)
          { const auto& g=work[child].gradient; gradient.x+=e*g.x; gradient.y+=e*g.y; gradient.z+=e*g.z; }
        }
        if constexpr(Checked)
          if(!(sum>=0 && sum<=n.count-1)) throw std::runtime_error("Invalid shifted LSE sum");
        // MIN: log(k)-log(1+sum) = log1p((k-1-sum)/(1+sum)).
        // This avoids subtracting nearly equal logarithms at ties. Both node
        // corrections are nonnegative; never clamp the result against hard q.
        double correction;
        if constexpr(Checked)
          correction=(is_max ? std::log1p(sum) : std::log1p((n.count-1-sum)/(1+sum)))/beta_;
        else if(is_max && sum==0)
          correction=0.0;
        else
          correction=(is_max ? std::log1p(sum) : std::log1p((n.count-1-sum)/(1+sum)))/beta_;
        result=anchor+correction;
        if constexpr(Checked)
          if(!(correction>=0 && result>=anchor)) throw std::runtime_error("Invalid LSE upper correction");
        if constexpr(Gradient)
        { gradient.x/=1+sum; gradient.y/=1+sum; gradient.z/=1+sum; }
      }
      if constexpr(Checked)
        if(!std::isfinite(result)) throw std::overflow_error("Nonfinite smooth field value");
      if constexpr(Gradient)
      {
        if constexpr(Checked)
          if(!std::isfinite(gradient.x) || !std::isfinite(gradient.y) || !std::isfinite(gradient.z))
            throw std::overflow_error("Nonfinite smooth gradient");
        work[id]={result,gradient};
      }
      else work[id]=result;
    }
    return scalar<Gradient>(work[program_->root_]);
  }
  ValueGradient evaluateWithGradientTrusted(const Vec3& p,std::vector<ValueGradient>& work) const
  {
    const double q=predecode_kernel_ ? upperValueFast(p,work)
      : (binary_kernel_ ? upperValue<true,true,false>(p,work) : upperValue<true,false,false>(p,work));
    const auto& g=work[program_->root_].gradient;
    return {-q,{-g.x,-g.y,-g.z}};
  }
  // One predecoded DAG instruction. Keeping each point's arithmetic order
  // unchanged also permits interleaving independent points at each node.
  [[gnu::always_inline]] inline void fastNodeValueGradient(std::size_t id,const Vec3& p,
                             std::vector<ValueGradient>& work) const
  {
      const auto& n=fast_nodes_[id];
      double result;
      Vec3 gradient;
      if(n.kind==FastNode::Kind::Leaf)
      {
        const auto& leaf=program_->planes_[n.first];
        result=leaf.normal.x*p.x+leaf.normal.y*p.y+leaf.normal.z*p.z-leaf.offset;
        gradient=leaf.normal;
      }
      else if(n.kind==FastNode::Kind::BinaryMin || n.kind==FastNode::Kind::BinaryMax)
      {
        const bool is_max=n.kind==FastNode::Kind::BinaryMax;
        const int first=n.first,second=n.second;
        const double a=work[first].value,b=work[second].value;
        const bool choose_second=is_max ? b>a : b<a;
        const double anchor=choose_second ? b : a,other=choose_second ? a : b;
        const double gap=is_max ? anchor-other : other-anchor;
        const double scaled_gap=beta_*gap;
        const double e=scaled_gap>exact_exp_zero_cutoff_ ? 0.0 : std::exp(-scaled_gap);
        double sum=0; sum+=e;
        const double correction=sum==0 ? (is_max ? 0.0 : binary_min_zero_correction_)
          : (is_max ? std::log1p(sum) : std::log1p((1-sum)/(1+sum)))/beta_;
        result=anchor+correction;
        gradient=work[choose_second ? second : first].gradient;
        const auto& g=work[choose_second ? first : second].gradient;
        gradient.x+=e*g.x; gradient.y+=e*g.y; gradient.z+=e*g.z;
        gradient.x/=1+sum; gradient.y/=1+sum; gradient.z/=1+sum;
      }
      else
      {
        const bool is_max=n.kind==FastNode::Kind::TernaryMax || n.kind==FastNode::Kind::OtherMax;
        int chosen=n.begin;
        double anchor;
        if(n.kind==FastNode::Kind::TernaryMax)
        {
          chosen=0;
          anchor=work[n.first].value;
          const double second=work[n.second].value;
          if(second>anchor) { anchor=second; chosen=1; }
          const double third=work[n.third].value;
          if(third>anchor) { anchor=third; chosen=2; }
          const int children[3]{n.first,n.second,n.third};
          gradient=work[children[chosen]].gradient;
          double sum=0;
          for(int k=0;k<3;++k) if(k!=chosen)
          {
            const auto child=children[k]; const double a=work[child].value;
            const double gap=anchor-a;
            const double scaled_gap=beta_*gap;
            const double e=scaled_gap>exact_exp_zero_cutoff_ ? 0.0 : std::exp(-scaled_gap);
            sum+=e;
            const auto& g=work[child].gradient;
            gradient.x+=e*g.x; gradient.y+=e*g.y; gradient.z+=e*g.z;
          }
          const double correction=sum==0 ? 0.0 : std::log1p(sum)/beta_;
          result=anchor+correction;
          gradient.x/=1+sum; gradient.y/=1+sum; gradient.z/=1+sum;
        }
        else
        {
          anchor=work[program_->children_[chosen]].value;
          for(int k=n.begin+1;k<n.begin+n.count;++k)
          {
            const double value=work[program_->children_[k]].value;
            if(is_max ? value>anchor : value<anchor) { anchor=value; chosen=k; }
          }
          gradient=work[program_->children_[chosen]].gradient;
          double sum=0;
          for(int k=n.begin;k<n.begin+n.count;++k) if(k!=chosen)
          {
            const auto child=program_->children_[k]; const double a=work[child].value;
            const double gap=is_max ? anchor-a : a-anchor;
            const double scaled_gap=beta_*gap;
            const double e=scaled_gap>exact_exp_zero_cutoff_ ? 0.0 : std::exp(-scaled_gap);
            sum+=e;
            const auto& g=work[child].gradient;
            gradient.x+=e*g.x; gradient.y+=e*g.y; gradient.z+=e*g.z;
          }
          const double correction=is_max && sum==0 ? 0.0
            : (is_max ? std::log1p(sum) : std::log1p((n.count-1-sum)/(1+sum)))/beta_;
          result=anchor+correction;
          gradient.x/=1+sum; gradient.y/=1+sum; gradient.z/=1+sum;
        }
      }
      work[id]={result,gradient};
  }
  // Predecoded hot kernel for the common binary/ternary node shapes. It uses
  // exactly the same ordered arithmetic as upperValue<true,true,false>.
  double upperValueFast(const Vec3& p,std::vector<ValueGradient>& work) const
  {
    for(std::size_t id=0;id<fast_nodes_.size();++id)
      fastNodeValueGradient(id,p,work);
    return work[program_->root_].value;
  }
  void evaluateWithGradientTrustedPair(
      const Vec3& a,const Vec3& b,
      std::vector<ValueGradient>& first,std::vector<ValueGradient>& second,
      ValueGradient& first_result,ValueGradient& second_result) const
  {
    for(std::size_t id=0;id<fast_nodes_.size();++id)
    {
      fastNodeValueGradient(id,a,first);
      fastNodeValueGradient(id,b,second);
    }
    const auto& f=first[program_->root_];
    const auto& s=second[program_->root_];
    first_result={-f.value,{-f.gradient.x,-f.gradient.y,-f.gradient.z}};
    second_result={-s.value,{-s.gradient.x,-s.gradient.y,-s.gradient.z}};
  }
public:
  // The program must remain immutable through all aliases for this field's
  // lifetime; construct it with make_shared<const CompiledTree>.
  SmoothTree(std::shared_ptr<const CompiledTree> program,const SmoothTreeOptions& options)
    : program_(std::move(program)),binary_kernel_(options.binary_kernel),
      predecode_kernel_(options.predecode_kernel),paired_batch_(options.paired_batch)
  {
    if(!program_)
      throw std::invalid_argument("Smooth tree requires a compiled program");
    if(program_->structuralSharingEnabled()!=options.share_subexpressions)
      throw std::invalid_argument("Smooth options must match compiled structural sharing");
    if(predecode_kernel_ && !binary_kernel_)
      throw std::invalid_argument("Predecoded smooth kernel requires binary kernel");
    if(paired_batch_ && !predecode_kernel_)
      throw std::invalid_argument("Paired smooth batch requires predecoded kernel");
    if(!std::isfinite(options.beta) || !std::isfinite(options.max_error) || options.beta<0 || options.max_error<0 ||
       (options.beta>0)==(options.max_error>0)) throw std::invalid_argument("Choose positive finite beta OR scalar smoothing error");
    std::vector<long double> weights(nodeCount(),0);
    for(std::size_t i=0;i<nodeCount();++i)
    {
      const auto& n=program_->instructions_[i];
      if(n.op==TreeOperator::LEAF)
      {
        const auto& normal=program_->planes_[n.leaf].normal;
        lipschitz_=std::max(lipschitz_,std::hypot(normal.x,normal.y,normal.z));
      }
      else
      {
        for(int k=n.begin;k<n.begin+n.count;++k) weights[i]=std::max(weights[i],weights[program_->children_[k]]);
        weights[i]+=std::log(static_cast<long double>(n.count));
      }
    }
    weight_=static_cast<double>(weights.back());
    beta_=options.beta>0 ? options.beta : std::nextafter(static_cast<double>((weight_>0 ? weights.back() : 1)/options.max_error),INFINITY);
    error_=std::nextafter(static_cast<double>(weights.back()/beta_),INFINITY);
    if(weight_==0) error_=0;
    if(!(beta_>0) || !std::isfinite(beta_) || !std::isfinite(1/beta_) || !std::isfinite(error_) || !std::isfinite(lipschitz_))
      throw std::invalid_argument("Unrepresentable smooth parameter/error/normal scale");
    binary_min_zero_correction_=std::log1p(1.0)/beta_;
    if(predecode_kernel_)
    {
      fast_nodes_.reserve(nodeCount());
      for(const auto& n:program_->instructions_)
      {
        FastNode fast{};
        if(n.op==TreeOperator::LEAF)
        { fast.kind=FastNode::Kind::Leaf; fast.first=n.leaf; }
        else if(n.count==2)
        {
          fast.kind=n.op==TreeOperator::MAX ? FastNode::Kind::BinaryMax : FastNode::Kind::BinaryMin;
          fast.first=program_->children_[n.begin]; fast.second=program_->children_[n.begin+1];
        }
        else if(n.op==TreeOperator::MAX && n.count==3)
        {
          fast.kind=FastNode::Kind::TernaryMax;
          fast.first=program_->children_[n.begin]; fast.second=program_->children_[n.begin+1];
          fast.third=program_->children_[n.begin+2];
        }
        else
        {
          fast.kind=n.op==TreeOperator::MAX ? FastNode::Kind::OtherMax : FastNode::Kind::OtherMin;
          fast.begin=n.begin; fast.count=n.count;
        }
        fast_nodes_.push_back(fast);
      }
    }
  }
  SmoothTree(const SerializedTree& tree,const std::vector<HalfspaceLeaf>& leaves,const SmoothTreeOptions& options)
    : SmoothTree(std::make_shared<const CompiledTree>(tree,leaves,options.share_subexpressions),options) {}
  SmoothTree(const ExpressionPtr& tree,const std::vector<HalfspaceLeaf>& leaves,const SmoothTreeOptions& options)
    : SmoothTree(std::make_shared<const CompiledTree>(tree,leaves,options.share_subexpressions),options) {}
  const CompiledTree& compiledProgram() const { return *program_; }
  std::size_t nodeCount() const { return program_->nodeCount(); }
  std::size_t sourceNodeCount() const { return program_->sourceNodeCount(); }
  std::size_t sharedNodeCount() const { return program_->sharedNodeCount(); }
  std::size_t edgeCount() const { return program_->edgeCount(); }
  std::size_t logCallsPerQuery() const
  { return std::count_if(program_->instructions_.begin(),program_->instructions_.end(),[](const auto& n) { return n.op!=TreeOperator::LEAF; }); }
  std::size_t expCallsPerQuery() const { return edgeCount()-logCallsPerQuery(); }
  bool binaryKernelEnabled() const { return binary_kernel_; }
  bool predecodedKernelEnabled() const { return predecode_kernel_; }
  bool pairedBatchEnabled() const { return paired_batch_; }
  std::size_t binaryNodeCount() const
  { return binary_kernel_ ? std::count_if(program_->instructions_.begin(),program_->instructions_.end(),[](const auto& n) { return n.op!=TreeOperator::LEAF && n.count==2; }) : 0; }
  std::size_t workPerQuery() const { return nodeCount()+3*program_->edgeCount(); }
  // Inclusive of the shared program; do not sum this with its programBytes().
  // Allocator/control-block overhead is not included.
  std::size_t programBytes() const
  { return sizeof(*this)+program_->programBytes()+fast_nodes_.capacity()*sizeof(FastNode); }
  double beta() const { return beta_; }
  double errorWeight() const { return weight_; }
  double errorBound() const { return error_; } // path-sum analytic bound; not interval-certified libm
  double lipschitzBound() const { return lipschitz_; }
  std::vector<double> makeValueWorkspace() const { return std::vector<double>(nodeCount()); }
  std::vector<ValueGradient> makeGradientWorkspace() const { return std::vector<ValueGradient>(nodeCount()); }
  ValueBatchWorkspace makeValueBatchWorkspace(std::size_t worker_count) const
  { return ValueBatchWorkspace(worker_count,nodeCount()); }
  void rebindGradientBatchWorkspace(GradientBatchWorkspace& workspace) const
  {
    if(workspace.workers_.empty() || !workspace.pool_)
      throw std::invalid_argument("Cannot rebind an empty smooth gradient batch workspace");
    for(auto& worker:workspace.workers_)
      worker.resize(nodeCount());
    if(paired_batch_)
    {
      if(workspace.paired_workers_.size()!=workspace.workers_.size())
        workspace.paired_workers_.resize(workspace.workers_.size());
      for(auto& worker:workspace.paired_workers_)
        worker.resize(nodeCount());
    }
    else workspace.paired_workers_.clear();
    workspace.values_per_worker_=nodeCount();
  }
  GradientBatchWorkspace makeGradientBatchWorkspace(std::size_t worker_count) const
  { return GradientBatchWorkspace(worker_count,nodeCount(),paired_batch_); }
  GradientBatchWorkspace makeGradientBatchWorkspace(const BatchExecutor& executor) const
  { return GradientBatchWorkspace(executor,nodeCount(),paired_batch_); }
  double evaluate(const Vec3& p,std::vector<double>& work) const
  { return -(binary_kernel_ ? upperValue<false,true,true>(p,work) : upperValue<false,false,true>(p,work)); }
  ValueGradient evaluateWithGradient(const Vec3& p,std::vector<ValueGradient>& work) const
  {
    const double q=binary_kernel_ ? upperValue<true,true,true>(p,work) : upperValue<true,false,true>(p,work);
    const auto& g=work[program_->root_].gradient;
    return {-q,{-g.x,-g.y,-g.z}};
  }
  // Ordered batches with private worker buffers, as for CompiledTree. Every
  // task finishes before return/rethrow; reusable workspace threads stay alive.
  // On failure outputs may be partial.
  void evaluateBatch(const std::vector<Vec3>& points,std::vector<double>& values,
                     ValueBatchWorkspace& workspace) const
  {
    if(values.size()!=points.size()) throw std::invalid_argument("Wrong smooth batch output size");
    if(workspace.workers_.empty() || !workspace.pool_ ||
       workspace.values_per_worker_!=nodeCount())
      throw std::invalid_argument("Wrong smooth batch workspace size");
    detail::evaluatePointBatch(points.size(),*workspace.pool_,[&](std::size_t worker,std::size_t begin,std::size_t end) {
      auto& work=workspace.workers_[worker];
      for(auto i=begin;i<end;++i) values[i]=evaluate(points[i],work);
    });
  }
  void evaluateWithGradientBatch(const std::vector<Vec3>& points,std::vector<ValueGradient>& values,
                                 GradientBatchWorkspace& workspace) const
  {
    if(values.size()!=points.size()) throw std::invalid_argument("Wrong smooth gradient batch output size");
    if(workspace.workers_.empty() || !workspace.pool_ ||
       workspace.values_per_worker_!=nodeCount())
      throw std::invalid_argument("Wrong smooth gradient batch workspace size");
    detail::evaluatePointBatch(points.size(),*workspace.pool_,[&](std::size_t worker,std::size_t begin,std::size_t end) {
      auto& work=workspace.workers_[worker];
      for(auto i=begin;i<end;++i) values[i]=evaluateWithGradient(points[i],work);
    });
  }
  // Hot path for callers that already validated all input points and validate
  // all returned samples. It retains the exact arithmetic/order of the checked
  // API while omitting redundant per-node finite/range tests.
  void evaluateWithGradientTrustedBatch(const std::vector<Vec3>& points,std::vector<ValueGradient>& values,
                                        GradientBatchWorkspace& workspace,bool dynamic_schedule=false) const
  {
    if(values.size()!=points.size()) throw std::invalid_argument("Wrong smooth gradient batch output size");
    if(workspace.workers_.empty() || !workspace.pool_ ||
       workspace.values_per_worker_!=nodeCount() ||
       (paired_batch_ &&
        (workspace.paired_workers_.size()!=workspace.workers_.size() ||
         workspace.paired_workers_[0].size()!=nodeCount())))
      throw std::invalid_argument("Wrong smooth gradient batch workspace size");
    const auto evaluate_range=[&](std::size_t worker,std::size_t begin,std::size_t end) {
      auto& work=workspace.workers_[worker];
      if(paired_batch_)
      {
        auto& second=workspace.paired_workers_[worker];
        auto i=begin;
        for(;i+1<end;i+=2)
          evaluateWithGradientTrustedPair(points[i],points[i+1],work,second,
                                         values[i],values[i+1]);
        if(i<end) values[i]=evaluateWithGradientTrusted(points[i],work);
      }
      else for(auto i=begin;i<end;++i) values[i]=evaluateWithGradientTrusted(points[i],work);
    };
    if(dynamic_schedule)
      detail::evaluatePointBatchDynamic(points.size(),*workspace.pool_,evaluate_range);
    else
      detail::evaluatePointBatch(points.size(),*workspace.pool_,evaluate_range);
  }
};
} // namespace rokae_demo::halfspace
