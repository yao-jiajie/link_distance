#pragma once
#include <rokae_demo/halfspace_compiled_tree.hpp>
#include <cstdint>
#include <cstring>

namespace rokae_demo::halfspace
{
// Exact finite-coordinate keys: preserve signed zero and adjacent doubles.
// Hashes choose buckets; full bit equality decides reuse.
struct ValidationPointKey
{
  std::uint64_t bits[3]{};
  explicit ValidationPointKey(const Vec3& p)
  {
    static_assert(sizeof(double)==sizeof(std::uint64_t));
    if(!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
      throw std::invalid_argument("Validation cache query must be finite");
    std::memcpy(&bits[0],&p.x,sizeof(double));
    std::memcpy(&bits[1],&p.y,sizeof(double));
    std::memcpy(&bits[2],&p.z,sizeof(double));
  }
  bool operator==(const ValidationPointKey& other) const
  { return bits[0]==other.bits[0] && bits[1]==other.bits[1] && bits[2]==other.bits[2]; }
};
struct ValidationPointHash
{
  std::size_t operator()(const ValidationPointKey& key) const
  {
    std::size_t hash=0;
    for(const auto bits:key.bits)
      hash^=std::hash<std::uint64_t>{}(bits)+std::size_t{0x9e3779b9}+(hash<<6)+(hash>>2);
    return hash;
  }
};

// One validation run/worker only. Retain the immutable snapshot so identity
// cannot accidentally match a new program allocated at a recycled address.
// Only evaluation can populate this cache; callers cannot inject stale values.
class ValidationValueCache
{
  std::shared_ptr<const CompiledTree> program_;
  std::unordered_map<ValidationPointKey,double,ValidationPointHash> values_;
  std::size_t limit_,hits_=0,evaluations_=0;
public:
  explicit ValidationValueCache(std::shared_ptr<const CompiledTree> program,std::size_t limit=200000)
    : program_(std::move(program)),limit_(limit)
  { if(!program_) throw std::invalid_argument("Validation cache requires a program"); }
  const CompiledTree& program() const { return *program_; }
  std::size_t size() const { return values_.size(); }
  std::size_t hits() const { return hits_; }
  std::size_t evaluations() const { return evaluations_; }
  double evaluate(const Vec3& point,std::vector<double>& workspace)
  {
    if(workspace.size()!=program_->nodeCount()) throw std::invalid_argument("Wrong validation workspace size");
    const ValidationPointKey key(point);
    const auto found=values_.find(key);
    if(found!=values_.end()) { ++hits_; return found->second; }
    const auto value=program_->evaluate(point,workspace); ++evaluations_;
    // Failed/nonfinite queries must not become successful cache hits.
    if(std::isfinite(value) && values_.size()<limit_) values_.emplace(key,value);
    return value;
  }
};
} // namespace rokae_demo::halfspace
