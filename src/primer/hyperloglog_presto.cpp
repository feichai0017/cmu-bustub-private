#include "primer/hyperloglog_presto.h"

#include <primer/hyperloglog.h>

namespace bustub {

template <typename KeyType>
HyperLogLogPresto<KeyType>::HyperLogLogPresto(int16_t n_leading_bits)
    : n_bits_(n_leading_bits),
      num_buckets_(1 << n_bits_),
      dense_bucket_(num_buckets_, std::bitset<4>(0)),
      cardinality_(0) {
  std::cout << "Initialized HLL Presto with " << num_buckets_ << " buckets\n";
}

template <typename KeyType>
auto HyperLogLogPresto<KeyType>::AddElem(KeyType val) -> void {
  // 1. 计算哈希值并转换为二进制
  hash_t hash = CalculateHash(val);

  // Get bucket index from high bits
  uint64_t bucket_idx = hash >> (64 - n_bits_);

  // Calculate leading zeros in remaining bits
  uint64_t remaining_hash = hash << n_bits_;  // Remove bucket bits
  uint64_t leading_zeros = 1;
  if (remaining_hash != 0) {
    leading_zeros += __builtin_clzll(remaining_hash) - n_bits_;
  }

  // Update buckets based on leading zeros
  if (leading_zeros <= DENSE_BUCKET_SIZE) {
    // Use dense bucket if value fits
    auto new_value = std::bitset<DENSE_BUCKET_SIZE>(leading_zeros);
    if (new_value.to_ullong() > dense_bucket_[bucket_idx].to_ullong()) {
      dense_bucket_[bucket_idx] = new_value;
    }
  } else {
    // Use overflow bucket for larger values
    dense_bucket_[bucket_idx] = std::bitset<DENSE_BUCKET_SIZE>((1 << DENSE_BUCKET_SIZE) - 1);
    auto overflow_value = leading_zeros - DENSE_BUCKET_SIZE;
    overflow_bucket_[bucket_idx] = std::bitset<OVERFLOW_BUCKET_SIZE>(overflow_value);
  }

  ComputeCardinality();
}

template <typename T>
auto HyperLogLogPresto<T>::ComputeCardinality() -> void {
  if (num_buckets_ == 0) {
    cardinality_ = 0;
    return;
  }

  // Calculate correction factor based on number of buckets
  double alpha;
  if (n_bits_ <= 4) {
    alpha = 0.673;
  } else if (n_bits_ <= 6) {
    alpha = 0.697;
  } else if (n_bits_ <= 8) {
    alpha = 0.709;
  } else {
    alpha = 0.7213 / (1.0 + 1.079 / num_buckets_);
  }

  // Calculate harmonic mean
  double sum = 0.0;
  size_t zeros = 0;

  for (size_t i = 0; i < num_buckets_; i++) {
    uint64_t rank = dense_bucket_[i].to_ullong();
    if (rank == 0) {
      zeros++;
    } else if (rank == (1ULL << DENSE_BUCKET_SIZE) - 1) {
      // Check overflow bucket
      auto it = overflow_bucket_.find(i);
      if (it != overflow_bucket_.end()) {
        rank = DENSE_BUCKET_SIZE + it->second.to_ullong();
      }
    }
    sum += std::pow(2.0, -static_cast<double>(rank));
  }

  double estimate = alpha * num_buckets_ * num_buckets_ / sum;

  // Apply corrections
  if (estimate <= 2.5 * num_buckets_) {
    if (zeros > 0) {
      estimate = num_buckets_ * std::log(static_cast<double>(num_buckets_) / zeros);
    }
  }

  cardinality_ = static_cast<uint64_t>(estimate);
}

template class HyperLogLogPresto<int64_t>;
template class HyperLogLogPresto<std::string>;
}  // namespace bustub
