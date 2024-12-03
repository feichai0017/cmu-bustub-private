#include "primer/hyperloglog.h"

namespace bustub {

template <typename KeyType>
HyperLogLog<KeyType>::HyperLogLog(int16_t n_bits)
    : n_bits_(n_bits), num_buckets_(1 << n_bits), registers_(num_buckets_, 0), cardinality_(0), mutex_() {
  if (n_bits < 0) {
    num_buckets_ = 0;
  }
  std::cout << "Initialized HLL with " << num_buckets_ << " buckets\n";
}

template <typename KeyType>
auto HyperLogLog<KeyType>::ComputeBinary(const hash_t &hash) const -> std::bitset<BITSET_CAPACITY> {
  auto binary = std::bitset<BITSET_CAPACITY>(hash);
  return binary;
}

template <typename KeyType>
auto HyperLogLog<KeyType>::PositionOfLeftmostOne(const std::bitset<BITSET_CAPACITY> &bset) const -> uint64_t {
  // 从最高位开始（跳过用于bucket的bits）
  int start = BITSET_CAPACITY - n_bits_ - 1;
  std::cout << "Starting position: " << start << std::endl;
  std::cout << "Looking at bits: ";
  for (int i = start; i >= 0; i--) {
    std::cout << bset[i];
  }
  std::cout << std::endl;

  // 从左向右扫描
  uint64_t pos = 1;
  for (int i = start; i >= 0; i--) {
    if (bset[i]) {
      pos = start - i + 1;
      std::cout << "Found 1 at position " << i << ", returning " << pos << std::endl;
      return pos;
    }
  }
  std::cout << "No 1 found, returning 1" << std::endl;
  return 1;
}

template <typename KeyType>
auto HyperLogLog<KeyType>::AddElem(KeyType val) -> void {
  std::cout << "\nAdding element: " << val << std::endl;

  // add write lock
  std::unique_lock lock(mutex_);

  hash_t calculate_hash = CalculateHash(val);
  std::cout << "Hash: " << calculate_hash << std::endl;

  std::bitset<BITSET_CAPACITY> compute_binary = ComputeBinary(calculate_hash);
  std::cout << "Binary: " << compute_binary << std::endl;

  size_t bucket_index = 0;
  std::cout << "First " << n_bits_ << " bits (for bucket index): ";
  for (int i = 0; i < n_bits_; i++) {
    if (compute_binary[BITSET_CAPACITY - 1 - i]) {
      bucket_index |= (1ULL << (n_bits_ - 1 - i));
      std::cout << "1";
    } else {
      std::cout << "0";
    }
  }
  std::cout << " (bucket: " << bucket_index << ")\n";

  uint64_t position = PositionOfLeftmostOne(compute_binary);
  std::cout << "Position of leftmost 1: " << position << std::endl;

  std::cout << "Old register value: " << (int)registers_[bucket_index] << std::endl;
  registers_[bucket_index] = std::max(registers_[bucket_index], static_cast<uint8_t>(position));
  std::cout << "New register value: " << (int)registers_[bucket_index] << std::endl;

  ComputeCardinality();
}

template <typename KeyType>
auto HyperLogLog<KeyType>::ComputeCardinality() -> void {
  if (num_buckets_ == 0) {
    cardinality_ = 0;
    return;
  }

  std::cout << "\nComputing cardinality:\n";
  double sum = 0.0;
  for (size_t i = 0; i < registers_.size(); i++) {
    double term = std::pow(2.0, -static_cast<double>(registers_[i]));
    sum += term;
    std::cout << "Register[" << i << "] = " << (int)registers_[i] << " -> 2^(-" << (int)registers_[i]
              << ") = " << std::fixed << std::setprecision(6) << term << std::endl;
  }
  std::cout << "Sum: " << sum << std::endl;

  if (sum == 0.0) {
    cardinality_ = 0;
    return;
  }

  double estimate = CONSTANT * num_buckets_ * num_buckets_ / sum;
  std::cout << "Estimate before floor: " << estimate << std::endl;
  cardinality_ = static_cast<size_t>(std::floor(estimate));
  std::cout << "Final cardinality: " << cardinality_ << "\n\n";
}

template class HyperLogLog<int64_t>;
template class HyperLogLog<std::string>;

}  // namespace bustub
