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
    // 计算哈希值
    hash_t hash = CalculateHash(val);

    // 获取桶索引（使用最高的 n_bits_ 位）
    size_t bucket_idx = hash >> (64 - n_bits_);

    // 计算最低有效位连续零的数量
    uint64_t remaining_bits = hash << n_bits_;  // 移除桶索引位
    uint32_t rho = 1;
    
    // 从最低位开始计算连续零的数量
    if (remaining_bits != 0) {
        rho += __builtin_clzll(remaining_bits);  // 使用前导零，因为我们已经左移了
        if (rho > TOTAL_BUCKET_SIZE) {
            rho = TOTAL_BUCKET_SIZE;
        }
    }

    // 更新桶
    if (rho <= DENSE_BUCKET_SIZE) {
        auto current_value = dense_bucket_[bucket_idx].to_ullong();
        if (rho > current_value) {
            dense_bucket_[bucket_idx] = std::bitset<DENSE_BUCKET_SIZE>(rho);
        }
    } else {
        dense_bucket_[bucket_idx].set();  // 设置为全1 (1111)
        overflow_bucket_[bucket_idx] = std::bitset<OVERFLOW_BUCKET_SIZE>(rho - DENSE_BUCKET_SIZE);
    }

    ComputeCardinality();
}

template <typename T>
auto HyperLogLogPresto<T>::ComputeCardinality() -> void {
    if (num_buckets_ == 0) {
        cardinality_ = 0;
        return;
    }
    
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

    double sum = 0.0;
    int zero_buckets = 0;
    
    for (size_t i = 0; i < num_buckets_; i++) {
        uint64_t value = dense_bucket_[i].to_ullong();
        if (value == 0) {
            zero_buckets++;
        } else if (value == (1ULL << DENSE_BUCKET_SIZE) - 1) {
            auto it = overflow_bucket_.find(i);
            if (it != overflow_bucket_.end()) {
                value = DENSE_BUCKET_SIZE + it->second.to_ullong();
            }
        }
        sum += std::pow(2.0, -static_cast<double>(value));
    }
    
    double estimate = alpha * num_buckets_ * num_buckets_ / sum;
    
    if (estimate <= 2.5 * num_buckets_ && zero_buckets > 0) {
        // 线性计数修正
        estimate = num_buckets_ * std::log(static_cast<double>(num_buckets_) / zero_buckets);
    }
    
    cardinality_ = static_cast<uint64_t>(std::floor(estimate));
}

template class HyperLogLogPresto<int64_t>;
template class HyperLogLogPresto<std::string>;
}  // namespace bustub
