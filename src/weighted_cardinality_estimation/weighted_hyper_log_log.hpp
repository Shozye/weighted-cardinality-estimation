#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include "sketch.hpp"
#include "hash_util.hpp"

// Weighted HyperLogLog based on Cohen, Katzir & Yehezkel (IPL 2015).
// Uses ExpSketch-style registers (min of -log(u)/w) with stochastic averaging:
//   H1(x) ~ {1..m} selects bucket, H2(x) ~ U(0,1) → g = -log(u)/w.
//   M_k = min{g(x_i) | H1(x_i) = k}.
//   Estimator: (K - 1) / sum(M_k for non-inf k) with Poisson correction.
template <typename T = double>
class WeightedHyperLogLogT : public Sketch, public MergeableMixin {
public:
    WeightedHyperLogLogT(std::size_t sketch_size, std::uint64_t master_seed)
        : Sketch(sketch_size, master_seed),
          M_(sketch_size, std::numeric_limits<T>::infinity()),
          seed1_(seeds_[0]),
          seed2_(seeds_[1]) {}

    WeightedHyperLogLogT(std::size_t sketch_size, std::uint64_t master_seed,
                        const std::vector<T>& registers)
        : Sketch(sketch_size, master_seed), M_(registers),
          seed1_(seeds_[0]),
          seed2_(seeds_[1]) {}

    void add(const std::string& elem, double weight = 1.0) override {
        validate_weight(weight);
        std::uint64_t h1 = murmur64(elem, seed1_);
        std::size_t k = h1 % size;
        std::uint64_t h2 = murmur64(elem, seed2_);
        double u = to_unit_interval(h2);
        T g = static_cast<T>(-std::log(u) / weight);
        if (g < M_[k]) M_[k] = g;
    }

    [[nodiscard]] double estimate() const override {
        std::size_t K = 0;
        double sum = 0.0;
        for (std::size_t i = 0; i < size; ++i) {
            if (M_[i] < std::numeric_limits<T>::infinity()) {
                sum += static_cast<double>(M_[i]);
                ++K;
            }
        }
        if (K == 0) return 0.0;
        if (sum <= 0.0) return std::numeric_limits<double>::infinity();

        double m = static_cast<double>(size);
        double Kd = static_cast<double>(K);
        double correction = (K > 1) ? (Kd - 1.0) / Kd : 0.5;
        return correction * Kd * m / sum;
    }

    void merge(const WeightedHyperLogLogT& other) {
        if (other.size != size)
            throw std::invalid_argument("Cannot merge sketches of different sizes.");
        for (std::size_t i = 0; i < size; ++i)
            M_[i] = std::min(M_[i], other.M_[i]);
    }

    [[nodiscard]] std::vector<T> get_registers() const { return M_; }
    [[nodiscard]] size_t memory_usage(uint64_t flags) const override {
        uint64_t f = resolve_flags(flags);
        size_t s = 0;
        if (f & MemoryFlag::REGISTERS) s += M_.capacity() * sizeof(T);
        if (f & MemoryFlag::ALL_WRITE_NO_REGISTERS) s += sizeof(size);
        if (f & MemoryFlag::SEEDS) s += seeds_.bytes();
        return s;
    }

private:
    std::vector<T> M_;
    std::uint32_t seed1_;
    std::uint32_t seed2_;
};

using WeightedHyperLogLog = WeightedHyperLogLogT<double>;
using WeightedHyperLogLogFloat32 = WeightedHyperLogLogT<float>;
