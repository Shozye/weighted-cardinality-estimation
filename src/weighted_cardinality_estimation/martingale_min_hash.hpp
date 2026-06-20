#pragma once
#include <limits>
#include <string>
#include <vector>
#include <cmath>
#include <cstdint>
#include "sketch.hpp"
#include "hash_util.hpp"

// Unweighted MinHash with martingale estimator (Pettie, Wang & Yin 2020).
// Same register update as MinHash but non-mergeable and no Jaccard support.
// Estimator: E = Σ 1[state changed] / P_k, where P_k = 1 - exp(-sum(M))
// is the probability that at least one register updates on the k-th distinct element.
class MartingaleMinHash : public UnweightedSketch {
public:
    MartingaleMinHash(std::size_t sketch_size, std::uint64_t master_seed)
        : UnweightedSketch(sketch_size, master_seed),
          M_(sketch_size, std::numeric_limits<double>::infinity()),
          E_(0.0) {}

    MartingaleMinHash(std::size_t sketch_size, std::uint64_t master_seed,
                      const std::vector<double>& registers, double E)
        : UnweightedSketch(sketch_size, master_seed), M_(registers), E_(E) {}

    void add(const std::string& elem) override {
        // Compute P_k = probability of state change BEFORE updating registers.
        // P = 1 - prod_i exp(-M[i]) = 1 - exp(-sum(M[i]))
        // When registers are infinity (empty sketch), P = 1.
        double sum_M = 0.0;
        bool all_finite = true;
        for (std::size_t i = 0; i < size; ++i) {
            if (std::isinf(M_[i])) { all_finite = false; break; }
            sum_M += M_[i];
        }
        double P = all_finite ? (1.0 - std::exp(-sum_M)) : 1.0;

        // Update registers
        bool changed = false;
        for (std::size_t i = 0; i < size; ++i) {
            double g = -std::log(to_unit_interval(murmur64(elem, seeds_[i])));
            if (g < M_[i]) { M_[i] = g; changed = true; }
        }

        // Accumulate martingale estimate
        if (changed) E_ += 1.0 / P;
    }

    [[nodiscard]] double estimate() const override { return E_; }

    [[nodiscard]] double get_E() const { return E_; }
    [[nodiscard]] const std::vector<double>& get_registers() const { return M_; }

    [[nodiscard]] size_t memory_usage(uint64_t flags) const override {
        uint64_t f = resolve_flags(flags);
        size_t s = 0;
        if (f & MemoryFlag::REGISTERS) s += M_.capacity() * sizeof(double);
        if (f & MemoryFlag::ALL_WRITE_NO_REGISTERS) s += sizeof(size) + sizeof(E_);
        if (f & MemoryFlag::SEEDS) s += seeds_.bytes();
        return s;
    }

private:
    std::vector<double> M_;
    double E_;
};
