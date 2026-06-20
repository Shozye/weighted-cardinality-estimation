#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>
#include "sketch.hpp"
#include "hash_util.hpp"

// HyperLogLog (Flajolet et al. 2007).
// Uses a single hash per element: first log2(m) bits select the register,
// remaining bits compute rho (position of leftmost 1-bit).
// Estimator: alpha_m * m^2 / sum(2^{-M[i]}).
// Supports merge (element-wise max). Unweighted, no Jaccard.
class HyperLogLog : public UnweightedSketch, public MergeableMixin {
public:
    HyperLogLog(std::size_t sketch_size, std::uint64_t master_seed)
        : UnweightedSketch(sketch_size, master_seed), M_(sketch_size, 0),
          seed_(seeds_[0]) {}

    HyperLogLog(std::size_t sketch_size, std::uint64_t master_seed,
                const std::vector<uint8_t>& registers)
        : UnweightedSketch(sketch_size, master_seed), M_(registers),
          seed_(seeds_[0]) {}

    void add(const std::string& elem) override {
        std::uint64_t h = murmur64(elem, seed_);
        // Bucket selection: use upper bits modulo m for non-power-of-2 support
        std::size_t j = h % size;
        // Remaining bits → rho: rehash with different seed for independent bits
        std::uint64_t h2 = murmur64(elem, seed_ + 1);
        uint8_t rho = count_leading_zeros(h2) + 1;
        if (rho > M_[j]) M_[j] = rho;
    }

    [[nodiscard]] double estimate() const override {
        double sum = 0.0;
        for (std::size_t i = 0; i < size; ++i)
            sum += std::pow(2.0, -static_cast<double>(M_[i]));

        double m = static_cast<double>(size);
        double raw = alpha(size) * m * m / sum;

        // Small-range correction (linear counting)
        if (raw <= 2.5 * m) {
            std::size_t zeros = 0;
            for (auto v : M_) if (v == 0) ++zeros;
            if (zeros > 0) return m * std::log(m / static_cast<double>(zeros));
        }
        return raw;
    }

    void merge(const HyperLogLog& other) {
        if (other.size != size)
            throw std::invalid_argument("Cannot merge sketches of different sizes.");
        for (std::size_t i = 0; i < size; ++i)
            M_[i] = std::max(M_[i], other.M_[i]);
    }

    [[nodiscard]] std::vector<uint8_t> get_registers() const { return M_; }

    [[nodiscard]] size_t memory_usage(uint64_t flags) const override {
        uint64_t f = resolve_flags(flags);
        size_t s = 0;
        if (f & MemoryFlag::REGISTERS) s += M_.capacity();
        if (f & MemoryFlag::ALL_WRITE_NO_REGISTERS) s += sizeof(size);
        if (f & MemoryFlag::SEEDS) s += seeds_.bytes();
        return s;
    }

private:
    std::vector<uint8_t> M_;
    std::uint32_t seed_;   // single hash seed

    static uint8_t count_leading_zeros(std::uint64_t x) {
        if (x == 0) return 64;
        return static_cast<uint8_t>(__builtin_clzll(x));
    }

    static double alpha(std::size_t m) {
        if (m == 16) return 0.673;
        if (m == 32) return 0.697;
        if (m == 64) return 0.709;
        return 0.7213 / (1.0 + 1.079 / static_cast<double>(m));
    }
};
