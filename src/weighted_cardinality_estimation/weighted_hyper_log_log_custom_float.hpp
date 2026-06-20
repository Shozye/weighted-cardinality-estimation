#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>
#include "sketch.hpp"
#include "hash_util.hpp"
#include "quantize_custom_float.hpp"

// Weighted HyperLogLog with CustomFloat register quantization.
// Combines WeightedHyperLogLog's stochastic averaging (2 seeds, h1 % m)
// with CustomFloat quantize-on-write register storage (ALL_NORMAL mode, positive-only).
class WeightedHyperLogLogCustomFloat : public Sketch, public MergeableMixin {
public:
    WeightedHyperLogLogCustomFloat(
        std::size_t sketch_size, std::uint64_t master_seed,
        int exp_bits, int mant_bits
    ) : Sketch(sketch_size, master_seed),
        exp_bits_(exp_bits), mant_bits_(mant_bits),
        M_(sketch_size, custom_float_max(exp_bits, mant_bits, kMode)),
        seed1_(seeds_[0]), seed2_(seeds_[1]) {}

    WeightedHyperLogLogCustomFloat(
        std::size_t sketch_size, std::uint64_t master_seed,
        int exp_bits, int mant_bits,
        const std::vector<double>& registers
    ) : Sketch(sketch_size, master_seed),
        exp_bits_(exp_bits), mant_bits_(mant_bits),
        M_(registers),
        seed1_(seeds_[0]), seed2_(seeds_[1]) {}

    void add(const std::string& elem, double weight = 1.0) override {
        validate_weight(weight);
        std::uint64_t h1 = murmur64(elem, seed1_);
        std::size_t k = h1 % size;
        std::uint64_t h2 = murmur64(elem, seed2_);
        double u = to_unit_interval(h2);
        double g = -std::log(u) / weight;
        double quantized = quantize_custom_float(g, 0, exp_bits_, mant_bits_, kMode);
        if (quantized < M_[k]) M_[k] = quantized;
    }

    [[nodiscard]] double estimate() const override {
        double max_val = custom_float_max(exp_bits_, mant_bits_, kMode);
        std::size_t K = 0;
        double sum = 0.0;
        for (std::size_t i = 0; i < size; ++i) {
            if (M_[i] < max_val) {
                sum += M_[i];
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

    void merge(const WeightedHyperLogLogCustomFloat& other) {
        if (other.size != size)
            throw std::invalid_argument("Cannot merge sketches of different sizes.");
        if (other.exp_bits_ != exp_bits_ || other.mant_bits_ != mant_bits_)
            throw std::invalid_argument("Cannot merge sketches with different float formats.");
        for (std::size_t i = 0; i < size; ++i)
            M_[i] = std::min(M_[i], other.M_[i]);
    }

    [[nodiscard]] const std::vector<double>& get_registers() const { return M_; }

    [[nodiscard]] int get_exp_bits() const { return exp_bits_; }
    [[nodiscard]] int get_mant_bits() const { return mant_bits_; }

    [[nodiscard]] size_t memory_usage(uint64_t flags) const override {
        uint64_t f = resolve_flags(flags);
        size_t s = 0;
        if (f & MemoryFlag::REGISTERS) s += M_.capacity() * sizeof(double);
        if (f & MemoryFlag::ALL_WRITE_NO_REGISTERS) s += sizeof(size);
        if (f & MemoryFlag::SEEDS) s += seeds_.bytes();
        return s;
    }

private:
    static constexpr QuantizationMode kMode = QuantizationMode::ALL_NORMAL;
    int exp_bits_;
    int mant_bits_;
    std::vector<double> M_;
    std::uint32_t seed1_;
    std::uint32_t seed2_;
};
