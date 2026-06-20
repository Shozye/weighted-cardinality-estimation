#pragma once
#include <algorithm>
#include <cmath>
#include <vector>
#include <string>
#include <cstdint>
#include "fisher_yates.hpp"
#include "hash_util.hpp"
#include "quantize_custom_float.hpp"
#include "rng_engine_type.hpp"
#include "sketch.hpp"

class FastExpSketchCustomFloat : public Sketch, public MergeableMixin, public JaccardMixin {
public:
    FastExpSketchCustomFloat(
        std::size_t sketch_size,
        std::uint64_t master_seed,
        int exp_bits, int mant_bits,
        RngEngine engine = kDefaultRngEngine
    );
    FastExpSketchCustomFloat(
        std::size_t sketch_size,
        std::uint64_t master_seed,
        int exp_bits, int mant_bits,
        const std::vector<double>& registers,
        RngEngine engine = kDefaultRngEngine
    );

    void add(const std::string& elem, double weight = 1.0) override;
    [[nodiscard]] double estimate() const override;
    [[nodiscard]] double jaccard_struct(const FastExpSketchCustomFloat& other) const;
    void merge(const FastExpSketchCustomFloat& other);

    const std::vector<double>& get_registers() const;
    [[nodiscard]] FastExpSketchCustomFloat clone_with(int exp_bits, int mant_bits) const;

    [[nodiscard]] int get_exp_bits() const { return exp_bits_; }
    [[nodiscard]] int get_mant_bits() const { return mant_bits_; }

    [[nodiscard]] size_t memory_usage(uint64_t flags) const override;

private:
    static constexpr QuantizationMode kMode = QuantizationMode::ALL_NORMAL;
    int exp_bits_;
    int mant_bits_;
    std::vector<double> M_;
    FisherYates fisher_yates;
    double max_;
};
