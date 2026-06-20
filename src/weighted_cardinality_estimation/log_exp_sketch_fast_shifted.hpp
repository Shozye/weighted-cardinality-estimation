#pragma once
#include "compact_vector.hpp"
#include "fisher_yates.hpp"
#include "rng_engine_type.hpp"
#include "sketch.hpp"
#include <cstdint>
#include <string>
#include <vector>

class LogExpSketchFastShifted : public Sketch, public MergeableMixin, public JaccardMixin {
public:
    LogExpSketchFastShifted(
        std::size_t sketch_size,
        std::uint64_t master_seed,
        std::uint8_t amount_bits,
        double v_max,
        RngEngine engine = kDefaultRngEngine
    );
    LogExpSketchFastShifted(
        std::size_t sketch_size,
        std::uint64_t master_seed,
        std::uint8_t amount_bits,
        double v_max,
        const std::vector<int>& registers,
        int offset,
        RngEngine engine = kDefaultRngEngine
    );

    void add(const std::string& elem, double weight = 1.0) override;
    [[nodiscard]] double estimate() const override;
    [[nodiscard]] double jaccard_struct(const LogExpSketchFastShifted& other) const;
    void merge(const LogExpSketchFastShifted& other);

    [[nodiscard]] std::vector<int> get_registers() const;
    [[nodiscard]] std::uint8_t get_amount_bits() const;
    [[nodiscard]] double get_v_max() const;
    [[nodiscard]] int get_offset() const;

    [[nodiscard]] size_t memory_usage(uint64_t flags) const override;

private:
    void shift_down();

    std::uint8_t amount_bits_;
    double v_max_;
    std::int32_t capacity_;
    std::int32_t offset_;
    std::int32_t num_maxed_;
    double min_value_;
    double log_r_;
    compact::vector<unsigned> M_;
    FisherYates fisher_yates;

    [[nodiscard]] double reconstruct(int abs_index) const;
    [[nodiscard]] int quantize(double value) const;
};
