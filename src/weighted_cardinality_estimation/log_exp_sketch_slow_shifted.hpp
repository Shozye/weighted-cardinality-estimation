#pragma once
#include "compact_vector.hpp"
#include "sketch.hpp"
#include <cstdint>
#include <string>
#include <vector>

class LogExpSketchSlowShifted : public Sketch, public MergeableMixin, public JaccardMixin {
public:
    LogExpSketchSlowShifted(
        std::size_t sketch_size,
        std::uint64_t master_seed,
        std::uint8_t amount_bits,
        double v_max
    );
    LogExpSketchSlowShifted(
        std::size_t sketch_size,
        std::uint64_t master_seed,
        std::uint8_t amount_bits,
        double v_max,
        const std::vector<int>& registers,
        int offset
    );

    void add(const std::string& elem, double weight = 1.0) override;
    [[nodiscard]] double estimate() const override;
    [[nodiscard]] double jaccard_struct(const LogExpSketchSlowShifted& other) const;
    void merge(const LogExpSketchSlowShifted& other);

    [[nodiscard]] std::vector<int> get_registers() const;
    [[nodiscard]] std::uint8_t get_amount_bits() const;
    [[nodiscard]] double get_v_max() const;
    [[nodiscard]] int get_offset() const;

    [[nodiscard]] size_t memory_usage(uint64_t flags) const override;

private:
    void shift_down();

    std::uint8_t amount_bits_;
    double v_max_;
    std::int32_t capacity_;   // (1 << amount_bits) - 1
    std::int32_t offset_;     // absolute index = M_[i] + offset_ (may be negative)
    std::int32_t num_maxed_;  // count of registers at capacity_ (sentinel "infinity")
    double log_r_;            // log ratio between adjacent grid points = ln(v_max)/capacity_
    compact::vector<unsigned> M_;

    [[nodiscard]] double reconstruct(int abs_index) const;
    [[nodiscard]] int quantize(double value) const;
};
