#pragma once
#include "compact_vector.hpp"
#include <vector>
#include <string>
#include <cstdint>
#include "sketch.hpp"

class LogExpSketchSlowNoShifted : public Sketch, public MergeableMixin, public JaccardMixin {
public:
    LogExpSketchSlowNoShifted(
        std::size_t sketch_size,
        std::uint64_t master_seed,
        std::uint8_t amount_bits,
        double v_max
    );
    LogExpSketchSlowNoShifted(
        std::size_t sketch_size,
        std::uint64_t master_seed,
        std::uint8_t amount_bits,
        double v_max,
        const std::vector<int>& registers
    );

    void add(const std::string& elem, double weight = 1.0) override;
    [[nodiscard]] double estimate() const override;
    [[nodiscard]] double jaccard_struct(const LogExpSketchSlowNoShifted& other) const;
    void merge(const LogExpSketchSlowNoShifted& other);

    [[nodiscard]] std::vector<int> get_registers() const;
    [[nodiscard]] std::uint8_t get_amount_bits() const;
    [[nodiscard]] double get_v_max() const;

    [[nodiscard]] size_t memory_usage(uint64_t flags) const override;

private:
    std::uint8_t amount_bits_;
    double v_max_;
    int r_max_;        // N-1 = max index (initial value)
    double min_value_; // smallest representable value on the grid
    double log_r_;     // log ratio between adjacent grid points
    compact::vector<unsigned> M_;

    [[nodiscard]] double reconstruct(int index) const;
    [[nodiscard]] int quantize(double value) const;
};
