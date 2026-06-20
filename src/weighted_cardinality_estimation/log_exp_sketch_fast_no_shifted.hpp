#pragma once
#include "compact_vector.hpp"
#include <vector>
#include <string>
#include <cstdint>
#include "fisher_yates.hpp"
#include "hash_util.hpp"
#include "rng_engine_type.hpp"
#include "sketch.hpp"

class LogExpSketchFastNoShifted : public Sketch, public MergeableMixin, public JaccardMixin {
public:
    LogExpSketchFastNoShifted(
        std::size_t sketch_size,
        std::uint64_t master_seed,
        std::uint8_t amount_bits,
        double v_max,
        RngEngine engine = kDefaultRngEngine
    );
    LogExpSketchFastNoShifted(
        std::size_t sketch_size,
        std::uint64_t master_seed,
        std::uint8_t amount_bits,
        double v_max,
        const std::vector<int>& registers,
        RngEngine engine = kDefaultRngEngine
    );

    void add(const std::string& elem, double weight = 1.0) override;
    [[nodiscard]] double estimate() const override;
    [[nodiscard]] double jaccard_struct(const LogExpSketchFastNoShifted& other) const;
    void merge(const LogExpSketchFastNoShifted& other);

    [[nodiscard]] std::vector<int> get_registers() const;
    [[nodiscard]] std::uint8_t get_amount_bits() const;
    [[nodiscard]] double get_v_max() const;

    [[nodiscard]] size_t memory_usage(uint64_t flags) const override;

private:
    std::uint8_t amount_bits_;
    double v_max_;
    int r_max_;
    double min_value_;
    double log_r_;
    compact::vector<unsigned> M_;
    FisherYates fisher_yates;
    int max_register_;

    [[nodiscard]] double reconstruct(int index) const;
    [[nodiscard]] int quantize(double value) const;
    void update_max_register();
};
