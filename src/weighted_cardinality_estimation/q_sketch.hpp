#pragma once
#include "compact_vector.hpp"
#include <vector>
#include <string>
#include <cstdint>
#include "fisher_yates.hpp"
#include "rng_engine_type.hpp"
#include "sketch.hpp"

class QSketch : public Sketch, public MergeableMixin {
// Paper: https://arxiv.org/abs/2406.19143v1
public:
    static constexpr double newton_max_error = 1e-6;
    static constexpr int newton_max_iterations = 5;
    QSketch(
        std::size_t sketch_size,
        std::uint64_t master_seed,
        std::uint8_t amount_bits,
        RngEngine engine = kDefaultRngEngine
    );
    QSketch(
        std::size_t sketch_size,
        std::uint64_t master_seed,
        std::uint8_t amount_bits,
        const std::vector<int>& registers,
        RngEngine engine = kDefaultRngEngine
    );

    void add(const std::string& elem, double weight = 1.0);
    [[nodiscard]] double estimate() const;
    std::uint8_t get_amount_bits() const;
    std::vector<int> get_registers() const;
    void merge(const QSketch& other);

    [[nodiscard]] size_t memory_usage(uint64_t flags) const override;
private:
    double initialValue() const;
    double ffunc_divided_by_dffunc(double w) const;
    double Newton(double c0) const;

    FisherYates fisher_yates;
    std::uint8_t amount_bits_;
    std::int32_t r_max; // maximum possible value in sketch due to amount of bits per register
    std::int32_t r_min; // minimum possible value in sketch due to amount of bits per register

    compact::vector<int> M_; // sketch structure with elements between < r_min ... r_max >
    uint32_t j_star;
};
