#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include "compact_vector.hpp"
#include "sketch.hpp"

class kQSketchRoundedDyn : public Sketch, public NewtonMixin {
public:
    kQSketchRoundedDyn(
        std::size_t sketch_size,
        std::uint64_t master_seed,
        std::uint8_t amount_bits,
        float logarithm_base,
        std::uint32_t g_seed = 42
    );
    kQSketchRoundedDyn(
        std::size_t sketch_size,
        std::uint8_t amount_bits,
        float logarithm_base,
        std::uint32_t g_seed,
        std::uint64_t master_seed,
        const std::vector<int>& registers,
        const std::vector<std::uint32_t>& t_histogram,
        double cardinality
    );

    void add(const std::string& elem, double weight = 1.0);
    [[nodiscard]] double estimate() const;
    [[nodiscard]] double estimate_direct() const;
    [[nodiscard]] double estimate_newton_cold() const;
    [[nodiscard]] double estimate_newton_warm() const;
    [[nodiscard]] int estimate_newton_cold_iterations() const;
    [[nodiscard]] int estimate_newton_warm_iterations() const;

    std::uint8_t get_amount_bits() const;
    float get_logarithm_base() const;
    std::uint32_t get_g_seed() const;
    std::vector<int> get_registers() const;
    std::vector<std::uint32_t> get_t_histogram() const;
    double get_cardinality() const;

    [[nodiscard]] size_t memory_usage(uint64_t flags) const override;

private:
    std::uint8_t amount_bits_;
    float logarithm_base_;
    std::int32_t r_min;
    std::int32_t r_max;
    std::uint32_t g_seed_;

    double cardinality_;
    double q_r_;
    compact::vector<int> R_;
    compact::vector<std::uint32_t> T_;
};
