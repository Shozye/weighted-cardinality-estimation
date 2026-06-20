#pragma once
#include "compact_vector.hpp"
#include <vector>
#include <string>
#include <cstdint>
#include <utility>
#include "fisher_yates.hpp"
#include "rng_engine_type.hpp"
#include "sketch.hpp"

class kQSketch : public Sketch, public MergeableMixin, public NewtonMixin {
public:
    static constexpr double newton_max_error = 1e-6;
    static constexpr int newton_max_iterations = 100;
    kQSketch(
        std::size_t sketch_size,
        std::uint64_t master_seed,
        std::uint8_t amount_bits,
        float logarithm_base,
        RngEngine engine = kDefaultRngEngine
    );
    kQSketch(
        std::size_t sketch_size,
        std::uint64_t master_seed,
        std::uint8_t amount_bits,
        float logarithm_base,
        const std::vector<int>& registers,
        RngEngine engine = kDefaultRngEngine
    );
    void add(const std::string& elem, double weight = 1.0);
    [[nodiscard]] double estimate() const;
    [[nodiscard]] double estimate_direct() const;
    [[nodiscard]] double estimate_newton_cold() const;
    [[nodiscard]] double estimate_newton_warm() const;
    [[nodiscard]] int estimate_newton_cold_iterations() const;
    [[nodiscard]] int estimate_newton_warm_iterations() const;

    std::uint8_t get_amount_bits() const;
    std::vector<int> get_registers() const;
    float get_logarithm_base() const;
    void merge(const kQSketch& other);

    [[nodiscard]] size_t memory_usage(uint64_t flags) const override;
private:
    double initialValue() const;
    double ffunc_divided_by_dffunc(double w) const;
    double Newton(double c0) const;
    std::pair<double, int> Newton_with_iterations(double c0) const;

    void update_treshold();

    FisherYates fisher_yates;
    std::uint8_t amount_bits_;
    float logarithm_base;
    std::int32_t r_max; // maximum possible value in sketch due to amount of bits per register
    std::int32_t r_min; // minimum possible value in sketch due to amount of bits per register

    compact::vector<int> M_; // sketch structure with elements between < r_min ... r_max >
    int min_sketch_value; 
    double min_value_to_change_sketch; // that's 2**{-min_sketch_value}
};
