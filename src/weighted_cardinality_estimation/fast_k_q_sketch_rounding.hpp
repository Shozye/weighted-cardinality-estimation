#pragma once
#include "compact_vector.hpp"
#include <vector>
#include <string>
#include <cstdint>
#include "fisher_yates.hpp"
#include "rng_engine_type.hpp"
#include "sketch.hpp"

class kQSketchRounding : public Sketch, public MergeableMixin, public NewtonMixin {
public:
    kQSketchRounding(
        std::size_t sketch_size,
        std::uint64_t master_seed,
        std::uint8_t amount_bits,
        float logarithm_base,
        RngEngine engine = kDefaultRngEngine
    );
    kQSketchRounding(
        std::size_t sketch_size,
        std::uint64_t master_seed,
        std::uint8_t amount_bits,
        float logarithm_base,
        const std::vector<int>& registers,
        RngEngine engine = kDefaultRngEngine
    );
    void add(const std::string& elem, double weight = 1.0);
    [[nodiscard]] double estimate() const;
    [[nodiscard]] double estimate_corrected() const;
    [[nodiscard]] double estimate_direct() const;
    [[nodiscard]] double estimate_newton_cold() const;
    [[nodiscard]] double estimate_newton_warm() const;
    [[nodiscard]] int estimate_newton_cold_iterations() const;
    [[nodiscard]] int estimate_newton_warm_iterations() const;

    std::uint8_t get_amount_bits() const;
    std::vector<int> get_registers() const;
    float get_logarithm_base() const;
    void merge(const kQSketchRounding& other);

    [[nodiscard]] size_t memory_usage(uint64_t flags) const override;
private:
    void update_treshold();

    FisherYates fisher_yates;
    std::uint8_t amount_bits_;
    float logarithm_base;
    std::int32_t r_max;
    std::int32_t r_min;

    compact::vector<int> M_;
    int min_sketch_value;
    double min_value_to_change_sketch;
};
