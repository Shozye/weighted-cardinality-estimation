#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include "fisher_yates.hpp"
#include "rng_engine_type.hpp"
#include "sketch.hpp"

class FastGMExpSketch : public Sketch, public MergeableMixin, public JaccardMixin {
// Paper: https://arxiv.org/abs/2302.05176 
public:
    FastGMExpSketch(std::size_t sketch_size, std::uint64_t master_seed, RngEngine engine = kDefaultRngEngine);
    FastGMExpSketch(std::size_t sketch_size, std::uint64_t master_seed, const std::vector<double>& registers, RngEngine engine = kDefaultRngEngine);
    
    void add(const std::string& elem, double weight = 1.0);
    [[nodiscard]] double estimate() const;
    [[nodiscard]] double jaccard_struct(const FastGMExpSketch& other) const;

    const std::vector<double>& get_registers() const;
    void merge(const FastGMExpSketch& other);

    [[nodiscard]] size_t memory_usage(uint64_t flags) const override;
private:
    std::vector<double> M_;
    FisherYates fisher_yates;

    uint32_t j_star;
    uint32_t k_star;
    bool flagFastPrune;
};
