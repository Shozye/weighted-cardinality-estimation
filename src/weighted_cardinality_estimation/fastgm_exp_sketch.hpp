#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include "seeds.hpp"
#include "fisher_yates.hpp"

class FastGMExpSketch {
// Paper: https://arxiv.org/abs/2302.05176 
public:
    FastGMExpSketch(std::size_t sketch_size, const std::vector<std::uint32_t>& seeds);
    void add(const std::string& elem, double weight = 1.0);
    void add_many(const std::vector<std::string>& elems, const std::vector<double>& weights);
    [[nodiscard]] double estimate() const;
    [[nodiscard]] double jaccard_struct(const FastGMExpSketch& other) const;

    FastGMExpSketch(std::size_t sketch_size, const std::vector<std::uint32_t>& seeds, const std::vector<double>& registers);

    std::size_t get_sketch_size() const;
    std::vector<std::uint32_t> get_seeds() const;
    const std::vector<double>& get_registers() const;

    [[nodiscard]] size_t memory_usage_total() const;
    [[nodiscard]] size_t memory_usage_write() const;
    [[nodiscard]] size_t memory_usage_estimate() const;
private:
    std::size_t size;
    std::vector<double> M_;

    Seeds seeds_;
    FisherYates fisher_yates;

    uint32_t j_star;
    uint32_t k_star;
    bool flagFastPrune;
};
