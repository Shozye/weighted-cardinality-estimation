#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include "seeds.hpp"
#include "fisher_yates.hpp"

class FastExpSketch {
public:
    FastExpSketch(std::size_t size, const std::vector<std::uint32_t>& seeds);
    void add(const std::string& elem, double weight = 1.0);
    void add_many(const std::vector<std::string>& elems, const std::vector<double>& weights);
    [[nodiscard]] double estimate() const;
    [[nodiscard]] double jaccard_struct(const FastExpSketch& other) const;

    FastExpSketch(std::size_t sketch_size, const std::vector<std::uint32_t>& seeds, const std::vector<double>& registers);

    std::size_t get_sketch_size() const;
    std::vector<std::uint32_t> get_seeds() const;
    const std::vector<double>& get_registers() const;

    [[nodiscard]] size_t memory_usage_total() const;
    [[nodiscard]] size_t memory_usage_write() const;
    [[nodiscard]] size_t memory_usage_estimate() const;
private:
    std::size_t size; 
    Seeds seeds_;
    std::vector<double> M_;
    FisherYates fisher_yates;
    double max;    
};
