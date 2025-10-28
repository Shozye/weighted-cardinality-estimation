#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include "sketch.hpp"

class ExpSketchFloat : public Sketch {
public:
    ExpSketchFloat(
        std::size_t sketch_size, 
        const std::vector<std::uint32_t>& seeds
    );
    ExpSketchFloat(
        std::size_t sketch_size, 
        const std::vector<std::uint32_t>& seeds, 
        const std::vector<float>& registers
    );

    void add(const std::string& elem, double weight = 1.0);
    [[nodiscard]] double estimate() const;
    [[nodiscard]] double jaccard_struct(const ExpSketchFloat& other) const;

    const std::vector<float>& get_registers() const;
    [[nodiscard]] size_t memory_usage_total() const;
    [[nodiscard]] size_t memory_usage_write() const;
    [[nodiscard]] size_t memory_usage_estimate() const;
private:
    std::vector<float> M_;
};
