#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include "compact_vector.hpp"
#include "seeds.hpp"

class BaseQSketch {
public:
    BaseQSketch(std::size_t sketch_size, const std::vector<std::uint32_t>& seeds, std::uint8_t amount_bits);
    void add(const std::string& elem, double weight = 1.0);
    void add_many(const std::vector<std::string>& elems, const std::vector<double>& weights);
    [[nodiscard]] double estimate() const;

    BaseQSketch(
        std::size_t sketch_size, 
        const std::vector<std::uint32_t>& seeds, 
        std::uint8_t amount_bits,
        const std::vector<int>& registers
    );

    std::size_t get_sketch_size() const;
    std::vector<std::uint32_t> get_seeds() const;
    std::vector<int> get_registers() const;
    std::uint8_t get_amount_bits() const;

    [[nodiscard]] size_t memory_usage_total() const;
    [[nodiscard]] size_t memory_usage_write() const;
    [[nodiscard]] size_t memory_usage_estimate() const;
private:
    std::size_t size;
    Seeds seeds_;
    std::uint8_t amount_bits_;
    std::int32_t r_max; // maximum possible value in sketch due to amount of bits per register
    std::int32_t r_min; // minimum possible value in sketch due to amount of bits per register
    compact::vector<int> M_; // sketch structure with elements between < r_min ... r_max >
    
    double initialValue() const;
    double ffunc(double w) const;
    double dffunc(double w) const;
    double Newton(double c0) const;
};
