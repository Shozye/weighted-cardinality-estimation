#pragma once
#include "compact_vector.hpp"
#include "sketch.hpp"
#include <sys/types.h>
#include <vector>
#include <string>
#include <cstdint>

class BaseShiftedLogExpSketch : public Sketch {
public:
    BaseShiftedLogExpSketch(
        std::size_t sketch_size, 
        const std::vector<std::uint32_t>& seeds, 
        std::uint8_t amount_bits,
        float logarithm_base
    );
    BaseShiftedLogExpSketch(
        std::size_t sketch_size, 
        const std::vector<std::uint32_t>& seeds, 
        std::uint8_t amount_bits, 
        float logarithm_base,
        const std::vector<uint32_t>& registers,
        std::int32_t offset
    );
    void add(const std::string& elem, double weight = 1.0);
    [[nodiscard]] double estimate() const ;

    std::uint8_t get_amount_bits() const;
    std::vector<uint32_t> get_registers() const;
    float get_logarithm_base() const;
    std::int32_t get_offset() const;

    [[nodiscard]] size_t memory_usage_total() const;
    [[nodiscard]] size_t memory_usage_write() const;
    [[nodiscard]] size_t memory_usage_estimate() const;
private:
    double initialValue() const;
    double ffunc_divided_by_dffunc(double w) const;
    double Newton(double c0) const;

    std::uint8_t amount_bits_;
    float logarithm_base;

    std::uint32_t r_max; // maximum possible value of register (0 to 2**amount_bits-1)
    std::int32_t offset; // this is a relative value for all registers. real_value = value + offset

    compact::vector<uint32_t> M_; // sketch structure with elements between < r_min ... r_max >
};
