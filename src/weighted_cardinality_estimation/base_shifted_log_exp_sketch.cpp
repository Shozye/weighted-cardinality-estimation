#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include "hash_util.hpp"
#include<cstring>
#include "utils.hpp"
#include"base_shifted_log_exp_sketch.hpp"

BaseShiftedLogExpSketch::BaseShiftedLogExpSketch(
    std::size_t sketch_size, 
    const std::vector<std::uint32_t>& seeds, 
    uint8_t amount_bits,
    float logarithm_base
)
    : Sketch(sketch_size, seeds),
      amount_bits_(amount_bits),
      logarithm_base(logarithm_base),
      r_max((1 << amount_bits) - 1),
      offset(-(1 << (amount_bits - 1)) + 1),
      M_(amount_bits, sketch_size)
{
    if (amount_bits == 0) { throw std::invalid_argument("Amount of bits 'b' must be positive."); }
    std::fill(M_.begin(), M_.end(), 0);
}

BaseShiftedLogExpSketch::BaseShiftedLogExpSketch(
    std::size_t sketch_size, 
    const std::vector<std::uint32_t>& seeds, 
    std::uint8_t amount_bits, 
    float logarithm_base,
    const std::vector<std::uint32_t>& registers,
    std::int32_t offset
)
    : Sketch(sketch_size, seeds),
      amount_bits_(amount_bits),
      logarithm_base(logarithm_base),
      r_max((1 << (amount_bits)) - 1),
      offset(offset),
      M_(amount_bits, sketch_size)
{
    if (amount_bits == 0) { throw std::invalid_argument("Amount of bits 'b' must be positive."); }
    if (registers.size() != sketch_size) { throw std::invalid_argument("Invalid state: registers vector size mismatch"); }
    for (std::size_t i = 0; i < size; ++i) {
        M_[i] = registers[i];
    }
}


size_t BaseShiftedLogExpSketch::memory_usage_total() const {
    size_t total_size = 0;
    total_size += sizeof(size);
    total_size += sizeof(amount_bits_);
    total_size += sizeof(r_max);
    total_size += sizeof(offset);
    total_size += sizeof(logarithm_base);
    total_size += seeds_.bytes();
    total_size += M_.bytes();
    return total_size;
}

size_t BaseShiftedLogExpSketch::memory_usage_write() const {
    size_t write_size = 0;
    write_size += M_.bytes();
    return write_size;
}

size_t BaseShiftedLogExpSketch::memory_usage_estimate() const {
    size_t estimate_size = M_.bytes();
    estimate_size += sizeof(logarithm_base);
    estimate_size += sizeof(offset);
    return estimate_size;
}

std::uint8_t BaseShiftedLogExpSketch::get_amount_bits() const { return amount_bits_; }
float BaseShiftedLogExpSketch::get_logarithm_base() const { return logarithm_base; }
std::int32_t BaseShiftedLogExpSketch::get_offset() const { return offset; }
std::vector<std::uint32_t> BaseShiftedLogExpSketch::get_registers() const {
    return std::vector<std::uint32_t>(M_.begin(), M_.end());
}

void BaseShiftedLogExpSketch::add(const std::string& elem, double weight){ 
    // print_vector(this->get_registers());
    for (std::size_t i = 0; i < size; ++i) {
        std::uint64_t h = murmur64(elem, seeds_[i]);
        double u = to_unit_interval(h);   
        double g = -std::log(u) / weight;
        int q = static_cast<int>(std::floor(-std::log(g)/std::log(logarithm_base)));
        if(q - offset < 0) {
            continue;
        }
        std::uint32_t possible = q - offset;
        if (possible == M_[i]) {
            continue; // jaccard usage
        } 
        // std::cout << "index i=" << i << " possible=" << possible << " rmax=" << r_max << " offset=" << offset << '\n';
        if (M_[i] < possible) { 
            if (possible > r_max) {
                std::uint32_t increase_offset = possible - r_max; // this is safe.
                for (std::size_t j = 0; j < M_.size(); j++) {
                    if (M_[j] >= increase_offset) {
                        M_[j] = M_[j] - increase_offset;
                    } else {
                        M_[j] = 0;
                    }
                }
                offset += (int)increase_offset; // theoretically it's not safe but who cares :)
                M_[i] = r_max;
            } else {
                // std::cout << "possible=" << possible << '\n';
                M_[i] = possible;
            }
        }
        // print_vector(this->get_registers());
    }
} 

double BaseShiftedLogExpSketch::initialValue() const {
    double tmp_sum = 0.0;
    for(std::uint32_t r: M_) { 
        tmp_sum += std::pow(logarithm_base, -(int)(r+offset));
    }
    return (double)(this->size-1) / tmp_sum;
}

double BaseShiftedLogExpSketch::ffunc_divided_by_dffunc(double w) const {
    double ffunc = 0;
    double dffunc = 0;
    for (std::uint32_t r: M_) {
        double x = std::pow(logarithm_base, -(int)(r+offset) - 1);;
        double ex = std::exp(w * x);
        ffunc += x * (2.0 - ex) / (ex - 1.0);
        dffunc += -x * x * ex * pow(ex - 1, -2);
    }
    return ffunc / dffunc;
}

double BaseShiftedLogExpSketch::Newton(double c0) const {
    double c1 = c0 - ffunc_divided_by_dffunc(c0);
    int it = 0;
    while (std::abs(c1 - c0) > NEWTON_MAX_ERROR) {
        c0 = c1;
        c1 = c0 - ffunc_divided_by_dffunc(c0);
        it += 1;
        if (it > NEWTON_MAX_ITERATIONS){ break; }
    }
    return c1;
}

double BaseShiftedLogExpSketch::estimate() const {
    return Newton(initialValue());
}
