#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "hash_util.hpp"
#include<cstring>
#include "utils.hpp"
#include"base_log_exp_sketch.hpp"

BaseLogExpSketch::BaseLogExpSketch(
    std::size_t sketch_size, 
    const std::vector<std::uint32_t>& seeds, 
    uint8_t amount_bits,
    float logarithm_base
)
    : size(sketch_size), 
      seeds_(seeds),
      amount_bits_(amount_bits),
      logarithm_base(logarithm_base),
      r_max((1 << (amount_bits - 1)) - 1),
      r_min(-(1 << (amount_bits - 1)) + 1),
      M_(amount_bits, sketch_size)
{
    if (sketch_size == 0) { throw std::invalid_argument("Sketch size 'm' must be positive."); }
    if (amount_bits == 0) { throw std::invalid_argument("Amount of bits 'b' must be positive."); }
    if ((!seeds.empty() && seeds.size() != size)) { 
        throw std::invalid_argument("Seeds must have length m or 0"); 
    }
    std::fill(M_.begin(), M_.end(), r_min);
}

BaseLogExpSketch::BaseLogExpSketch(
    std::size_t sketch_size, 
    const std::vector<std::uint32_t>& seeds, 
    std::uint8_t amount_bits, 
    float logarithm_base,
    const std::vector<int>& registers
)
    : size(sketch_size),
      seeds_(seeds),
      amount_bits_(amount_bits),
      logarithm_base(logarithm_base),
      r_max((1 << (amount_bits - 1)) - 1),
      r_min(-(1 << (amount_bits - 1)) + 1),
      M_(amount_bits, sketch_size)
{
    if (sketch_size == 0) { throw std::invalid_argument("Sketch size 'm' must be positive."); }
    if (amount_bits == 0) { throw std::invalid_argument("Amount of bits 'b' must be positive."); }
    if ((!seeds.empty() && seeds.size() != size)) { 
        throw std::invalid_argument("Seeds must have length m or 0"); 
    }
    if (M_.size() != sketch_size) { throw std::invalid_argument("Invalid state: registers vector size mismatch"); }
    for (std::size_t i = 0; i < size; ++i) {
        M_[i] = registers[i];
    }
}


size_t BaseLogExpSketch::memory_usage_total() const {
    size_t total_size = 0;
    total_size += sizeof(size);
    total_size += sizeof(amount_bits_);
    total_size += sizeof(r_max);
    total_size += sizeof(r_min);
    total_size += sizeof(logarithm_base);
    total_size += seeds_.bytes();
    total_size += M_.bytes();
    return total_size;
}

size_t BaseLogExpSketch::memory_usage_write() const {
    size_t write_size = 0;
    write_size += M_.bytes();
    return write_size;
}

size_t BaseLogExpSketch::memory_usage_estimate() const {
    size_t estimate_size = M_.bytes();
    estimate_size += sizeof(logarithm_base);
    return estimate_size;
}

std::size_t BaseLogExpSketch::get_sketch_size() const { return size; }
std::vector<std::uint32_t> BaseLogExpSketch::get_seeds() const { return seeds_.toVector(); }
std::uint8_t BaseLogExpSketch::get_amount_bits() const { return amount_bits_; }
float BaseLogExpSketch::get_logarithm_base() const { return logarithm_base; }
std::vector<int> BaseLogExpSketch::get_registers() const {
    return std::vector<int>(M_.begin(), M_.end());
}

void BaseLogExpSketch::add(const std::string& elem, double weight){ 
    for (std::size_t i = 0; i < size; ++i) {
        std::uint64_t h = murmur64(elem, seeds_[i]);
        double u = to_unit_interval(h);   
        double g = -std::log(u) / weight;
        int q = static_cast<int>(std::floor(-std::log(g)/std::log(logarithm_base)));
        q = std::min(q, r_max);
        if (q > M_[i]){
            M_[i] = q;
        }
    }
} 

void BaseLogExpSketch::add_many(
    const std::vector<std::string>& elems,
    const std::vector<double>& weights
) {
    if (elems.size() != weights.size()){
        throw std::invalid_argument("add_many: elems and weights size mismatch");
    }
    for (std::size_t i = 0; i < elems.size(); ++i) {
        this->add(elems[i], weights[i]);
    }
}

double BaseLogExpSketch::initialValue() const {
    double tmp_sum = 0.0;
    for(int r: M_) { 
        tmp_sum += std::pow(logarithm_base, -r);
    }
    return (double)(this->size-1) / tmp_sum;
}

double BaseLogExpSketch::ffunc_divided_by_dffunc(double w) const {
    double ffunc = 0;
    double dffunc = 0;
    for (int r: M_) {
        double x = std::pow(logarithm_base, -r - 1);;
        double ex = std::exp(w * x);
        ffunc += x * (2.0 - ex) / (ex - 1.0);
        dffunc += -x * x * ex * pow(ex - 1, -2);
    }
    return ffunc / dffunc;
}

double BaseLogExpSketch::Newton(double c0) const {
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

double BaseLogExpSketch::estimate() {
    return Newton(initialValue());
}
