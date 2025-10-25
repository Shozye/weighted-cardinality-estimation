#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "hash_util.hpp"
#include<cstring>
#include "utils.hpp"
#include"base_log_exp_sketch_jacc.hpp"

BaseLogExpSketchJacc::BaseLogExpSketchJacc(
    std::size_t sketch_size, 
    const std::vector<std::uint32_t>& seeds, 
    uint8_t amount_bits,
    float logarithm_base,
    std::uint8_t amount_bits_jaccard
)
    : Sketch(sketch_size, seeds),
      amount_bits_(amount_bits),
      logarithm_base(logarithm_base),
      r_max((1 << (amount_bits - 1)) - 1),
      r_min(-(1 << (amount_bits - 1)) + 1),
      amount_bits_jaccard(amount_bits_jaccard),
      M_(amount_bits, sketch_size),
      H_(amount_bits_jaccard, sketch_size)
{
    if (amount_bits == 0) { throw std::invalid_argument("Amount of bits 'b' must be positive."); }
    std::fill(M_.begin(), M_.end(), r_min);
}

BaseLogExpSketchJacc::BaseLogExpSketchJacc(
    std::size_t sketch_size, 
    const std::vector<std::uint32_t>& seeds, 
    std::uint8_t amount_bits, 
    float logarithm_base,
    std::uint8_t amount_bits_jaccard,
    const std::vector<int>& registers,
    const std::vector<std::uint32_t>& h_registers
)
    : Sketch(sketch_size, seeds),
      amount_bits_(amount_bits),
      logarithm_base(logarithm_base),
      r_max((1 << (amount_bits - 1)) - 1),
      r_min(-(1 << (amount_bits - 1)) + 1),
      amount_bits_jaccard(amount_bits_jaccard),
      M_(amount_bits, sketch_size),
      H_(amount_bits_jaccard, sketch_size, h_registers)
{
    if (amount_bits == 0) { throw std::invalid_argument("Amount of bits 'b' must be positive."); }
    if (registers.size() != sketch_size) { throw std::invalid_argument("Invalid state: registers vector size mismatch"); }
    for (std::size_t i = 0; i < size; ++i) {
        M_[i] = registers[i];
    }
}


size_t BaseLogExpSketchJacc::memory_usage_total() const {
    size_t size = 0;
    size += sizeof(size); // 8
    size += seeds_.bytes(); // m * 4
    size += M_.bytes(); // mb/8
    size += sizeof(amount_bits_); // 1
    size += sizeof(r_max); // 4
    size += sizeof(r_min); // 4
    size += sizeof(logarithm_base); // 4
    size += H_.memory_usage_total(); // mb/8 + 1
    return size; 
}

size_t BaseLogExpSketchJacc::memory_usage_write() const {
    size_t size = 0;
    size += M_.bytes(); // mb/8
    size += H_.memory_usage_write(); // mb/8
    return size; 
}

size_t BaseLogExpSketchJacc::memory_usage_estimate() const {
    size_t size = M_.bytes(); // mb/8
    size += sizeof(logarithm_base); // 4
    return size; // mb/8 + 4
}

double BaseLogExpSketchJacc::jaccard_struct(const BaseLogExpSketchJacc& other) const {
    return this->H_.compute_jaccard(other.H_);
}

std::uint8_t BaseLogExpSketchJacc::get_amount_bits() const { return amount_bits_; }
float BaseLogExpSketchJacc::get_logarithm_base() const { return logarithm_base; }
std::vector<int> BaseLogExpSketchJacc::get_registers() const {
    return std::vector<int>(M_.begin(), M_.end());
}
std::vector<std::uint32_t> BaseLogExpSketchJacc::get_h_registers() const {
    return H_.get_h_registers();
}
std::uint8_t BaseLogExpSketchJacc::get_amount_bits_jaccard() const {
    return amount_bits_jaccard;
}

void BaseLogExpSketchJacc::add(const std::string& elem, double weight){ 
    for (std::size_t i = 0; i < size; ++i) {
        std::uint64_t h = murmur64(elem, seeds_[i]);
        double u = to_unit_interval(h);   
        double g = -std::log(u) / weight;
        int q = static_cast<int>(std::floor(-std::log(g)/std::log(logarithm_base)));
        q = std::min(q, r_max);
        if (q > M_[i]){
            M_[i] = q;
            H_.set_elem(i, elem);
        }
    }
} 

double BaseLogExpSketchJacc::initialValue() const {
    double tmp_sum = 0.0;
    for(int r: M_) { 
        tmp_sum += std::pow(logarithm_base, -r);
    }
    return (double)(this->size-1) / tmp_sum;
}

double BaseLogExpSketchJacc::ffunc_divided_by_dffunc(double w) const {
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

double BaseLogExpSketchJacc::Newton(double c0) const {
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

double BaseLogExpSketchJacc::estimate() const {
    return Newton(initialValue());
}
