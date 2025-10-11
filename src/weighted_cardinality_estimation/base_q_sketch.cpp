#include "base_q_sketch.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "hash_util.hpp"
#include "compact_vector.hpp"

#define NEWTON_MAX_ITERATIONS 5

BaseQSketch::BaseQSketch(std::size_t sketch_size, const std::vector<std::uint32_t>& seeds, uint8_t amount_bits)
    : size(sketch_size), 
      seeds_(seeds),
      amount_bits_(amount_bits),
      r_max((1 << (amount_bits - 1)) - 1),
      r_min(-(1 << (amount_bits - 1)) + 1),
      M_(amount_bits, sketch_size)
{
    if (seeds.size() != size) { throw std::invalid_argument("Seeds vector must have length m"); }
    for (std::size_t i = 0; i < size; ++i) {
        M_[i] = r_min;
    }
}

void BaseQSketch::add(const std::string& elem, double weight)
{ 
    auto inv_weight = 1.0 / weight;
    for (std::size_t i = 0; i < size; ++i) {
        std::uint64_t h = murmur64(elem, seeds_[i], hash_answer);
        double u = to_unit_interval(h);   
        double g = -std::log(u) * inv_weight;
        int q = static_cast<int>(std::floor(-std::log2(g)));
        q = std::min(q, r_max);
        if (q > M_[i]){
            M_[i] = q;
        }
    }
} 

size_t BaseQSketch::memory_usage_total() const {
    size_t total_size = 0;
    total_size += sizeof(size);
    total_size += seeds_.bytes();
    total_size += M_.bytes();
    total_size += sizeof(r_min);
    total_size += sizeof(r_max);
    total_size += sizeof(amount_bits_);
    return total_size;
}

size_t BaseQSketch::memory_usage_write() const {
    return M_.bytes();
}

size_t BaseQSketch::memory_usage_estimate() const {
    return M_.bytes();
}

void BaseQSketch::add_many(const std::vector<std::string>& elems,
                                  const std::vector<double>& weights) {
    if (elems.size() != weights.size()){
        throw std::invalid_argument("add_many: elems and weights size mismatch");
    }
    for (std::size_t i = 0; i < elems.size(); ++i) {
        this->add(elems[i], weights[i]);
    }
}


BaseQSketch::BaseQSketch(
    std::size_t sketch_size, 
    const std::vector<std::uint32_t>& seeds, 
    std::uint8_t amount_bits, 
    const std::vector<int>& registers
)
    : size(sketch_size), 
    seeds_(seeds), 
    amount_bits_(amount_bits),
    r_max((1 << (amount_bits - 1)) - 1),
    r_min(-(1 << (amount_bits - 1)) + 1),
    M_(amount_bits, sketch_size)
{
    for (std::size_t i = 0; i < size; ++i) {
        M_[i] = registers[i];
    }
}

std::size_t BaseQSketch::get_sketch_size() const {
    return size;
}

std::vector<std::uint32_t> BaseQSketch::get_seeds() const {
    return seeds_.toVector();
}

std::vector<int> BaseQSketch::get_registers() const {
    return std::vector<int>(M_.begin(), M_.end());
}

std::uint8_t BaseQSketch::get_amount_bits() const {
    return amount_bits_;
}


double BaseQSketch::initialValue() const {
    double c0 = 0.0;
    double tmp_sum = 0.0;
 
    for(size_t i=0; i<this->size; i++) { 
        tmp_sum += std::ldexp(1.0, -M_[i]);
    }

    c0 = (double)(this->size-1) / tmp_sum;
    return c0;
}

double BaseQSketch::ffunc(double w) const {
    double res = 0;
    for (size_t i = 0; i < size; ++i) {
        double x = std::ldexp(1.0, -M_[i] - 1);
        double ex = std::exp(w * x);
        res += x * (2.0 - ex) / (ex - 1.0);
    }
    return res;
}

double BaseQSketch::dffunc(double w) const {
    double res = 0;
    for (size_t i = 0; i < size; ++i) {
        double x = std::ldexp(1.0, -M_[i] - 1);
        double ex = std::exp(w * x);
        res += -x * x * ex * pow(ex - 1, -2);
    }
    return res;
}

double BaseQSketch::Newton(double c0) const {
    double err = 1e-5;
    double c1 = c0 - (ffunc(c0) / dffunc(c0));
    int it = 0;
    while (std::abs(c1 - c0) > err) {
        c0 = c1;
        c1 = c0 - ffunc(c0) / dffunc(c0);
        it += 1;
        if (it > NEWTON_MAX_ITERATIONS){ break; }
    }
    return c1;
}

double BaseQSketch::estimate() const {
    return Newton(initialValue());
}
