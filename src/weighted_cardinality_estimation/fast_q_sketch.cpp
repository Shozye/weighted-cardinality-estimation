#include "fast_q_sketch.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "fast_exp_sketch.hpp"
#include "hash_util.hpp"
#include<cstring>

#define NEWTON_MAX_ITERATIONS 5

FastQSketch::FastQSketch(std::size_t sketch_size, const std::vector<std::uint32_t>& seeds, uint8_t amount_bits)
    : size(sketch_size), 
      seeds_(seeds),
      fisher_yates(size),
      amount_bits_(amount_bits),
      r_max((1 << (amount_bits - 1)) - 1),
      r_min(-(1 << (amount_bits - 1)) + 1),
      M_(amount_bits, sketch_size)
{
    if (sketch_size == 0) { throw std::invalid_argument("Sketch size 'm' must be positive."); }
    if (amount_bits == 0) { throw std::invalid_argument("Amount of bits 'b' must be positive."); }
    if ((!seeds.empty() && seeds.size() != size)) { 
        throw std::invalid_argument("Seeds must have length m or 0"); 
    }
    for (std::size_t i = 0; i < size; ++i) {
        M_[i] = r_min;
    }
    update_treshold();
}

FastQSketch::FastQSketch(std::size_t sketch_size, const std::vector<std::uint32_t>& seeds, std::uint8_t amount_bits, const std::vector<int>& registers)
    : size(sketch_size),
      seeds_(seeds),
      fisher_yates(size),
      amount_bits_(amount_bits),
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
    update_treshold();
}


size_t FastQSketch::memory_usage_total() const {
    size_t total_size = 0;
    total_size += sizeof(size);
    total_size += sizeof(amount_bits_);
    total_size += sizeof(r_max);
    total_size += sizeof(r_min);
    total_size += sizeof(min_sketch_value);
    total_size += sizeof(min_value_to_change_sketch);
    total_size += seeds_.bytes();
    total_size += M_.bytes();
    total_size += fisher_yates.bytes_total();
    return total_size;
}

size_t FastQSketch::memory_usage_write() const {
    size_t write_size = 0;
    write_size += sizeof(min_sketch_value);
    write_size += sizeof(min_value_to_change_sketch);
    write_size += M_.bytes();
    write_size += fisher_yates.bytes_write();
    return write_size;
}

size_t FastQSketch::memory_usage_estimate() const {
    size_t estimate_size = M_.bytes();
    return estimate_size;
}

std::size_t FastQSketch::get_sketch_size() const { return size; }
std::vector<std::uint32_t> FastQSketch::get_seeds() const { return seeds_.toVector(); }
std::uint8_t FastQSketch::get_amount_bits() const { return amount_bits_; }
std::vector<int> FastQSketch::get_registers() const {
    return std::vector<int>(M_.begin(), M_.end());
}

void FastQSketch::update_treshold(){
    this->min_sketch_value = *std::min_element(this->M_.begin(), this->M_.end());
    this->min_value_to_change_sketch = std::pow(2, -this->min_sketch_value);
}

void FastQSketch::add(const std::string& elem, double weight){ 
    double S = 0;
    bool touched_min = false; 

    fisher_yates.initialize(murmur64(elem, 1, hash_answer));
    for (size_t k = 0; k < this->size; ++k){
        std::uint64_t hashed = murmur64(elem, seeds_[k], hash_answer); 
        double unit_interval_hash = to_unit_interval(hashed); 
        double exponential_variable = -std::log(unit_interval_hash) / weight; 
        S += exponential_variable/(double)(this->size-k); 

        if ( S >= this->min_value_to_change_sketch ) { break; } 

        auto j = fisher_yates.get_fisher_yates_element(k);

        int q = static_cast<int>(std::floor(-std::log2(S)));
        q = std::min(q, r_max);

        if (q > this->M_[j]){
            if (this->M_[j] == min_sketch_value){
                touched_min = true;
            }
            this->M_[j] = q;
        }
    }

    if(touched_min){
        this->update_treshold();
    }
} 

void FastQSketch::add_many(const std::vector<std::string>& elems,
                                  const std::vector<double>& weights) {
    if (elems.size() != weights.size()){
        throw std::invalid_argument("add_many: elems and weights size mismatch");
    }
    for (std::size_t i = 0; i < elems.size(); ++i) {
        this->add(elems[i], weights[i]);
    }
}


double FastQSketch::initialValue(){
    double c0 = 0.0;
    double tmp_sum = 0.0;
 
    for(size_t i=0; i<this->size; i++) { 
        tmp_sum += std::ldexp(1.0, -M_[i]);
    }

    c0 = (double)(this->size-1) / tmp_sum;
    return c0;
}

double FastQSketch::ffunc(double w) {
    double res = 0;
    for (size_t i = 0; i < size; ++i) {
        double x = std::ldexp(1.0, -M_[i] - 1);
        double ex = std::exp(w * x);
        res += x * (2.0 - ex) / (ex - 1.0);
    }
    return res;
}

double FastQSketch::dffunc(double w) {
    double res = 0;
    for (size_t i = 0; i < size; ++i) {
        double x = std::ldexp(1.0, -M_[i] - 1);
        double ex = std::exp(w * x);
        res += -x * x * ex * pow(ex - 1, -2);
    }
    return res;
}

double FastQSketch::Newton(double c0) {
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

double FastQSketch::estimate() {
    return Newton(initialValue());
}
