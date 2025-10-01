#include "fast_q_sketch.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "hash_util.hpp"
#include<cstring>
#include"utils.hpp"

#define NEWTON_MAX_ITERATIONS 5

FastQSketch::FastQSketch(std::size_t sketch_size, const std::vector<std::uint32_t>& seeds, uint8_t amount_bits)
    : size(sketch_size), 
      seeds_(seeds),
      amount_bits_(amount_bits),
      r_max(mypow(2, amount_bits-1) - 1),
      r_min(-mypow(2, amount_bits-1) + 1),
      M_(sketch_size, r_min), 
      permInit(sketch_size),
      permWork(sketch_size),
      rng_seed(0)
{
    if (seeds_.size() != sketch_size) {
        throw std::invalid_argument("Seeds vector must have length m");
    }
    for(size_t i = 0; i < sketch_size; i++){
        permInit[i] = i+1;
    }
    update_treshold();
}

FastQSketch::FastQSketch(std::size_t sketch_size, const std::vector<std::uint32_t>& seeds, std::uint8_t amount_bits, const std::vector<int>& registers)
    : size(sketch_size),
      seeds_(seeds),
      amount_bits_(amount_bits),
      r_max(mypow(2, amount_bits-1) - 1),
      r_min(-mypow(2, amount_bits-1) + 1),
      M_(registers),
      permInit(sketch_size),
      permWork(sketch_size),
      rng_seed(0) // add nadpisuje stan rng_seed
{
    if (seeds_.size() != sketch_size) { throw std::invalid_argument("Invalid state: seeds vector size mismatch"); }
    if (M_.size() != sketch_size) { throw std::invalid_argument("Invalid state: registers vector size mismatch"); }

    for(size_t i = 0; i < sketch_size; i++){
        permInit[i] = i+1;
    }
    update_treshold();
}


size_t FastQSketch::memory_usage_total() const {
    size_t total_size = 0;
    total_size += sizeof(size);
    total_size += sizeof(amount_bits_);
    total_size += sizeof(r_max);
    total_size += sizeof(r_min);
    total_size += sizeof(rng_seed);
    total_size += sizeof(min_sketch_value);
    total_size += sizeof(min_value_to_change_sketch);
    total_size += seeds_.capacity() * sizeof(uint32_t);
    total_size += M_.capacity() * sizeof(int);
    total_size += permInit.capacity() * sizeof(uint32_t);
    total_size += permWork.capacity() * sizeof(uint32_t);
    return total_size;
}

size_t FastQSketch::memory_usage_write() const {
    size_t write_size = 0;
    write_size += sizeof(rng_seed);
    write_size += sizeof(min_sketch_value);
    write_size += sizeof(min_value_to_change_sketch);
    write_size += M_.capacity() * sizeof(int);
    write_size += permWork.capacity() * sizeof(uint32_t);
    return write_size;
}

size_t FastQSketch::memory_usage_estimate() const {
    size_t estimate_size = M_.capacity() * sizeof(double);
    return estimate_size;
}

std::size_t FastQSketch::get_sketch_size() const { return size; }
const std::vector<std::uint32_t>& FastQSketch::get_seeds() const { return seeds_; }
std::uint8_t FastQSketch::get_amount_bits() const { return amount_bits_; }
const std::vector<int>& FastQSketch::get_registers() const { return M_; }

uint32_t FastQSketch::rand(uint32_t min, uint32_t max){
    this->rng_seed = this->rng_seed * 1103515245 + 12345;
    auto temp = (unsigned)(this->rng_seed/65536) % 32768;
    return (temp % (max-min)) + min;
}

void FastQSketch::update_treshold(){
    this->min_sketch_value = *std::min_element(this->M_.begin(), this->M_.end());
    this->min_value_to_change_sketch = std::pow(2, -this->min_sketch_value);
}

void FastQSketch::add(const std::string& elem, double weight){ 
    double S = 0;
    bool touched_min = false; 

    this->rng_seed = murmur64(elem, 1, hash_answer); 
    permWork = permInit; 
    auto inv_weight = 1.0 / weight;
    for (size_t k = 0; k < this->size; ++k){
        std::uint64_t hashed = murmur64(elem, seeds_[k], hash_answer); 
        double unit_interval_hash = to_unit_interval(hashed); 
        double exponential_variable = -std::log(unit_interval_hash) * inv_weight; 
        S += exponential_variable/(double)(this->size-k); 

        if ( S >= this->min_value_to_change_sketch ) { break; } 

        uint32_t r = rand(k, this->size);
        std::swap(permWork[k], permWork[r]);
        auto j = permWork[k] - 1;

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
    double e = 2.718282;
    for (size_t i = 0; i < size; ++i) {
        double x = std::ldexp(1.0, -M_[i] - 1);
        double ex = pow(e, w * x);
        res += x * (2.0 - ex) / (ex - 1.0);
    }
    return res;
}

double FastQSketch::dffunc(double w) {
    double res = 0;
    double e = 2.718282;
    for (size_t i = 0; i < size; ++i) {
        double x = std::ldexp(1.0, -M_[i] - 1);
        double ex = pow(e, w * x);
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
