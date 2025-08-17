#include "fast_q_sketch.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "hash_util.hpp"
#include<cstring>
#include<iostream>
#include"utils.hpp"

FastQSketch::FastQSketch(std::size_t m, const std::vector<std::uint32_t>& seeds, uint8_t amount_bits)
    : m_(m), 
      seeds_(seeds),
      r_max(mypow(2, amount_bits-1) - 1),
      r_min(-mypow(2, amount_bits-1) + 1),
      M_(m, r_min), // we are putting minimums because we are storing floor(- log_2 x ) 
      permInit(m),
      permWork(m),
      rng_seed(0)
{
    if (seeds_.size() != m_) {
        throw std::invalid_argument("Seeds vector must have length m");
    }
    for(size_t i = 0; i < m_; i++){
        permInit[i] = i+1;
    }
    update_treshold();
    // std::cout << "End of initialization: r_min=" << r_min << ", r_max=" << r_max << "\n";
    // print_vector(M_);
}

int FastQSketch::rand(int min, int max){
    this->rng_seed = this->rng_seed * 1103515245 + 12345;
    auto temp = (unsigned)(this->rng_seed/65536) % 32768;
    return (temp % (max-min)) + min;
}

void FastQSketch::update_treshold(){
    this->min_sketch_value = *std::min_element(this->M_.begin(), this->M_.end());
    this->min_value_to_change_sketch = std::pow(2, -this->min_sketch_value);
}

void FastQSketch::add(const std::string& x, double weight){ 
    uint64_t hash_answer[2];
    double S = 0;
    bool touched_min = false; 

    this->rng_seed = murmur64(x, 1, hash_answer); 
    permWork = permInit; 
    for (size_t k = 0; k < this->m_; ++k){
        std::uint64_t hashed = murmur64(x, seeds_[k], hash_answer); 
        double U = to_unit_interval(hashed); 
        double E = -std::log(U) / weight; 
        S += E/(double)(this->m_-k); 

        // std::cout << "I am during addition of element S=" << S << ", this->min_value_to_change_sketch=" << this->min_value_to_change_sketch << '\n'; 
        if ( S >= this->min_value_to_change_sketch ) { break; } 

        uint32_t r = rand(k, m_);
        auto swap = permWork[k];
        permWork[k] = permWork[r];
        permWork[r] = swap;
        auto j = permWork[k] - 1;

        int q = static_cast<int>(std::floor(-std::log2(S)));
        q = std::min(q, r_max);
        // std::cout << "I am during addition of element S=" << S << ", q=" << q << '\n'; 

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


double FastQSketch::initialValue(){
    double c0 = 0.0;
    double tmp_sum = 0.0;
 
    for(size_t i=0; i<this->m_; i++) { 
        tmp_sum += pow(2, -M_[i]); 
    }

    c0 = (double)(this->m_-1) / tmp_sum;
    return c0;
}

double FastQSketch::ffunc(double w) {
    double res = 0;
    double e = 2.718282;
    for (size_t i = 0; i < m_; ++i) {
        double x = pow(2.0, -M_[i] - 1);
        double ex = pow(e, w * x);
        res += x * (2.0 - ex) / (ex - 1.0);
    }
    return res;
}

double FastQSketch::dffunc(double w) {
    double res = 0;
    double e = 2.718282;
    for (size_t i = 0; i < m_; ++i) {
        double x = pow(2.0, -M_[i] - 1);
        double ex = pow(e, w * x);
        res += -x * x * ex * pow(ex - 1, -2);
    }
    return res;
}

double FastQSketch::Newton(double c0) {
    std::cout << "ESTIMATE TIME:" << std::endl;
    // print_vector(this->M_);
    double err = 1e-5;
    double c1 = c0 - ffunc(c0) / dffunc(c0);
    int it = 0;
    while (abs (c1 - c0) > err) {
        // std::cout << "i am doing iteration for estimation " << it << "and c0=" << c0 << ", c1=" << c1 << std::endl;

        c0 = c1;
        c1 = c0 - (double)ffunc(c0) / dffunc(c0);
        it += 1;
        if (it > 5){
            break;
        }
    }
    return c1;
}

double FastQSketch::estimate() {
    return Newton(initialValue());
}
