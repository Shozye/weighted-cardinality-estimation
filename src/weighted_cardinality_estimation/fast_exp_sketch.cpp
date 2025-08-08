#include "fast_exp_sketch.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
#include "hash_util.hpp"
#include<cstring>
#include<iostream>

FastExpSketch::FastExpSketch(std::size_t m, const std::vector<std::uint32_t>& seeds)
    : m_(m), seeds_(seeds),
      M_(m, std::numeric_limits<double>::infinity()),
      permInit(m),
      permWork(m),
      rng_seed(0),
      max(std::numeric_limits<double>::infinity())
{
    if (seeds_.size() != m_)
        throw std::invalid_argument("Seeds vector must have length m");
    for(int i = 0; i < m_; i++){
        permInit[i] = i+1;
    }
}

int FastExpSketch::rand(int min, int max){
    this->rng_seed = this->rng_seed * 1103515245 + 12345;
    auto temp = (unsigned)(this->rng_seed/65536) % 32768;
    return (temp % (max-min)) + min;
}

void FastExpSketch::add(const std::string& x, double weight)
{ 
    std::uint64_t hash_answer[2];
    double S = 0;
    bool updateMax = false; 

    this->rng_seed = murmur64(x, 1, hash_answer); 
    permWork = permInit; 
    int k;
    for (k = 0; k < this->m_; ++k){
        std::uint64_t hashed = murmur64(x, seeds_[k], hash_answer); 
        double U = to_unit_interval(hashed); 
        double E = -std::log(U) / weight; 

        S += E/(double)(this->m_-k); 
        if ( S >= this->max ) break; 

        uint32_t r = rand(k, m_);
        auto swap = permWork[k];
        permWork[k] = permWork[r];
        permWork[r] = swap;
        auto j = permWork[k] - 1;

        if (this->M_[j] == this->max ) updateMax = true;
        if (this->M_[j] > S) this->M_[j] = S;
    }

    if(updateMax){
        this->max = this->M_[0];
        for(int k = 0; k < this->m_; ++k){
            if (this->M_[k] > this->max){
                this->max = this->M_[k];
            }
        }
    }
} 

double FastExpSketch::estimate() const
{
    double sum = 0.0;
    for (double v : M_) sum += v;
    return (m_ - 1.0) / sum;
}

double FastExpSketch::jaccard_struct(const FastExpSketch& other) const
{
    if (other.m_ != m_) return 0.0;
    std::size_t equal = 0;
    for (std::size_t i = 0; i < m_; ++i)
        if (M_[i] == other.M_[i]) ++equal;
    return static_cast<double>(equal) / static_cast<double>(m_);
}
