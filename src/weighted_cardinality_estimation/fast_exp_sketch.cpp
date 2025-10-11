#include "fast_exp_sketch.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include "hash_util.hpp"
#include <cstring>

FastExpSketch::FastExpSketch(std::size_t sketch_size, const std::vector<std::uint32_t>& seeds)
    : size(sketch_size), 
      seeds_(seeds),
      M_(sketch_size, std::numeric_limits<double>::infinity()),
      fisher_yates(sketch_size),
      max(std::numeric_limits<double>::infinity())
{
    if (sketch_size == 0) { throw std::invalid_argument("Sketch size 'm' must be positive."); }
    if ((!seeds.empty() && seeds.size() != size)) { 
        throw std::invalid_argument("Seeds must have length m or 0"); 
    }
}



void FastExpSketch::add(const std::string& elem, double weight)
{ 
    double S = 0;
    bool updateMax = false; 

    fisher_yates.initialize(murmur64(elem, 1, hash_answer)); 
    

    auto inv_weight = 1.0 / weight;
    for (size_t k = 0; k < this->size; ++k){
        std::uint64_t hashed = murmur64(elem, seeds_[k], hash_answer); 
        double U = to_unit_interval(hashed); 
        double E = -std::log(U) * inv_weight; 

        S += E/(double)(this->size-k); 
        if ( S >= this->max ) { break; }

        auto j = fisher_yates.get_fisher_yates_element(k);

        if (this->M_[j] == this->max ) { updateMax = true; }
        this->M_[j] = std::min(this->M_[j], S);
    }

    if(updateMax){
        this->max = *std::max_element(this->M_.begin(), this->M_.end());
    }
} 

size_t FastExpSketch::memory_usage_total() const {
    size_t total_size = 0;
    total_size += sizeof(this->size);
    total_size += fisher_yates.bytes_total();
    total_size += sizeof(max);
    total_size += seeds_.bytes();
    total_size += M_.capacity() * sizeof(double);
    return total_size;
}

size_t FastExpSketch::memory_usage_write() const {
    size_t write_size = 0;
    write_size += sizeof(max);
    write_size += fisher_yates.bytes_write();
    write_size += M_.capacity() * sizeof(double);
    return write_size;
}

size_t FastExpSketch::memory_usage_estimate() const {
    size_t estimate_size = M_.capacity() * sizeof(double);
    return estimate_size;
}


void FastExpSketch::add_many(const std::vector<std::string>& elems,
                                  const std::vector<double>& weights) {
    if (elems.size() != weights.size()){
        throw std::invalid_argument("add_many: elems and weights size mismatch");
    }
    for (std::size_t i = 0; i < elems.size(); ++i) {
        this->add(elems[i], weights[i]);
    }
}

double FastExpSketch::estimate() const
{
    double total = 0.0;
    for (double val : M_) { total += val;}
    return ((double)this->size - 1.0) / total;
}

double FastExpSketch::jaccard_struct(const FastExpSketch& other) const
{
    if (other.size != this->size) { return 0.0; }
    std::size_t equal = 0;
    for (std::size_t i = 0; i < this->size; ++i) {
        if (M_[i] == other.M_[i]) { ++equal; } 
    }
    return static_cast<double>(equal) / static_cast<double>(this->size);
}

FastExpSketch::FastExpSketch(
    std::size_t sketch_size,
    const std::vector<std::uint32_t>& seeds,
    const std::vector<double>& registers)
:   size(sketch_size),
    seeds_(seeds),
    M_(registers),
    fisher_yates(sketch_size)
{
    if (size == 0) { throw std::invalid_argument("Sketch size 'm' must be positive."); }
    max = *std::max_element(M_.begin(), M_.end());
}
std::size_t FastExpSketch::get_sketch_size() const { return this->size; }
std::vector<std::uint32_t> FastExpSketch::get_seeds() const { return seeds_.toVector(); }
const std::vector<double>& FastExpSketch::get_registers() const { return M_; }
