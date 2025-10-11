#include "fastgm_exp_sketch.hpp"
#include "fisher_yates.hpp"
#include "hash_util.hpp"
#include <cmath>
#include <cstdint>
#include"utils.hpp"

FastGMExpSketch::FastGMExpSketch(std::size_t sketch_size, const std::vector<std::uint32_t>& seeds):   
    size(sketch_size),
    M_(std::vector<double>(sketch_size, -1)),
    seeds_(seeds),
    fisher_yates(FisherYates(sketch_size)),
    j_star(1),
    k_star(sketch_size),
    flagFastPrune(false)
{
} 

FastGMExpSketch::FastGMExpSketch(
    std::size_t sketch_size, 
    const std::vector<std::uint32_t>& seeds, 
    const std::vector<double>& registers
):   
    size(sketch_size),
    M_(registers),
    seeds_(seeds),
    fisher_yates(FisherYates(sketch_size)),
    j_star(argmax(M_))
    {
        k_star = size;
        for(double elem: M_){
            if (elem >= 0){
                k_star--;
            }
        }
        flagFastPrune = k_star == 0;
    } // TODO:

void FastGMExpSketch::add(const std::string& elem, double weight)
{ 
    // TODO: Get to know why in original paper there is s_vec
    double b = 0;
    fisher_yates.initialize(murmur64(elem, 1)); 

    for(uint32_t l = 1; l < size; ++l){
        std::uint64_t hashed = murmur64(elem, seeds_[l-1]); 
        double U = to_unit_interval(hashed); 
        b = b - (1/weight)*(std::log(U)/(double)(size-l+1));
        uint32_t c = fisher_yates.get_fisher_yates_element(l-1);

        if (!flagFastPrune){
            if (M_[c] < 0){
                M_[c] = b;
                k_star--;
                if(k_star == 0){ 
                    flagFastPrune = true;
                    j_star = argmax(M_);
                }
            } else if(b < M_[c]){
                M_[c] = b;
            }
        } else if (flagFastPrune) {
            if ( b > M_[j_star]){
                break;
            }
            if ( b < M_[c]){
                M_[c] = b;
                if ( c == j_star){
                    j_star = argmax(M_);
                }
            }
        }
    }
} 

size_t FastGMExpSketch::memory_usage_total() const {
    size_t total_size = 0;
    total_size += sizeof(size);
    total_size += sizeof(k_star);
    total_size += sizeof(j_star);
    total_size += M_.capacity() * sizeof(double);
    total_size += seeds_.bytes();
    total_size += fisher_yates.bytes_total();
    total_size += sizeof(flagFastPrune);
    return total_size; 
}

size_t FastGMExpSketch::memory_usage_write() const {
    size_t total_size = 0;
    total_size += sizeof(k_star);
    total_size += sizeof(j_star);
    total_size += M_.capacity() * sizeof(double);
    total_size += fisher_yates.bytes_write();
    total_size += sizeof(flagFastPrune);
    return total_size; 
}

size_t FastGMExpSketch::memory_usage_estimate() const {
    size_t estimate_size = M_.capacity() * sizeof(double);
    return estimate_size; 
}

double FastGMExpSketch::estimate() const {
    double total = 0.0;
    for (double value : M_) { total += value; }
    return ((double)this->size - 1.0) / total;
}

double FastGMExpSketch::jaccard_struct(const FastGMExpSketch& other) const {
    if (other.size != size) { return 0.0; }
    std::size_t equal = 0;
    for (std::size_t i = 0; i < size; ++i) {
        if (M_[i] == other.M_[i]) { ++equal; }
    }
    return static_cast<double>(equal) / static_cast<double>(size);
}

std::size_t FastGMExpSketch::get_sketch_size() const {
    return size;
}

std::vector<std::uint32_t> FastGMExpSketch::get_seeds() const {
    return seeds_.toVector();
}


const std::vector<double>& FastGMExpSketch::get_registers() const {
    return M_;
}

void FastGMExpSketch::add_many(const std::vector<std::string>& elems,
                                  const std::vector<double>& weights) {
    if (elems.size() != weights.size()){
        throw std::invalid_argument("add_many: elems and weights size mismatch");
    }
    for (std::size_t i = 0; i < elems.size(); ++i) {
        this->add(elems[i], weights[i]);
    }
}
