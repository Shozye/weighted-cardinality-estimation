#include "fastgm_exp_sketch.hpp"
#include "fisher_yates.hpp"
#include "hash_util.hpp"
#include <cmath>
#include <cstdint>
#include"utils.hpp"

FastGMExpSketch::FastGMExpSketch(
    std::size_t sketch_size, 
    std::uint64_t master_seed,
    RngEngine engine
):  Sketch(sketch_size, master_seed),
    M_(std::vector<double>(sketch_size, -1)),
    fisher_yates(FisherYates(sketch_size, engine)),
    j_star(1),
    k_star(sketch_size),
    flagFastPrune(false)
{}

FastGMExpSketch::FastGMExpSketch(
    std::size_t sketch_size, 
    std::uint64_t master_seed, 
    const std::vector<double>& registers,
    RngEngine engine
):  Sketch(sketch_size, master_seed),
    M_(registers),
    fisher_yates(FisherYates(sketch_size, engine)),
    j_star(argmax(M_))
    {
        k_star = size;
        for(double elem: M_){
            if (elem >= 0){
                k_star--;
            }
        }
        flagFastPrune = k_star == 0;
    }

void FastGMExpSketch::add(const std::string& elem, double weight)
{ 
    validate_weight(weight);
    // TODO: Get to know why in original paper there is s_vec
    double b = 0;
    fisher_yates.initialize(elem); 

    for(uint32_t t = 0; t < size; ++t){
        std::uint64_t hashed = murmur64(elem, seeds_[t]); 
        double U = to_unit_interval(hashed); 
        b = b - ((1/weight)*(std::log(U)/(double)(size-t)));
        uint32_t c = fisher_yates.get_fisher_yates_element(t);

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

size_t FastGMExpSketch::memory_usage(uint64_t flags) const {
    uint64_t f = resolve_flags(flags);
    size_t s = 0;
    if (f & MemoryFlag::REGISTERS) s += M_.capacity() * sizeof(double);
    if (f & MemoryFlag::ALL_WRITE_NO_REGISTERS) s += sizeof(size) + sizeof(k_star) + sizeof(j_star) + sizeof(flagFastPrune);
    s += fisher_yates.memory_usage(f);
    if (f & MemoryFlag::SEEDS) s += seeds_.bytes();
    return s;
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


const std::vector<double>& FastGMExpSketch::get_registers() const {
    return M_;
}

void FastGMExpSketch::merge(const FastGMExpSketch& other) {
    if (other.size != size) { throw std::invalid_argument("Cannot merge sketches of different sizes."); }
    for (std::size_t i = 0; i < size; ++i) {
        if (M_[i] < 0 && other.M_[i] < 0) { continue; }          // both unfilled
        if (M_[i] < 0) { M_[i] = other.M_[i]; }                   // this unfilled
        else if (other.M_[i] >= 0) { M_[i] = std::min(M_[i], other.M_[i]); }  // both filled
    }
    // Rebuild fast-prune state
    k_star = 0;
    for (std::size_t i = 0; i < size; ++i) { if (M_[i] < 0) { k_star++; } }
    flagFastPrune = (k_star == 0);
    if (flagFastPrune) { j_star = argmax(M_); }
}
