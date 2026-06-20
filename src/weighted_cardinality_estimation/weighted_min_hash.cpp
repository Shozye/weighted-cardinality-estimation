#include "weighted_min_hash.hpp"
#include <cmath>
#include "hash_util.hpp"

WeightedMinHash::WeightedMinHash(std::size_t sketch_size, std::uint64_t master_seed)
    : Sketch(sketch_size, master_seed), M_(sketch_size, 0.0)
{}

WeightedMinHash::WeightedMinHash(std::size_t sketch_size, std::uint64_t master_seed, const std::vector<double>& registers)
    : Sketch(sketch_size, master_seed), M_(registers)
{}

void WeightedMinHash::add(const std::string& elem, double weight)
{
    validate_weight(weight);
    for (std::size_t i = 0; i < size; ++i) {
        std::uint64_t h = murmur64(elem, seeds_[i]);
        double u = to_unit_interval(h);
        // Beta(w,1) transform: u^(1/w) ~ Beta(w,1)
        double beta_sample = std::pow(u, 1.0 / weight);
        if (beta_sample > M_[i]) {
            M_[i] = beta_sample;
        }
    }
}

double WeightedMinHash::estimate() const
{
    double sum = 0.0;
    for (double v : M_) { sum += (1.0 - v); }
    // Method-1 estimator (m independent hash functions, Algorithm 2):
    // h+_k ~ Beta(w,1), E[1-h+_k] = 1/(w+1), MLE: w_hat = (m-1)/sum(1-h+_k)
    return static_cast<double>(size - 1) / sum;
}

size_t WeightedMinHash::memory_usage(uint64_t flags) const {
    uint64_t f = resolve_flags(flags);
    size_t s = 0;
    if (f & MemoryFlag::REGISTERS) s += M_.capacity() * sizeof(double);
    if (f & MemoryFlag::ALL_WRITE_NO_REGISTERS) s += sizeof(size);
    if (f & MemoryFlag::SEEDS) s += seeds_.bytes();
    return s;
}

void WeightedMinHash::merge(const WeightedMinHash& other) {
    if (other.size != size) { throw std::invalid_argument("Cannot merge sketches of different sizes."); }
    for (std::size_t i = 0; i < size; ++i) {
        M_[i] = std::max(M_[i], other.M_[i]);
    }
}
