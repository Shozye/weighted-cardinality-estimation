#include "q_sketch.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "hash_util.hpp"
#include<cstring>
#include "utils.hpp"

QSketch::QSketch(std::size_t sketch_size, std::uint64_t master_seed, uint8_t amount_bits, RngEngine engine)
    : Sketch(sketch_size, master_seed),
      fisher_yates(size, engine),
      amount_bits_(amount_bits),
      r_max((1 << (amount_bits - 1)) - 1),
      r_min(-(1 << (amount_bits - 1)) + 1),
      M_(amount_bits, sketch_size),
      j_star(0)
{
    if (amount_bits < 2) { throw std::invalid_argument("Amount of bits 'b' must be >= 2."); }
    for (std::size_t i = 0; i < size; ++i) {
        M_[i] = r_min;
    }
}

QSketch::QSketch(std::size_t sketch_size, std::uint64_t master_seed, std::uint8_t amount_bits, const std::vector<int>& registers, RngEngine engine)
    : Sketch(sketch_size, master_seed),
      fisher_yates(size, engine),
      amount_bits_(amount_bits),
      r_max((1 << (amount_bits - 1)) - 1),
      r_min(-(1 << (amount_bits - 1)) + 1),
      M_(amount_bits, sketch_size),
      j_star(argmin(registers))
{
    if (amount_bits < 2) { throw std::invalid_argument("Amount of bits 'b' must be >= 2."); }
    if (registers.size() != sketch_size) { throw std::invalid_argument("Invalid state: registers vector size mismatch"); }
    for (std::size_t i = 0; i < size; ++i) {
        M_[i] = registers[i];
    }
}


size_t QSketch::memory_usage(uint64_t flags) const {
    uint64_t f = resolve_flags(flags);
    size_t s = 0;
    if (f & MemoryFlag::REGISTERS) s += M_.bytes();
    if (f & MemoryFlag::ALL_WRITE_NO_REGISTERS) s += sizeof(size) + sizeof(amount_bits_) + sizeof(r_max) + sizeof(r_min) + sizeof(j_star);
    s += fisher_yates.memory_usage(f);
    if (f & MemoryFlag::SEEDS) s += seeds_.bytes();
    return s;
}

std::uint8_t QSketch::get_amount_bits() const { return amount_bits_; }
std::vector<int> QSketch::get_registers() const {
    return std::vector<int>(M_.begin(), M_.end());
}

void QSketch::add(const std::string& elem, double weight){ 
    validate_weight(weight);
    double r = 0;

    fisher_yates.initialize(elem); 
    for (size_t k = 0; k < this->size; ++k){
        std::uint64_t hashed = murmur64(elem, seeds_[k]); 
        double unit_interval_hash = to_unit_interval(hashed); 
        r -= (std::log(unit_interval_hash) / (weight*(double)(size - k))); 
        int y = static_cast<int>(std::floor(-std::log2(r)));

        if ( y <= M_[j_star] ) { break; } 

        auto j = fisher_yates.get_fisher_yates_element(k);

        if (y > this->M_[j]){
            M_[j] = std::min(std::max(y, r_min), r_max);
            if (j == j_star){
                j_star = argmin(M_);
            }
        }
    }

} 

double QSketch::initialValue() const {
    double tmp_sum = 0.0;
    for(int r: M_) { 
        tmp_sum += std::ldexp(1.0, -r);
    }
    return (double)(this->size-1) / tmp_sum;
}

double QSketch::ffunc_divided_by_dffunc(double w) const {
    double ffunc = 0;
    double dffunc = 0;
    for (int r: M_) {
        double x = std::ldexp(1.0, -r - 1);
        double ex = std::exp(w * x);
        ffunc += x * (2.0 - ex) / (ex - 1.0);
        dffunc += -x * x * ex * pow(ex - 1, -2);
    }
    return ffunc / dffunc;
}

double QSketch::Newton(double c0) const {
    double c1 = c0 - ffunc_divided_by_dffunc(c0);
    int it = 0;
    while (std::abs(c1 - c0) > newton_max_error) {
        c0 = c1;
        c1 = c0 - ffunc_divided_by_dffunc(c0);
        it += 1;
        if (it > newton_max_iterations){ break; }
    }
    return c1;
}

double QSketch::estimate() const {
    return Newton(initialValue());
}

void QSketch::merge(const QSketch& other) {
    if (other.size != size) { throw std::invalid_argument("Cannot merge sketches of different sizes."); }
    for (std::size_t i = 0; i < size; ++i) {
        M_[i] = std::max((int)M_[i], (int)other.M_[i]);
    }
    j_star = argmin(M_);
}
