#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "hash_util.hpp"
#include "fast_k_q_sketch_rounding.hpp"

kQSketchRounding::kQSketchRounding(
    std::size_t sketch_size,
    std::uint64_t master_seed,
    uint8_t amount_bits,
    float logarithm_base,
    RngEngine engine
)
    : Sketch(sketch_size, master_seed),
      fisher_yates(size, engine),
      amount_bits_(amount_bits),
      logarithm_base(logarithm_base),
      r_max((1 << (amount_bits - 1)) - 1),
      r_min(-(1 << (amount_bits - 1)) + 1),
      M_(amount_bits, sketch_size)
{
    if (amount_bits < 2) { throw std::invalid_argument("Amount of bits 'b' must be >= 2."); }
    std::fill(M_.begin(), M_.end(), r_min);
    update_treshold();
}

kQSketchRounding::kQSketchRounding(
    std::size_t sketch_size,
    std::uint64_t master_seed,
    std::uint8_t amount_bits,
    float logarithm_base,
    const std::vector<int>& registers,
    RngEngine engine
)
    : Sketch(sketch_size, master_seed),
      fisher_yates(size, engine),
      amount_bits_(amount_bits),
      logarithm_base(logarithm_base),
      r_max((1 << (amount_bits - 1)) - 1),
      r_min(-(1 << (amount_bits - 1)) + 1),
      M_(amount_bits, sketch_size)
{
    if (amount_bits < 2) { throw std::invalid_argument("Amount of bits 'b' must be >= 2."); }
    if (registers.size() != sketch_size) { throw std::invalid_argument("Invalid state: registers vector size mismatch"); }
    for (std::size_t i = 0; i < size; ++i) { M_[i] = registers[i]; }
    update_treshold();
}

size_t kQSketchRounding::memory_usage(uint64_t flags) const {
    uint64_t f = resolve_flags(flags);
    size_t s = 0;
    if (f & MemoryFlag::REGISTERS) s += M_.bytes();
    if (f & MemoryFlag::ALL_WRITE_NO_REGISTERS) s += sizeof(size) + sizeof(amount_bits_) + sizeof(r_max) + sizeof(r_min) + sizeof(logarithm_base) + sizeof(min_sketch_value) + sizeof(min_value_to_change_sketch);
    s += fisher_yates.memory_usage(f);
    if (f & MemoryFlag::SEEDS) s += seeds_.bytes();
    return s;
}

std::uint8_t kQSketchRounding::get_amount_bits() const { return amount_bits_; }
float kQSketchRounding::get_logarithm_base() const { return logarithm_base; }
std::vector<int> kQSketchRounding::get_registers() const {
    return std::vector<int>(M_.begin(), M_.end());
}

void kQSketchRounding::update_treshold() {
    this->min_sketch_value = *std::min_element(this->M_.begin(), this->M_.end());
    this->min_value_to_change_sketch = std::pow(logarithm_base, -this->min_sketch_value);
}

void kQSketchRounding::add(const std::string& elem, double weight) {
    validate_weight(weight);
    double S = 0;
    bool touched_min = false;

    fisher_yates.initialize(elem);
    for (size_t k = 0; k < this->size; ++k) {
        std::uint64_t hashed = murmur64(elem, seeds_[k]);
        double unit_interval_hash = to_unit_interval(hashed);
        double exponential_variable = -std::log(unit_interval_hash) / weight;
        S += exponential_variable / (double)(this->size - k);

        if (S >= this->min_value_to_change_sketch) { break; }

        auto j = fisher_yates.get_fisher_yates_element(k);
        int q = static_cast<int>(std::round(-std::log(S) / std::log(logarithm_base)));

        q = std::min(q, r_max);
        if (q > this->M_[j]) {
            if (this->M_[j] == min_sketch_value) { touched_min = true; }
            this->M_[j] = q;
        }
    }

    if (touched_min) { this->update_treshold(); }
}

double kQSketchRounding::estimate_direct() const {
    double tmp_sum = 0.0;
    for (int r : M_) { tmp_sum += std::pow(logarithm_base, -r); }
    return (double)(this->size - 1) / tmp_sum;
}

double kQSketchRounding::estimate() const {
    return estimate_direct();
}

double kQSketchRounding::estimate_corrected() const {
    const double lnk = std::log(logarithm_base);
    const double C_k = 2.0 * std::sinh(lnk / 2.0) / lnk;
    return estimate_direct() * C_k;
}

double kQSketchRounding::estimate_newton_cold() const { return estimate_direct(); }
double kQSketchRounding::estimate_newton_warm() const { return estimate_direct(); }
int kQSketchRounding::estimate_newton_cold_iterations() const { return 0; }
int kQSketchRounding::estimate_newton_warm_iterations() const { return 0; }

void kQSketchRounding::merge(const kQSketchRounding& other) {
    if (other.size != size) { throw std::invalid_argument("Cannot merge sketches of different sizes."); }
    for (std::size_t i = 0; i < size; ++i) {
        M_[i] = std::max((int)M_[i], (int)other.M_[i]);
    }
    update_treshold();
}
