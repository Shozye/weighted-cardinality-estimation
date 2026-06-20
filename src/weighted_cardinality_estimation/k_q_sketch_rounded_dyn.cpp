#include "k_q_sketch_rounded_dyn.hpp"
#include <stdexcept>
#include <cmath>
#include "hash_util.hpp"

kQSketchRoundedDyn::kQSketchRoundedDyn(
    std::size_t sketch_size, std::uint64_t master_seed,
    std::uint8_t amount_bits, float logarithm_base, std::uint32_t g_seed)
    : Sketch(sketch_size, master_seed),
      amount_bits_(amount_bits),
      logarithm_base_(logarithm_base),
      r_min(-(1 << (amount_bits - 1)) + 1),
      r_max((1 << (amount_bits - 1)) - 1),
      g_seed_(g_seed),
      cardinality_(0.0),
      q_r_(0.0),
      R_(amount_bits, sketch_size),
      T_(std::ceil(std::log2(sketch_size)), 1 << amount_bits)
{
    if (amount_bits < 2) { throw std::invalid_argument("Amount of bits 'b' must be >= 2."); }
    for (size_t i = 0; i < T_.size(); i++) { T_[i] = 0; }
    T_[0] = sketch_size;
    for (std::size_t i = 0; i < size; ++i) { R_[i] = r_min; }
}

kQSketchRoundedDyn::kQSketchRoundedDyn(
    std::size_t sketch_size, std::uint8_t amount_bits, float logarithm_base,
    std::uint32_t g_seed, std::uint64_t master_seed,
    const std::vector<int>& registers,
    const std::vector<std::uint32_t>& t_histogram, double cardinality)
    : Sketch(sketch_size, master_seed),
      amount_bits_(amount_bits),
      logarithm_base_(logarithm_base),
      r_min(-(1 << (amount_bits - 1)) + 1),
      r_max((1 << (amount_bits - 1)) - 1),
      g_seed_(g_seed),
      cardinality_(cardinality),
      q_r_(0.0),
      R_(amount_bits, sketch_size),
      T_(std::ceil(std::log2(sketch_size)), 1 << amount_bits)
{
    if (amount_bits < 2) { throw std::invalid_argument("Amount of bits 'b' must be >= 2."); }
    for (std::size_t i = 0; i < t_histogram.size(); ++i) { T_[i] = t_histogram[i]; }
    for (std::size_t i = 0; i < size; ++i) { R_[i] = registers[i]; }
}

void kQSketchRoundedDyn::add(const std::string& elem, double weight) {
    validate_weight(weight);
    const uint64_t g_hash = murmur64(elem, g_seed_);
    const size_t j = g_hash % size;

    const uint64_t u_hash = murmur64(elem, seeds_[j]);
    const double u = to_unit_interval(u_hash);
    if (u == 0.0) { return; }
    const double r = -std::log(u) / weight;

    // kQSketchRounding quantization: round(-log_k(r))
    const int y = static_cast<int>(std::round(-std::log(r) / std::log(logarithm_base_)));

    if (y <= R_[j]) { return; }

    // Compute q_r from PRE-update state (correct martingale property)
    double q_val_sum = 0.0;
    for (size_t k = 0; k < T_.size(); ++k) {
        if (T_[k] > 0) {
            double v = static_cast<double>(k) + r_min;
            double threshold = std::pow(logarithm_base_, -(v + 0.5));
            q_val_sum += static_cast<double>(T_[k]) * std::exp(-weight * threshold);
        }
    }
    this->q_r_ = 1.0 - (q_val_sum / (double)size);
    cardinality_ += weight / this->q_r_;

    // Now update registers and histogram
    const int old_r_val = R_[j];
    const int new_r_val = std::min(y, r_max);

    const size_t old_t_idx = std::max(0, old_r_val - r_min);
    const size_t new_t_idx = std::max(0, new_r_val - r_min);

    if (old_t_idx < T_.size() && T_[old_t_idx] > 0) { T_[old_t_idx] = T_[old_t_idx] - 1; }
    if (new_t_idx < T_.size()) { T_[new_t_idx] = T_[new_t_idx] + 1; }
    R_[j] = new_r_val;
}

double kQSketchRoundedDyn::estimate() const { return cardinality_; }
double kQSketchRoundedDyn::estimate_direct() const { return cardinality_; }
double kQSketchRoundedDyn::estimate_newton_cold() const { return cardinality_; }
double kQSketchRoundedDyn::estimate_newton_warm() const { return cardinality_; }
int kQSketchRoundedDyn::estimate_newton_cold_iterations() const { return 0; }
int kQSketchRoundedDyn::estimate_newton_warm_iterations() const { return 0; }

std::uint8_t kQSketchRoundedDyn::get_amount_bits() const { return amount_bits_; }
float kQSketchRoundedDyn::get_logarithm_base() const { return logarithm_base_; }
std::uint32_t kQSketchRoundedDyn::get_g_seed() const { return g_seed_; }
std::vector<int> kQSketchRoundedDyn::get_registers() const {
    return std::vector<int>(R_.begin(), R_.end());
}
std::vector<std::uint32_t> kQSketchRoundedDyn::get_t_histogram() const {
    return std::vector<std::uint32_t>(T_.begin(), T_.end());
}
double kQSketchRoundedDyn::get_cardinality() const { return cardinality_; }

size_t kQSketchRoundedDyn::memory_usage(uint64_t flags) const {
    uint64_t f = resolve_flags(flags);
    size_t s = 0;
    if (f & MemoryFlag::REGISTERS) s += R_.bytes() + T_.bytes();
    if (f & MemoryFlag::ALL_WRITE_NO_REGISTERS) s += sizeof(size) + sizeof(amount_bits_) + sizeof(logarithm_base_) + sizeof(r_min) + sizeof(r_max) + sizeof(g_seed_) + sizeof(cardinality_) + sizeof(q_r_);
    if (f & MemoryFlag::SEEDS) s += seeds_.bytes();
    return s;
}
