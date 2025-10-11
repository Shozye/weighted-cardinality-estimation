#include "q_sketch_dyn.hpp"
#include <stdexcept>
#include <cmath>
#include <numeric> 
#include "hash_util.hpp"

QSketchDyn::QSketchDyn(std::size_t m, const std::vector<std::uint32_t>& seeds, std::uint8_t amount_bits, std::uint32_t g_seed)
    : m_(m),
      amount_bits_(amount_bits),
      r_min(-(1 << (amount_bits - 1))),
      r_max((1 << (amount_bits - 1)) - 1),
      seeds_(seeds),
      g_seed_(g_seed),
      k_idx_(1 << amount_bits),
      cardinality_(0.0),
      q_r_(0.0),
      R_(amount_bits, m),
      T_(1 << amount_bits, 0)
{
    if (m == 0) { throw std::invalid_argument("Sketch size 'm' must be positive."); }
    if (amount_bits == 0) { throw std::invalid_argument("Amount of bits 'b' must be positive."); }
    if (seeds.size() != m) { throw std::invalid_argument("Seeds vector must have length m"); }

    if (!T_.empty()) {
        T_[0] = static_cast<int>(m);
    }

    for (std::size_t i = 0; i < m_; ++i) {
        R_[i] = r_min;
    }
    
    std::iota(k_idx_.begin(), k_idx_.end(), 0);
}

QSketchDyn::QSketchDyn(
    std::size_t m, std::uint8_t amount_bits, std::uint32_t g_seed,
    const std::vector<std::uint32_t>& seeds, const std::vector<int>& registers,
    const std::vector<int>& t_histogram, double cardinality)
    : m_(m),
      amount_bits_(amount_bits),
      r_min(-(1 << (amount_bits - 1))),
      r_max((1 << (amount_bits - 1)) - 1),
      seeds_(seeds),
      g_seed_(g_seed),
      k_idx_(1 << amount_bits),
      cardinality_(cardinality),
      q_r_(0.0),
      R_(amount_bits, m),
      T_(t_histogram)
{
    std::iota(k_idx_.begin(), k_idx_.end(), 0);
    for (std::size_t i = 0; i < m_; ++i) {
        R_[i] = registers[i];
    }
}

void QSketchDyn::add(const std::string& elem, double weight) {
    const uint64_t g_hash = murmur64(elem, g_seed_, hash_answer);
    const size_t j = g_hash % m_;

    const uint64_t u_hash = murmur64(elem, seeds_[j], hash_answer);
    const double u = to_unit_interval(u_hash);
    if (u == 0.0) { return; }
    const double r = -std::log(u) / weight;
    const int y = static_cast<int>(std::floor(-std::log2(r)));

    if (y <= R_[j]) {
        return;
    }

    const int old_r_val = R_[j];
    const int new_r_val = std::min(y, r_max);

    const size_t old_t_idx = std::max(0, old_r_val - r_min);
    const size_t new_t_idx = std::max(0, new_r_val - r_min);
    
    if (old_t_idx < T_.size() && T_[old_t_idx] > 0) {
        T_[old_t_idx]--;
    }
    if (new_t_idx < T_.size()) {
        T_[new_t_idx]++;
    }
    R_[j] = new_r_val;

    double q_val_sum = 0.0;
    for (size_t k = 0; k < T_.size(); ++k) {
        if (T_[k] > 0) {
            double exponent = -(static_cast<double>(k_idx_[k]) + r_min + 1.0);
            double exp2_val = std::ldexp(1.0, static_cast<int>(exponent));
            double rate = std::exp(-weight * exp2_val);
            q_val_sum += static_cast<double>(T_[k]) * rate;
        }
    }

    this->q_r_ = 1.0 - (q_val_sum / (double)m_);
    cardinality_ += weight / this->q_r_;
}

void QSketchDyn::add_many(const std::vector<std::string>& elems, const std::vector<double>& weights) {
    if (elems.size() != weights.size()) {
        throw std::invalid_argument("add_many: elems and weights size mismatch");
    }
    for (size_t i = 0; i < elems.size(); ++i) {
        this->add(elems[i], weights[i]);
    }
}

double QSketchDyn::estimate() const {
    return cardinality_;
}

std::size_t QSketchDyn::get_m() const { return m_; }
std::uint8_t QSketchDyn::get_amount_bits() const { return amount_bits_; }
std::uint32_t QSketchDyn::get_g_seed() const { return g_seed_; }
std::vector<std::uint32_t> QSketchDyn::get_seeds() const { return seeds_.toVector(); }
std::vector<int> QSketchDyn::get_registers() const {
    return std::vector<int>(R_.begin(), R_.end());
}
const std::vector<int>& QSketchDyn::get_t_histogram() const { return T_; }
double QSketchDyn::get_cardinality() const { return cardinality_; }

size_t QSketchDyn::memory_usage_total() const {
    size_t total = 0;
    total += sizeof(m_);
    total += sizeof(amount_bits_);
    total += sizeof(r_min);
    total += sizeof(r_max);
    total += sizeof(g_seed_);
    total += sizeof(cardinality_);
    total += sizeof(q_r_);
    total += seeds_.bytes();
    total += k_idx_.capacity() * sizeof(int);
    total += R_.bytes();
    total += T_.capacity() * sizeof(int);
    return total;
}

size_t QSketchDyn::memory_usage_write() const {
    size_t write_size = 0;
    write_size += sizeof(cardinality_);
    write_size += sizeof(q_r_);
    write_size += R_.bytes();
    write_size += T_.capacity() * sizeof(int);
    return write_size;
}

size_t QSketchDyn::memory_usage_estimate() const {
    return R_.bytes() + (T_.capacity() * sizeof(int));
}
