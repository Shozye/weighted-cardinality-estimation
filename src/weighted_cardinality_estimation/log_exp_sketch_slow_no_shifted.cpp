#include "log_exp_sketch_slow_no_shifted.hpp"
#include "hash_util.hpp"
#include "quantize_custom_float.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

LogExpSketchSlowNoShifted::LogExpSketchSlowNoShifted(
    std::size_t sketch_size,
    std::uint64_t master_seed,
    std::uint8_t amount_bits,
    double v_max
) : Sketch(sketch_size, master_seed),
    amount_bits_(amount_bits),
    v_max_(v_max),
    r_max_((1 << amount_bits) - 1),
    min_value_(1.0 / v_max),
    log_r_(std::log(v_max / min_value_) / r_max_),
    M_(amount_bits, sketch_size)
{
    if (amount_bits < 2) {
        throw std::invalid_argument("amount_bits must be >= 2.");
    }
    if (v_max <= 1.0) {
        throw std::invalid_argument("v_max must be > 1.");
    }
    for (std::size_t i = 0; i < size; ++i) {
        M_[i] = r_max_;
    }
}

LogExpSketchSlowNoShifted::LogExpSketchSlowNoShifted(
    std::size_t sketch_size,
    std::uint64_t master_seed,
    std::uint8_t amount_bits,
    double v_max,
    const std::vector<int>& registers
) : Sketch(sketch_size, master_seed),
    amount_bits_(amount_bits),
    v_max_(v_max),
    r_max_((1 << amount_bits) - 1),
    min_value_(1.0 / v_max),
    log_r_(std::log(v_max / min_value_) / r_max_),
    M_(amount_bits, sketch_size)
{
    if (amount_bits < 2) {
        throw std::invalid_argument("amount_bits must be >= 2.");
    }
    for (std::size_t i = 0; i < size; ++i) {
        M_[i] = static_cast<unsigned>(registers[i]);
    }
}

int LogExpSketchSlowNoShifted::quantize(double value) const {
    if (value <= min_value_) { return 0; }
    int idx = static_cast<int>(std::round(std::log(value / min_value_) / log_r_));
    return std::min(idx, r_max_);
}

double LogExpSketchSlowNoShifted::reconstruct(int index) const {
    return min_value_ * std::exp(index * log_r_);
}

void LogExpSketchSlowNoShifted::add(const std::string& elem, double weight) {
    validate_weight(weight);
    for (std::size_t i = 0; i < size; ++i) {
        std::uint64_t h = murmur64(elem, seeds_[i]);
        double u = to_unit_interval(h);
        double g = -std::log(u) / weight;
        int idx = quantize(g);
        if (idx < static_cast<int>(M_[i])) {
            M_[i] = static_cast<unsigned>(idx);
        }
    }
}

double LogExpSketchSlowNoShifted::estimate() const {
    double total = 0.0;
    for (std::size_t i = 0; i < size; ++i) {
        total += reconstruct(static_cast<int>(M_[i]));
    }
    return (static_cast<double>(size) - 1.0) / total;
}

double LogExpSketchSlowNoShifted::jaccard_struct(const LogExpSketchSlowNoShifted& other) const {
    if (other.size != size) { return 0.0; }
    std::size_t equal = 0;
    for (std::size_t i = 0; i < size; ++i) {
        if (M_[i] == other.M_[i]) { ++equal; }
    }
    return static_cast<double>(equal) / static_cast<double>(size);
}

void LogExpSketchSlowNoShifted::merge(const LogExpSketchSlowNoShifted& other) {
    if (other.size != size) {
        throw std::invalid_argument("Cannot merge sketches of different sizes.");
    }
    if (other.amount_bits_ != amount_bits_ || other.v_max_ != v_max_) {
        throw std::invalid_argument("Cannot merge sketches with different LNS parameters.");
    }
    for (std::size_t i = 0; i < size; ++i) {
        if (other.M_[i] < M_[i]) {
            M_[i] = static_cast<unsigned>(other.M_[i]);
        }
    }
}

std::vector<int> LogExpSketchSlowNoShifted::get_registers() const {
    return std::vector<int>(M_.begin(), M_.end());
}

std::uint8_t LogExpSketchSlowNoShifted::get_amount_bits() const { return amount_bits_; }
double LogExpSketchSlowNoShifted::get_v_max() const { return v_max_; }

size_t LogExpSketchSlowNoShifted::memory_usage(uint64_t flags) const {
    uint64_t f = resolve_flags(flags);
    size_t s = 0;
    if (f & MemoryFlag::REGISTERS) s += M_.bytes();
    if (f & MemoryFlag::ALL_WRITE_NO_REGISTERS) s += sizeof(size) + sizeof(amount_bits_) + sizeof(v_max_) + sizeof(r_max_) + sizeof(min_value_) + sizeof(log_r_);
    if (f & MemoryFlag::SEEDS) s += seeds_.bytes();
    return s;
}
