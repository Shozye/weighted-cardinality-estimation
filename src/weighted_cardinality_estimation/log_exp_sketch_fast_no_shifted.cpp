#include "log_exp_sketch_fast_no_shifted.hpp"
#include "quantize_custom_float.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

LogExpSketchFastNoShifted::LogExpSketchFastNoShifted(
    std::size_t sketch_size,
    std::uint64_t master_seed,
    std::uint8_t amount_bits,
    double v_max,
    RngEngine engine
) : Sketch(sketch_size, master_seed),
    amount_bits_(amount_bits),
    v_max_(v_max),
    r_max_((1 << amount_bits) - 1),
    min_value_(1.0 / v_max),
    log_r_(std::log(v_max / min_value_) / r_max_),
    M_(amount_bits, sketch_size),
    fisher_yates(sketch_size, engine),
    max_register_(r_max_)
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

LogExpSketchFastNoShifted::LogExpSketchFastNoShifted(
    std::size_t sketch_size,
    std::uint64_t master_seed,
    std::uint8_t amount_bits,
    double v_max,
    const std::vector<int>& registers,
    RngEngine engine
) : Sketch(sketch_size, master_seed),
    amount_bits_(amount_bits),
    v_max_(v_max),
    r_max_((1 << amount_bits) - 1),
    min_value_(1.0 / v_max),
    log_r_(std::log(v_max / min_value_) / r_max_),
    M_(amount_bits, sketch_size),
    fisher_yates(sketch_size, engine),
    max_register_(0)
{
    for (std::size_t i = 0; i < size; ++i) {
        M_[i] = static_cast<unsigned>(registers[i]);
        if (registers[i] > max_register_) { max_register_ = registers[i]; }
    }
}

int LogExpSketchFastNoShifted::quantize(double value) const {
    if (value <= min_value_) { return 0; }
    int idx = static_cast<int>(std::round(std::log(value / min_value_) / log_r_));
    return std::min(idx, r_max_);
}

double LogExpSketchFastNoShifted::reconstruct(int index) const {
    return min_value_ * std::exp(index * log_r_);
}

void LogExpSketchFastNoShifted::update_max_register() {
    max_register_ = 0;
    for (std::size_t i = 0; i < size; ++i) {
        int val = static_cast<int>(M_[i]);
        if (val > max_register_) { max_register_ = val; }
    }
}

void LogExpSketchFastNoShifted::add(const std::string& elem, double weight) {
    validate_weight(weight);
    double S = 0;
    bool update_max = false;

    double max_threshold = reconstruct(max_register_);

    fisher_yates.initialize(elem);
    for (std::size_t k = 0; k < size; ++k) {
        std::uint64_t hashed = murmur64(elem, seeds_[k]);
        double U = to_unit_interval(hashed);
        double E = -std::log(U) / weight;

        S += E / static_cast<double>(size - k);
        if (S >= max_threshold) { break; }

        std::uint32_t j = fisher_yates.get_fisher_yates_element(k);
        int idx = quantize(S);

        if (idx < static_cast<int>(M_[j])) {
            if (static_cast<int>(M_[j]) == max_register_) { update_max = true; }
            M_[j] = static_cast<unsigned>(idx);
        }
    }

    if (update_max) {
        update_max_register();
    }
}

double LogExpSketchFastNoShifted::estimate() const {
    double total = 0.0;
    for (std::size_t i = 0; i < size; ++i) {
        total += reconstruct(static_cast<int>(M_[i]));
    }
    return (static_cast<double>(size) - 1.0) / total;
}

double LogExpSketchFastNoShifted::jaccard_struct(const LogExpSketchFastNoShifted& other) const {
    if (other.size != size) { return 0.0; }
    std::size_t equal = 0;
    for (std::size_t i = 0; i < size; ++i) {
        if (M_[i] == other.M_[i]) { ++equal; }
    }
    return static_cast<double>(equal) / static_cast<double>(size);
}

void LogExpSketchFastNoShifted::merge(const LogExpSketchFastNoShifted& other) {
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
    update_max_register();
}

std::vector<int> LogExpSketchFastNoShifted::get_registers() const {
    return std::vector<int>(M_.begin(), M_.end());
}

std::uint8_t LogExpSketchFastNoShifted::get_amount_bits() const { return amount_bits_; }
double LogExpSketchFastNoShifted::get_v_max() const { return v_max_; }

size_t LogExpSketchFastNoShifted::memory_usage(uint64_t flags) const {
    uint64_t f = resolve_flags(flags);
    size_t s = 0;
    if (f & MemoryFlag::REGISTERS) s += M_.bytes();
    if (f & MemoryFlag::ALL_WRITE_NO_REGISTERS) s += sizeof(size) + sizeof(amount_bits_) + sizeof(v_max_) + sizeof(r_max_) + sizeof(min_value_) + sizeof(log_r_) + sizeof(max_register_);
    s += fisher_yates.memory_usage(f);
    if (f & MemoryFlag::SEEDS) s += seeds_.bytes();
    return s;
}
