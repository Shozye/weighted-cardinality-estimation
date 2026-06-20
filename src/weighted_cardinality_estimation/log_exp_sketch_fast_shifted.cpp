#include "log_exp_sketch_fast_shifted.hpp"
#include "hash_util.hpp"
#include "quantize_custom_float.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

LogExpSketchFastShifted::LogExpSketchFastShifted(
    std::size_t sketch_size,
    std::uint64_t master_seed,
    std::uint8_t amount_bits,
    double v_max,
    RngEngine engine
) : Sketch(sketch_size, master_seed),
    amount_bits_(amount_bits),
    v_max_(v_max),
    capacity_((1 << amount_bits) - 1),
    offset_(capacity_),
    num_maxed_(static_cast<int>(sketch_size)),
    min_value_(1.0 / v_max),
    log_r_(std::log(v_max / min_value_) / capacity_),
    M_(amount_bits, sketch_size),
    fisher_yates(sketch_size, engine)
{
    if (amount_bits < 2) { throw std::invalid_argument("amount_bits must be >= 2."); }
    if (v_max <= 1.0) { throw std::invalid_argument("v_max must be > 1."); }
    for (std::size_t i = 0; i < size; ++i) {
        M_[i] = capacity_;
    }
}

LogExpSketchFastShifted::LogExpSketchFastShifted(
    std::size_t sketch_size,
    std::uint64_t master_seed,
    std::uint8_t amount_bits,
    double v_max,
    const std::vector<int>& registers,
    int offset,
    RngEngine engine
) : Sketch(sketch_size, master_seed),
    amount_bits_(amount_bits),
    v_max_(v_max),
    capacity_((1 << amount_bits) - 1),
    offset_(offset),
    num_maxed_(0),
    min_value_(1.0 / v_max),
    log_r_(std::log(v_max / min_value_) / capacity_),
    M_(amount_bits, sketch_size),
    fisher_yates(sketch_size, engine)
{
    if (amount_bits < 2) { throw std::invalid_argument("amount_bits must be >= 2."); }
    if (v_max <= 1.0) { throw std::invalid_argument("v_max must be > 1."); }
    if (registers.size() != sketch_size) {
        throw std::invalid_argument("registers vector size mismatch");
    }
    for (std::size_t i = 0; i < size; ++i) {
        M_[i] = static_cast<unsigned>(registers[i]);
        if (static_cast<int>(M_[i]) == capacity_) ++num_maxed_;
    }
}

int LogExpSketchFastShifted::quantize(double value) const {
    if (value <= min_value_) { return 0; }
    return static_cast<int>(std::round(std::log(value / min_value_) / log_r_));
}

double LogExpSketchFastShifted::reconstruct(int abs_index) const {
    return min_value_ * std::exp(abs_index * log_r_);
}

void LogExpSketchFastShifted::shift_down() {
    int min_val = static_cast<int>(*std::min_element(M_.begin(), M_.end()));
    if (min_val <= 0) {
        num_maxed_ = static_cast<int>(std::count(M_.begin(), M_.end(), static_cast<unsigned>(capacity_)));
        return;
    }
    offset_ -= min_val;
    for (std::size_t i = 0; i < size; ++i) {
        M_[i] = static_cast<int>(M_[i]) - min_val;
    }
    num_maxed_ = static_cast<int>(std::count(M_.begin(), M_.end(), static_cast<unsigned>(capacity_)));
}

void LogExpSketchFastShifted::add(const std::string& elem, double weight) {
    validate_weight(weight);
    double S = 0;
    bool triggered_shift = false;

    // Early-exit threshold: if S exceeds this, no register can be updated
    double max_threshold = reconstruct(capacity_ + offset_);

    fisher_yates.initialize(elem);
    for (std::size_t k = 0; k < size; ++k) {
        std::uint64_t hashed = murmur64(elem, seeds_[k]);
        double U = to_unit_interval(hashed);
        double E = -std::log(U) / weight;

        S += E / static_cast<double>(size - k);
        if (S >= max_threshold) {
            // If sentinel registers remain, shift window up to cover S
            if (num_maxed_ > 0) {
                int q_abs = quantize(S);
                int delta = q_abs - (offset_ + capacity_ - 1);
                if (delta > 0) {
                    offset_ += delta;
                    num_maxed_ = 0;
                    for (std::size_t ii = 0; ii < size; ++ii) {
                        if (static_cast<int>(M_[ii]) == capacity_) {
                            ++num_maxed_;
                        } else {
                            int v = static_cast<int>(M_[ii]) - delta;
                            M_[ii] = (v > 0) ? v : 0;
                        }
                    }
                    max_threshold = reconstruct(capacity_ + offset_);
                }
            } else {
                break;
            }
        }

        std::uint32_t j = fisher_yates.get_fisher_yates_element(k);
        int q_abs = quantize(S);
        int rel = q_abs - offset_;

        // Shift down if value underflows the current window
        if (rel < 0) {
            int delta = -rel;
            offset_ -= delta;
            num_maxed_ = 0;
            for (std::size_t ii = 0; ii < size; ++ii) {
                int v = static_cast<int>(M_[ii]) + delta;
                M_[ii] = (v < capacity_) ? v : capacity_;
                if (static_cast<int>(M_[ii]) == capacity_) ++num_maxed_;
            }
            rel = 0;
            max_threshold = reconstruct(capacity_ + offset_);
        }

        if (rel < static_cast<int>(M_[j])) {
            if (static_cast<int>(M_[j]) == capacity_) {
                if (--num_maxed_ == 0) { triggered_shift = true; }
            }
            M_[j] = rel;
        }
    }

    if (triggered_shift) { shift_down(); }
}

double LogExpSketchFastShifted::estimate() const {
    double total = 0.0;
    for (std::size_t i = 0; i < size; ++i) {
        int abs_idx = static_cast<int>(M_[i]) + offset_;
        total += reconstruct(abs_idx);
    }
    return (static_cast<double>(size) - 1.0) / total;
}

double LogExpSketchFastShifted::jaccard_struct(const LogExpSketchFastShifted& other) const {
    if (other.size != size) { return 0.0; }
    std::size_t equal = 0;
    for (std::size_t i = 0; i < size; ++i) {
        int abs_this = static_cast<int>(M_[i]) + offset_;
        int abs_other = static_cast<int>(other.M_[i]) + other.offset_;
        if (abs_this == abs_other) { ++equal; }
    }
    return static_cast<double>(equal) / static_cast<double>(size);
}

void LogExpSketchFastShifted::merge(const LogExpSketchFastShifted& other) {
    if (other.size != size) {
        throw std::invalid_argument("Cannot merge sketches of different sizes.");
    }
    if (other.amount_bits_ != amount_bits_ || other.v_max_ != v_max_) {
        throw std::invalid_argument("Cannot merge sketches with different LNS parameters.");
    }
    int min_abs = std::numeric_limits<int>::max();
    for (std::size_t i = 0; i < size; ++i) {
        int abs_this = static_cast<int>(M_[i]) + offset_;
        int abs_other = static_cast<int>(other.M_[i]) + other.offset_;
        int abs_min = std::min(abs_this, abs_other);
        if (abs_min < min_abs) min_abs = abs_min;
    }
    int new_offset = min_abs;
    num_maxed_ = 0;
    for (std::size_t i = 0; i < size; ++i) {
        int abs_this = static_cast<int>(M_[i]) + offset_;
        int abs_other = static_cast<int>(other.M_[i]) + other.offset_;
        int abs_min = std::min(abs_this, abs_other);
        int rel = abs_min - new_offset;
        M_[i] = (rel < capacity_) ? rel : capacity_;
        if (static_cast<int>(M_[i]) == capacity_) ++num_maxed_;
    }
    offset_ = new_offset;
}

std::vector<int> LogExpSketchFastShifted::get_registers() const {
    return std::vector<int>(M_.begin(), M_.end());
}

std::uint8_t LogExpSketchFastShifted::get_amount_bits() const { return amount_bits_; }
double LogExpSketchFastShifted::get_v_max() const { return v_max_; }
int LogExpSketchFastShifted::get_offset() const { return offset_; }

size_t LogExpSketchFastShifted::memory_usage(uint64_t flags) const {
    uint64_t f = resolve_flags(flags);
    size_t s = 0;
    if (f & MemoryFlag::REGISTERS) s += M_.bytes();
    if (f & MemoryFlag::ALL_WRITE_NO_REGISTERS) s += sizeof(size) + sizeof(amount_bits_) + sizeof(v_max_) + sizeof(capacity_) + sizeof(offset_) + sizeof(num_maxed_) + sizeof(min_value_) + sizeof(log_r_);
    s += fisher_yates.memory_usage(f);
    if (f & MemoryFlag::SEEDS) s += seeds_.bytes();
    return s;
}
