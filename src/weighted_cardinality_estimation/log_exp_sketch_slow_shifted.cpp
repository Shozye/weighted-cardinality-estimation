#include "log_exp_sketch_slow_shifted.hpp"
#include "hash_util.hpp"
#include "quantize_custom_float.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

LogExpSketchSlowShifted::LogExpSketchSlowShifted(
    std::size_t sketch_size,
    std::uint64_t master_seed,
    std::uint8_t amount_bits,
    double v_max
) : Sketch(sketch_size, master_seed),
    amount_bits_(amount_bits),
    v_max_(v_max),
    capacity_((1 << amount_bits) - 1),
    offset_(0),
    num_maxed_(static_cast<int>(sketch_size)),
    log_r_(std::log(v_max) / capacity_),
    M_(amount_bits, sketch_size)
{
    if (amount_bits < 2) { throw std::invalid_argument("amount_bits must be >= 2."); }
    if (v_max <= 1.0) { throw std::invalid_argument("v_max must be > 1."); }
    // v_max is the *relative* window width: the capacity_ registers span a value
    // ratio of v_max, so the grid slides freely (offset may go negative) and the
    // absolute range is unbounded. v_max only needs to exceed the spread of the m
    // register minima (~m*ln m); smaller v_max means a finer grid for the same b.
    // All registers start at capacity_ (relative); the window slides on the first
    // inserts to cover the data (down when a value underflows, up while sentinels
    // remain when a value overflows).
    for (std::size_t i = 0; i < size; ++i) {
        M_[i] = capacity_;
    }
}

LogExpSketchSlowShifted::LogExpSketchSlowShifted(
    std::size_t sketch_size,
    std::uint64_t master_seed,
    std::uint8_t amount_bits,
    double v_max,
    const std::vector<int>& registers,
    int offset
) : Sketch(sketch_size, master_seed),
    amount_bits_(amount_bits),
    v_max_(v_max),
    capacity_((1 << amount_bits) - 1),
    offset_(offset),
    num_maxed_(0),
    log_r_(std::log(v_max) / capacity_),
    M_(amount_bits, sketch_size)
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

int LogExpSketchSlowShifted::quantize(double value) const {
    // Absolute log grid with reference value 1.0 at index 0; index may be negative.
    // No fixed floor: the sliding offset (not v_max) sets the representable scale.
    return static_cast<int>(std::round(std::log(value) / log_r_));
}

double LogExpSketchSlowShifted::reconstruct(int abs_index) const {
    return std::exp(abs_index * log_r_);
}

void LogExpSketchSlowShifted::shift_down() {
    // In a min-sketch, shift_down subtracts the global min from all registers
    // (the min relative value becomes 0, offset decreases)
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

void LogExpSketchSlowShifted::add(const std::string& elem, double weight) {
    validate_weight(weight);
    bool triggered_shift = false;

    for (std::size_t i = 0; i < size; ++i) {
        std::uint64_t h = murmur64(elem, seeds_[i]);
        double u = to_unit_interval(h);
        double g = -std::log(u) / weight;
        int q_abs = quantize(g);
        int rel = q_abs - offset_;

        // If the quantized value underflows (below current window), shift offset down.
        if (rel < 0) {
            int delta = -rel;
            offset_ -= delta;
            num_maxed_ = 0;
            for (std::size_t j = 0; j < size; ++j) {
                int v = static_cast<int>(M_[j]) + delta;
                M_[j] = (v < capacity_) ? v : capacity_;
                if (static_cast<int>(M_[j]) == capacity_) ++num_maxed_;
            }
            rel = 0;
        }
        // If the value overflows the window while sentinel registers remain, slide the
        // window up to cover it. Set registers stay put; only unset (sentinel) ones are
        // still waiting for a value, so this just positions the window on early inserts.
        else if (rel > capacity_ && num_maxed_ > 0) {
            int delta = rel - capacity_;
            offset_ += delta;
            num_maxed_ = 0;
            for (std::size_t j = 0; j < size; ++j) {
                if (static_cast<int>(M_[j]) == capacity_) {
                    ++num_maxed_;
                } else {
                    int v = static_cast<int>(M_[j]) - delta;
                    M_[j] = (v > 0) ? v : 0;
                }
            }
            rel = capacity_;
        }

        if (rel < static_cast<int>(M_[i])) {
            if (static_cast<int>(M_[i]) == capacity_) {
                if (--num_maxed_ == 0) { triggered_shift = true; }
            }
            M_[i] = rel;
        }
    }

    if (triggered_shift) { shift_down(); }
}

double LogExpSketchSlowShifted::estimate() const {
    double total = 0.0;
    for (std::size_t i = 0; i < size; ++i) {
        int abs_idx = static_cast<int>(M_[i]) + offset_;
        total += reconstruct(abs_idx);
    }
    return (static_cast<double>(size) - 1.0) / total;
}

double LogExpSketchSlowShifted::jaccard_struct(const LogExpSketchSlowShifted& other) const {
    if (other.size != size) { return 0.0; }
    std::size_t equal = 0;
    for (std::size_t i = 0; i < size; ++i) {
        int abs_this = static_cast<int>(M_[i]) + offset_;
        int abs_other = static_cast<int>(other.M_[i]) + other.offset_;
        if (abs_this == abs_other) { ++equal; }
    }
    return static_cast<double>(equal) / static_cast<double>(size);
}

void LogExpSketchSlowShifted::merge(const LogExpSketchSlowShifted& other) {
    if (other.size != size) {
        throw std::invalid_argument("Cannot merge sketches of different sizes.");
    }
    if (other.amount_bits_ != amount_bits_ || other.v_max_ != v_max_) {
        throw std::invalid_argument("Cannot merge sketches with different LNS parameters.");
    }
    // Union = pointwise min of absolute indices
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

std::vector<int> LogExpSketchSlowShifted::get_registers() const {
    return std::vector<int>(M_.begin(), M_.end());
}

std::uint8_t LogExpSketchSlowShifted::get_amount_bits() const { return amount_bits_; }
double LogExpSketchSlowShifted::get_v_max() const { return v_max_; }
int LogExpSketchSlowShifted::get_offset() const { return offset_; }

size_t LogExpSketchSlowShifted::memory_usage(uint64_t flags) const {
    uint64_t f = resolve_flags(flags);
    size_t s = 0;
    if (f & MemoryFlag::REGISTERS) s += M_.bytes();
    if (f & MemoryFlag::ALL_WRITE_NO_REGISTERS) s += sizeof(size) + sizeof(amount_bits_) + sizeof(v_max_) + sizeof(capacity_) + sizeof(offset_) + sizeof(num_maxed_) + sizeof(log_r_);
    if (f & MemoryFlag::SEEDS) s += seeds_.bytes();
    return s;
}
