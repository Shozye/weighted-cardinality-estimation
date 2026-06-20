#include "k_q_sketch_shifted.hpp"
#include "hash_util.hpp"
#include "utils.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

kQSketchShifted::kQSketchShifted(
    std::size_t sketch_size,
    std::uint64_t master_seed,
    std::uint8_t amount_bits,
    float logarithm_base,
    RngEngine engine)
    : Sketch(sketch_size, master_seed),
      fisher_yates(sketch_size, engine),
      amount_bits_(amount_bits),
      logarithm_base(logarithm_base),
      capacity_((1 << amount_bits) - 1),
      offset_(-1022),
      num_zeros_(static_cast<int>(sketch_size)),
      threshold_(std::pow(logarithm_base, 1022)),
      M_(amount_bits, sketch_size)
{
    if (amount_bits < 2) { throw std::invalid_argument("Amount of bits 'b' must be >= 2."); }
    std::fill(M_.begin(), M_.end(), 0);
}

kQSketchShifted::kQSketchShifted(
    std::size_t sketch_size,
    std::uint64_t master_seed,
    std::uint8_t amount_bits,
    float logarithm_base,
    const std::vector<int>& registers,
    int offset,
    RngEngine engine)
    : Sketch(sketch_size, master_seed),
      fisher_yates(sketch_size, engine),
      amount_bits_(amount_bits),
      logarithm_base(logarithm_base),
      capacity_((1 << amount_bits) - 1),
      offset_(offset),
      num_zeros_(0),
      threshold_(std::pow(logarithm_base, -offset)),
      M_(amount_bits, sketch_size)
{
    if (amount_bits < 2) { throw std::invalid_argument("Amount of bits 'b' must be >= 2."); }
    if (registers.size() != sketch_size) {
        throw std::invalid_argument("Invalid state: registers vector size mismatch");
    }
    for (std::size_t i = 0; i < size; ++i) {
        M_[i] = registers[i];
        if (registers[i] == 0) num_zeros_++;
    }
}

void kQSketchShifted::shift_up() {
    int min_val = *std::min_element(M_.begin(), M_.end());
    if (min_val <= 0) {
        num_zeros_ = static_cast<int>(std::count(M_.begin(), M_.end(), 0));
        threshold_ = std::pow(logarithm_base, -offset_);
        return;
    }
    offset_ += min_val;
    for (std::size_t i = 0; i < size; ++i) {
        M_[i] = static_cast<int>(M_[i]) - min_val;
    }
    num_zeros_ = static_cast<int>(std::count(M_.begin(), M_.end(), 0));
    threshold_ = std::pow(logarithm_base, -offset_);
}

void kQSketchShifted::add(const std::string& elem, double weight) {
    validate_weight(weight);
    double S = 0;
    bool triggered_shift = false;

    fisher_yates.initialize(elem);
    for (std::size_t k = 0; k < size; ++k) {
        std::uint64_t h = murmur64(elem, seeds_[k]);
        double u = to_unit_interval(h);
        double g = -std::log(u) / weight;
        S += g / (double)(size - k);

        if (S >= threshold_) { break; }

        auto j = fisher_yates.get_fisher_yates_element(k);
        int q_abs = static_cast<int>(std::floor(-std::log(S) / std::log(logarithm_base)));
        int rel = q_abs - offset_;

        // Shift-on-overflow: if the value doesn't fit, shift offset up immediately
        if (rel > capacity_) {
            int delta = rel - capacity_;
            offset_ += delta;
            threshold_ = std::pow(logarithm_base, -offset_);
            num_zeros_ = 0;
            for (std::size_t i = 0; i < size; ++i) {
                int v = static_cast<int>(M_[i]) - delta;
                M_[i] = (v > 0) ? v : 0;
                if (M_[i] == 0) ++num_zeros_;
            }
            rel = capacity_;
        }

        if (rel > static_cast<int>(M_[j])) {
            if (static_cast<int>(M_[j]) == 0) {
                if (--num_zeros_ == 0) { triggered_shift = true; }
            }
            M_[j] = rel;
        }
    }

    if (triggered_shift) { shift_up(); }
}

// ─── Estimation (absolute values = M_[i] + offset_) ─────────────────────────

double kQSketchShifted::initialValue() const {
    double tmp_sum = 0.0;
    for (int r : M_) {
        tmp_sum += std::pow(logarithm_base, -(r + offset_));
    }
    return (double)(this->size - 1) / tmp_sum;
}

double kQSketchShifted::ffunc_divided_by_dffunc(double w) const {
    double ffunc = 0;
    double dffunc = 0;
    const double k = logarithm_base;
    for (int r : M_) {
        double c = std::pow(k, -(r + offset_));
        double a = std::exp(-w * c);
        double b = std::exp(-w * c / k);
        double diff = b - a;
        ffunc += (c * a - (c / k) * b) / diff;
        dffunc += -(c - c / k) * (c - c / k) * a * b / (diff * diff);
    }
    return ffunc / dffunc;
}

double kQSketchShifted::Newton(double c0) const {
    double c1 = c0 - ffunc_divided_by_dffunc(c0);
    int it = 0;
    while (std::abs(c1 - c0) / std::abs(c1) > newton_max_error) {
        c0 = c1;
        c1 = c0 - ffunc_divided_by_dffunc(c0);
        if (++it > newton_max_iterations) { throw std::runtime_error("Newton-Raphson did not converge within max iterations"); }
    }
    return c1;
}

std::pair<double, int> kQSketchShifted::Newton_with_iterations(double c0) const {
    double c1 = c0 - ffunc_divided_by_dffunc(c0);
    int it = 0;
    while (std::abs(c1 - c0) / std::abs(c1) > newton_max_error) {
        c0 = c1;
        c1 = c0 - ffunc_divided_by_dffunc(c0);
        if (++it > newton_max_iterations) { throw std::runtime_error("Newton-Raphson did not converge within max iterations"); }
    }
    return {c1, it};
}

double kQSketchShifted::estimate() const {
    return Newton(initialValue());
}

double kQSketchShifted::estimate_direct() const {
    double tmp_sum = 0.0;
    for (int r : M_) {
        tmp_sum += std::pow(logarithm_base, -(r + offset_));
    }
    const double m = (double)this->size;
    const double k = logarithm_base;
    return (k - 1) * m / (std::log(k) * tmp_sum);
}

double kQSketchShifted::estimate_newton_cold() const {
    return Newton(std::pow(logarithm_base, offset_));
}

double kQSketchShifted::estimate_newton_warm() const {
    return Newton(estimate_direct());
}

int kQSketchShifted::estimate_newton_cold_iterations() const {
    return Newton_with_iterations(std::pow(logarithm_base, offset_)).second;
}

int kQSketchShifted::estimate_newton_warm_iterations() const {
    return Newton_with_iterations(estimate_direct()).second;
}

// ─── Merge ───────────────────────────────────────────────────────────────────

void kQSketchShifted::merge(const kQSketchShifted& other) {
    if (other.size != size) { throw std::invalid_argument("Cannot merge sketches of different sizes."); }

    int max_abs = std::numeric_limits<int>::min();
    for (std::size_t i = 0; i < size; ++i) {
        int abs_this = static_cast<int>(M_[i]) + offset_;
        int abs_other = static_cast<int>(other.M_[i]) + other.offset_;
        int abs_max = std::max(abs_this, abs_other);
        if (abs_max > max_abs) max_abs = abs_max;
    }
    int new_offset = max_abs - capacity_;
    for (std::size_t i = 0; i < size; ++i) {
        int abs_this = static_cast<int>(M_[i]) + offset_;
        int abs_other = static_cast<int>(other.M_[i]) + other.offset_;
        int abs_max = std::max(abs_this, abs_other);
        M_[i] = abs_max - new_offset;
    }
    offset_ = new_offset;
    shift_up();
}

// ─── Accessors ───────────────────────────────────────────────────────────────

std::uint8_t kQSketchShifted::get_amount_bits() const { return amount_bits_; }
float kQSketchShifted::get_logarithm_base() const { return logarithm_base; }
int kQSketchShifted::get_offset() const { return offset_; }
std::vector<int> kQSketchShifted::get_registers() const {
    return std::vector<int>(M_.begin(), M_.end());
}

size_t kQSketchShifted::memory_usage(uint64_t flags) const {
    uint64_t f = resolve_flags(flags);
    size_t s = 0;
    if (f & MemoryFlag::REGISTERS) s += M_.bytes();
    if (f & MemoryFlag::ALL_WRITE_NO_REGISTERS) s += sizeof(size) + sizeof(amount_bits_) + sizeof(logarithm_base) + sizeof(capacity_) + sizeof(offset_) + sizeof(num_zeros_) + sizeof(threshold_);
    s += fisher_yates.memory_usage(f);
    if (f & MemoryFlag::SEEDS) s += seeds_.bytes();
    return s;
}
