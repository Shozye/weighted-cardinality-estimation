#include "fast_exp_sketch_custom_float.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

FastExpSketchCustomFloat::FastExpSketchCustomFloat(
    std::size_t sketch_size,
    std::uint64_t master_seed,
    int exp_bits, int mant_bits,
    RngEngine engine
) : Sketch(sketch_size, master_seed),
    exp_bits_(exp_bits), mant_bits_(mant_bits),
    M_(sketch_size, custom_float_max(exp_bits, mant_bits, kMode)),
    fisher_yates(sketch_size, engine),
    max_(custom_float_max(exp_bits, mant_bits, kMode)) {}

FastExpSketchCustomFloat::FastExpSketchCustomFloat(
    std::size_t sketch_size,
    std::uint64_t master_seed,
    int exp_bits, int mant_bits,
    const std::vector<double>& registers,
    RngEngine engine
) : Sketch(sketch_size, master_seed),
    exp_bits_(exp_bits), mant_bits_(mant_bits),
    M_(registers),
    fisher_yates(sketch_size, engine),
    max_(*std::max_element(registers.begin(), registers.end())) {}

void FastExpSketchCustomFloat::add(const std::string& elem, double weight) {
    validate_weight(weight);
    double S = 0;
    bool update_max = false;

    fisher_yates.initialize(elem);
    for (std::size_t k = 0; k < size; ++k) {
        std::uint64_t hashed = murmur64(elem, seeds_[k]);
        double U = to_unit_interval(hashed);
        double E = -std::log(U) / weight;

        S += E / static_cast<double>(size - k);
        if (S >= max_) { break; }

        std::uint32_t j = fisher_yates.get_fisher_yates_element(k);
        double quantized = quantize_custom_float(S, 0, exp_bits_, mant_bits_, kMode);

        if (quantized < M_[j]) {
            if (M_[j] == max_) { update_max = true; }
            M_[j] = quantized;
        }
    }

    if (update_max) {
        max_ = *std::max_element(M_.begin(), M_.end());
    }
}

double FastExpSketchCustomFloat::estimate() const {
    double total = 0.0;
    for (double val : M_) { total += val; }
    return (static_cast<double>(size) - 1.0) / total;
}

double FastExpSketchCustomFloat::jaccard_struct(const FastExpSketchCustomFloat& other) const {
    if (other.size != size) { return 0.0; }
    std::size_t equal = 0;
    for (std::size_t i = 0; i < size; ++i) {
        if (M_[i] == other.M_[i]) { ++equal; }
    }
    return static_cast<double>(equal) / static_cast<double>(size);
}

const std::vector<double>& FastExpSketchCustomFloat::get_registers() const { return M_; }

void FastExpSketchCustomFloat::merge(const FastExpSketchCustomFloat& other) {
    if (other.size != size)
        throw std::invalid_argument("Cannot merge sketches of different sizes.");
    if (other.exp_bits_ != exp_bits_ || other.mant_bits_ != mant_bits_)
        throw std::invalid_argument("Cannot merge sketches with different float formats.");
    for (std::size_t i = 0; i < size; ++i)
        M_[i] = std::min(M_[i], other.M_[i]);
    max_ = *std::max_element(M_.begin(), M_.end());
}

FastExpSketchCustomFloat FastExpSketchCustomFloat::clone_with(int exp_bits, int mant_bits) const {
    if (exp_bits > exp_bits_ || mant_bits > mant_bits_)
        throw std::invalid_argument("clone_with: new format must not exceed original precision");
    std::vector<double> new_regs(size);
    for (std::size_t i = 0; i < size; ++i)
        new_regs[i] = quantize_custom_float(M_[i], 0, exp_bits, mant_bits, kMode);
    return FastExpSketchCustomFloat(size, get_master_seed(), exp_bits, mant_bits, new_regs);
}

size_t FastExpSketchCustomFloat::memory_usage(uint64_t flags) const {
    uint64_t f = resolve_flags(flags);
    size_t s = 0;
    if (f & MemoryFlag::REGISTERS) s += M_.capacity() * sizeof(double);
    if (f & MemoryFlag::ALL_WRITE_NO_REGISTERS) s += sizeof(size) + sizeof(max_);
    s += fisher_yates.memory_usage(f);
    if (f & MemoryFlag::SEEDS) s += seeds_.bytes();
    return s;
}
