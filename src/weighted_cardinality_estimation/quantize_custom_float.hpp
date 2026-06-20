#pragma once
#include <cmath>
#include <cstdint>
#include <stdexcept>

// ╔══════════════════════════════════════════════════════════════════════════════╗
// ║              CUSTOM FLOATING-POINT QUANTIZATION MODES                       ║
// ║                                                                             ║
// ║  Both modes define an UNSIGNED (sign_bits=0) or SIGNED (sign_bits=1)        ║
// ║  floating-point format with:                                                ║
// ║    - exp_bits bits of exponent                                              ║
// ║    - mant_bits bits of mantissa (fractional part)                           ║
// ║    - NO infinities, NO NaN — all bit patterns encode finite values          ║
// ║    - Truncation rounding (round toward zero)                                ║
// ║    - Saturation on overflow/underflow (clamp to max/min representable)      ║
// ║                                                                             ║
// ║  The format has bias = 2^(exp_bits - 1) - 1  (standard IEEE convention).   ║
// ╚══════════════════════════════════════════════════════════════════════════════╝

enum class QuantizationMode : int {
    ALL_NORMAL = 0,
    WITH_SUBNORMALS = 1,
    LINEAR = 2,
    LOGARITHMIC = 3,
};

// ─────────────────────────────────────────────────────────────────────────────
// MODE 1: ALL-NORMAL (Log-Uniform Distribution)
// ─────────────────────────────────────────────────────────────────────────────
//
// ┌─────────────────────────────────────────────────────────────────────────┐
// │ KEY IDEA: Every exponent pattern (including 0) encodes a NORMAL number  │
// │ with an implicit leading 1. This gives CONSTANT RELATIVE PRECISION      │
// │ across the entire representable range.                                   │
// └─────────────────────────────────────────────────────────────────────────┘
//
// A normal floating-point number has the form:
//
//      value = 1.mmm...m × 2^(biased_exp - bias)
//
// where "1.mmm...m" means an implicit leading 1 followed by mant_bits of
// fractional binary digits.
//
// REPRESENTABLE VALUES:
//   For a given biased exponent e ∈ {0, 1, ..., 2^exp_bits - 1}:
//     value = (1 + k/2^mant_bits) × 2^(e - bias)     for k ∈ {0, ..., 2^mant_bits - 1}
//
// SPACING PROPERTY (why it's "log-uniform"):
//   In the interval [2^n, 2^(n+1)) there are exactly 2^mant_bits representable values.
//   The step size within this interval is:
//       step = 2^n / 2^mant_bits = 2^(n - mant_bits)
//
//   One exponent lower, in [2^(n-1), 2^n), the step size halves:
//       step = 2^(n-1) / 2^mant_bits = 2^(n - 1 - mant_bits)
//
//   The RELATIVE step (step / value) is always ≈ 1/2^mant_bits, independent of
//   magnitude. This means equal precision at every scale — "log-uniform".
//
// RANGE:
//   min = 1.0 × 2^(0 - bias) = 2^(-bias)
//   max = (1 + (2^mant_bits - 1)/2^mant_bits) × 2^(2^exp_bits - 1 - bias)
//       = (2 - 2^(-mant_bits)) × 2^(2^exp_bits - 1 - bias)
//
// TOTAL REPRESENTABLE VALUES (unsigned): 2^exp_bits × 2^mant_bits = 2^(exp_bits + mant_bits)
//   Every single bit pattern maps to a distinct value. Maximum information density.
//
// TRADEOFF vs WITH_SUBNORMALS:
//   ✓ Uniform relative precision everywhere (every value has full mant_bits precision)
//   ✗ Larger minimum value (cannot represent values below 2^(-bias))
//
inline double quantize_all_normal(double value, int sign_bits, int exp_bits, int mant_bits) {
    if (exp_bits < 1 || mant_bits < 0 || sign_bits < 0 || sign_bits > 1) {
        throw std::invalid_argument("Invalid custom float format");
    }

    const int bias = (1 << (exp_bits - 1)) - 1;
    const int max_biased_exp = (1 << exp_bits) - 1;
    const double mant_scale = static_cast<double>(1ULL << mant_bits);
    // min = 2^(0 - bias) = 2^(-bias)  [biased_exp = 0, mantissa = 0, implicit leading 1]
    const double min_value = std::ldexp(1.0, -bias);
    // max = (1 + (2^mant_bits - 1) / 2^mant_bits) × 2^(max_biased_exp - bias)
    const double max_value = std::ldexp(1.0 + (mant_scale - 1.0) / mant_scale, max_biased_exp - bias);

    // Handle special IEEE doubles: clamp to max
    if (std::isnan(value) || std::isinf(value)) {
        return (value < 0.0 && sign_bits) ? -max_value : max_value;
    }

    // Handle sign
    double sign = 1.0;
    if (value < 0.0) {
        if (sign_bits == 0) return min_value;  // unsigned: clamp negatives to min
        sign = -1.0;
        value = -value;
    }

    // Clamp to [min_value, max_value]
    if (value <= 0.0 || value < min_value) return sign * min_value;
    if (value >= max_value) return sign * max_value;

    // Decompose: value = frac × 2^exp_unbiased, where frac ∈ [1.0, 2.0)
    int exp_unbiased;
    double frac = std::frexp(value, &exp_unbiased);
    frac *= 2.0;         // now frac ∈ [1.0, 2.0)
    exp_unbiased -= 1;   // adjust so value = frac × 2^exp_unbiased

    int biased_exp = exp_unbiased + bias;

    // Saturate (should be caught above, but defensive)
    if (biased_exp > max_biased_exp) return sign * max_value;
    if (biased_exp < 0) return sign * min_value;

    // Truncate mantissa to mant_bits (round toward zero)
    // mantissa = frac - 1.0 ∈ [0, 1), multiply by 2^mant_bits and floor
    double mantissa = frac - 1.0;
    double quantized_mant = std::floor(mantissa * mant_scale) / mant_scale;
    return sign * std::ldexp(1.0 + quantized_mant, exp_unbiased);
}

inline double all_normal_max(int exp_bits, int mant_bits) {
    const int bias = (1 << (exp_bits - 1)) - 1;
    const int max_biased_exp = (1 << exp_bits) - 1;
    const double mant_scale = static_cast<double>(1ULL << mant_bits);
    return std::ldexp(1.0 + (mant_scale - 1.0) / mant_scale, max_biased_exp - bias);
}

inline double all_normal_min(int exp_bits, [[maybe_unused]] int mant_bits) {
    const int bias = (1 << (exp_bits - 1)) - 1;
    return std::ldexp(1.0, -bias);
}

// ─────────────────────────────────────────────────────────────────────────────
// MODE 2: WITH SUBNORMALS (IEEE 754-style gradual underflow)
// ─────────────────────────────────────────────────────────────────────────────
//
// ┌─────────────────────────────────────────────────────────────────────────┐
// │ KEY IDEA: When biased_exp = 0, the implicit leading bit becomes 0      │
// │ instead of 1. This creates LINEARLY-SPACED values below the normal     │
// │ range, extending the minimum representable value at the cost of         │
// │ reduced precision for small values.                                     │
// └─────────────────────────────────────────────────────────────────────────┘
//
// NORMAL NUMBERS (biased_exp ≥ 1):
//   Same as all-normal mode:
//      value = (1 + k/2^mant_bits) × 2^(biased_exp - bias)
//
//   Step within [2^n, 2^(n+1)):  step = 2^(n - mant_bits)
//   These have constant relative precision.
//
// SUBNORMAL NUMBERS (biased_exp = 0):
//   The implicit leading bit is 0 (not 1):
//      value = (0 + k/2^mant_bits) × 2^(1 - bias)  =  k × 2^(1 - bias - mant_bits)
//
//   for k ∈ {1, 2, ..., 2^mant_bits - 1}  (k=0 would be zero, which we exclude)
//
//   These are LINEARLY SPACED with constant step:
//      step = 2^(1 - bias - mant_bits)   (= min_subnormal)
//
//   This means:
//     - The smallest subnormal (k=1) has ZERO mantissa bits of precision
//     - The largest subnormal (k=2^mant_bits - 1) has almost mant_bits precision
//     - Relative precision DEGRADES as values approach zero
//
// RANGE:
//   min = 1 × 2^(1 - bias - mant_bits) = 2^(1 - bias - mant_bits)   [smallest subnormal]
//   max = (1 + (2^mant_bits - 1)/2^mant_bits) × 2^(2^exp_bits - 1 - bias)
//       (same max as all-normal since max exponent is unchanged)
//
// TOTAL REPRESENTABLE VALUES (unsigned):
//   Normals: (2^exp_bits - 1) × 2^mant_bits     [biased_exp 1..max_exp]
//   Subnormals: 2^mant_bits - 1                   [k = 1..2^mant_bits - 1]
//   Total: 2^(exp_bits + mant_bits) - 1
//   (One fewer than all-normal because subnormal k=0 is excluded as zero)
//
// TRADEOFF vs ALL_NORMAL:
//   ✓ Smaller minimum value: 2^(1-bias-mant_bits) vs 2^(-bias)
//     (ratio = 2^(mant_bits-1), so with 10 mant bits, min is 512× smaller!)
//   ✗ Non-uniform relative precision in subnormal region (wastes bits)
//   ✗ One fewer distinct representable value (loss of k=0 pattern)
//
// EXAMPLE with exp_bits=3, mant_bits=2, bias=3:
//   Normal step at exponent 0: 2^(0-2) = 0.25
//   Subnormal step (constant): 2^(1-3-2) = 2^(-4) = 0.0625
//
//   Subnormals: 0.0625, 0.125, 0.1875  (k=1,2,3; all with step 0.0625)
//   First normals (biased_exp=1): 0.25, 0.3125, 0.375, 0.4375
//   Note: subnormal 0.0625 has relative precision 0.0625/0.0625 = 100% error
//         while normal 0.25 has relative precision 0.0625/0.25 = 25% error
//
inline double quantize_with_subnormals(double value, int sign_bits, int exp_bits, int mant_bits) {
    if (exp_bits < 1 || mant_bits < 0 || sign_bits < 0 || sign_bits > 1) {
        throw std::invalid_argument("Invalid custom float format");
    }

    const int bias = (1 << (exp_bits - 1)) - 1;
    const int max_biased_exp = (1 << exp_bits) - 1;
    const double mant_scale = static_cast<double>(1ULL << mant_bits);
    // Smallest subnormal: k=1, so value = 1 × 2^(1 - bias - mant_bits)
    const double min_subnormal = std::ldexp(1.0, 1 - bias - mant_bits);
    const double max_value = std::ldexp(1.0 + (mant_scale - 1.0) / mant_scale, max_biased_exp - bias);

    // Handle special IEEE doubles
    if (std::isnan(value) || std::isinf(value)) {
        return (value < 0.0 && sign_bits) ? -max_value : max_value;
    }

    // Handle sign
    double sign = 1.0;
    if (value < 0.0) {
        if (sign_bits == 0) return min_subnormal;  // unsigned: clamp to min
        sign = -1.0;
        value = -value;
    }

    // Clamp to [min_subnormal, max_value]
    if (value <= 0.0 || value < min_subnormal) return sign * min_subnormal;
    if (value >= max_value) return sign * max_value;

    // Decompose: value = frac × 2^exp_unbiased, where frac ∈ [1.0, 2.0)
    int exp_unbiased;
    double frac = std::frexp(value, &exp_unbiased);
    frac *= 2.0;
    exp_unbiased -= 1;

    int biased_exp = exp_unbiased + bias;

    if (biased_exp > max_biased_exp) {
        return sign * max_value;
    }

    if (biased_exp >= 1) {
        // Normal number: truncate mantissa to mant_bits
        double mantissa = frac - 1.0;
        double quantized_mant = std::floor(mantissa * mant_scale) / mant_scale;
        return sign * std::ldexp(1.0 + quantized_mant, exp_unbiased);
    }

    // Subnormal: biased_exp <= 0
    // Quantize to nearest smaller multiple of min_subnormal (truncation)
    double quantized = std::floor(value / min_subnormal) * min_subnormal;
    if (quantized < min_subnormal) return sign * min_subnormal;
    return sign * quantized;
}

inline double with_subnormals_max(int exp_bits, int mant_bits) {
    const int bias = (1 << (exp_bits - 1)) - 1;
    const int max_biased_exp = (1 << exp_bits) - 1;
    const double mant_scale = static_cast<double>(1ULL << mant_bits);
    return std::ldexp(1.0 + (mant_scale - 1.0) / mant_scale, max_biased_exp - bias);
}

inline double with_subnormals_min(int exp_bits, int mant_bits) {
    const int bias = (1 << (exp_bits - 1)) - 1;
    return std::ldexp(1.0, 1 - bias - mant_bits);
}

// ─────────────────────────────────────────────────────────────────────────────
// MODE 3: LINEAR (Uniform Linear Distribution)
// ─────────────────────────────────────────────────────────────────────────────
//
// ┌─────────────────────────────────────────────────────────────────────────┐
// │ KEY IDEA: All N = 2^(exp_bits + mant_bits) bit patterns are used as     │
// │ indices into a linearly-spaced grid [min, max].                         │
// │ Every adjacent pair of representable values has the SAME absolute gap.  │
// └─────────────────────────────────────────────────────────────────────────┘
//
// REPRESENTABLE VALUES:
//   values[i] = min + i × step    for i ∈ {0, 1, ..., N-1}
//   step = (max - min) / (N - 1)
//
// Min/max reuse ALL_NORMAL's range for consistency and to avoid infinities.
//
// TRADEOFFS:
//   ✓ Intuitive, constant absolute precision
//   ✗ Poor relative precision for small values (step is large relative to value)
//   ✗ Wasted resolution — most of the bits near zero are wasted for large values
//
inline double quantize_linear(double value, int sign_bits, int exp_bits, int mant_bits) {
    if (exp_bits < 1 || mant_bits < 0 || sign_bits < 0 || sign_bits > 1) {
        throw std::invalid_argument("Invalid custom float format");
    }

    const double min_value = all_normal_min(exp_bits, mant_bits);
    const double max_value = all_normal_max(exp_bits, mant_bits);
    const double N = static_cast<double>(1ULL << (exp_bits + mant_bits));  // total levels
    const double step = (max_value - min_value) / (N - 1.0);

    if (std::isnan(value) || std::isinf(value)) {
        return (value < 0.0 && sign_bits) ? -max_value : max_value;
    }

    double sign = 1.0;
    if (value < 0.0) {
        if (sign_bits == 0) return min_value;
        sign = -1.0;
        value = -value;
    }

    if (value <= min_value) return sign * min_value;
    if (value >= max_value) return sign * max_value;

    // Truncate toward zero: floor to nearest lower grid point
    double index = std::floor((value - min_value) / step);
    return sign * (min_value + index * step);
}

// ─────────────────────────────────────────────────────────────────────────────
// MODE 4: LOGARITHMIC (Fully Log-Uniform / Geometric Distribution)
// ─────────────────────────────────────────────────────────────────────────────
//
// ┌─────────────────────────────────────────────────────────────────────────┐
// │ KEY IDEA: All N = 2^(exp_bits + mant_bits) bit patterns map to values   │
// │ that are GEOMETRICALLY spaced: values[i] = min × r^i                   │
// │ where r = (max/min)^(1/(N-1)).                                          │
// │ Every adjacent pair has the SAME relative ratio — perfectly uniform     │
// │ relative precision across the entire range.                             │
// └─────────────────────────────────────────────────────────────────────────┘
//
// Compare with ALL_NORMAL:
//   ALL_NORMAL changes step at every power-of-two boundary (octave), but
//   within each octave the spacing is linear (the mantissa bits are linear).
//   So relative precision dips slightly at the low end of each octave and
//   rises at the top — it's log-uniform per octave, not globally.
//
//   LOGARITHMIC mode places every single step as a fixed ratio, giving
//   PERFECTLY uniform relative precision:
//     step_ratio = r - 1 = (max/min)^(1/(N-1)) - 1   (constant everywhere)
//
// Quantization (truncation toward zero):
//   index = floor(log(value/min) / log(r))
//   quantized = min × r^index
//
inline double quantize_logarithmic(double value, int sign_bits, int exp_bits, int mant_bits) {
    if (exp_bits < 1 || mant_bits < 0 || sign_bits < 0 || sign_bits > 1) {
        throw std::invalid_argument("Invalid custom float format");
    }

    const double min_value = all_normal_min(exp_bits, mant_bits);
    const double max_value = all_normal_max(exp_bits, mant_bits);
    const double N = static_cast<double>(1ULL << (exp_bits + mant_bits));  // total levels
    // log_r = log(max/min) / (N-1) — work in log-space for numerical stability
    const double log_r = std::log(max_value / min_value) / (N - 1.0);

    if (std::isnan(value) || std::isinf(value)) {
        return (value < 0.0 && sign_bits) ? -max_value : max_value;
    }

    double sign = 1.0;
    if (value < 0.0) {
        if (sign_bits == 0) return min_value;
        sign = -1.0;
        value = -value;
    }

    if (value <= min_value) return sign * min_value;
    if (value >= max_value) return sign * max_value;

    // Truncate toward zero: floor to nearest lower grid point in log-space
    double index = std::floor(std::log(value / min_value) / log_r);
    return sign * min_value * std::exp(index * log_r);
}

// ─────────────────────────────────────────────────────────────────────────────
// DISPATCH: quantize_custom_float selects mode at runtime
// ─────────────────────────────────────────────────────────────────────────────

inline double quantize_custom_float(double value, int sign_bits, int exp_bits, int mant_bits,
                                    QuantizationMode mode = QuantizationMode::WITH_SUBNORMALS) {
    switch (mode) {
        case QuantizationMode::ALL_NORMAL:
            return quantize_all_normal(value, sign_bits, exp_bits, mant_bits);
        case QuantizationMode::WITH_SUBNORMALS:
            return quantize_with_subnormals(value, sign_bits, exp_bits, mant_bits);
        case QuantizationMode::LINEAR:
            return quantize_linear(value, sign_bits, exp_bits, mant_bits);
        case QuantizationMode::LOGARITHMIC:
            return quantize_logarithmic(value, sign_bits, exp_bits, mant_bits);
    }
    return quantize_with_subnormals(value, sign_bits, exp_bits, mant_bits);
}

inline double custom_float_max(int exp_bits, int mant_bits,
                               QuantizationMode mode = QuantizationMode::WITH_SUBNORMALS) {
    switch (mode) {
        case QuantizationMode::ALL_NORMAL:    return all_normal_max(exp_bits, mant_bits);
        case QuantizationMode::WITH_SUBNORMALS: return with_subnormals_max(exp_bits, mant_bits);
        case QuantizationMode::LINEAR:        return all_normal_max(exp_bits, mant_bits);
        case QuantizationMode::LOGARITHMIC:   return all_normal_max(exp_bits, mant_bits);
    }
    return with_subnormals_max(exp_bits, mant_bits);
}

inline double custom_float_min(int exp_bits, int mant_bits,
                               QuantizationMode mode = QuantizationMode::WITH_SUBNORMALS) {
    switch (mode) {
        case QuantizationMode::ALL_NORMAL:    return all_normal_min(exp_bits, mant_bits);
        case QuantizationMode::WITH_SUBNORMALS: return with_subnormals_min(exp_bits, mant_bits);
        case QuantizationMode::LINEAR:        return all_normal_min(exp_bits, mant_bits);
        case QuantizationMode::LOGARITHMIC:   return all_normal_min(exp_bits, mant_bits);
    }
    return with_subnormals_min(exp_bits, mant_bits);
}
