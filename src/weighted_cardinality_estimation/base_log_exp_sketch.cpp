#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "hash_util.hpp"
#include<cstring>
#include "utils.hpp"
#include"base_log_exp_sketch.hpp"

BaseLogExpSketch::BaseLogExpSketch(
    std::size_t sketch_size, 
    const std::vector<std::uint32_t>& seeds, 
    uint8_t amount_bits,
    float logarithm_base
)
    : Sketch(sketch_size, seeds),
      amount_bits_(amount_bits),
      logarithm_base(logarithm_base),
      r_max((1 << (amount_bits - 1)) - 1),
      r_min(-(1 << (amount_bits - 1)) + 1),
      M_(amount_bits, sketch_size)
{
    if (amount_bits == 0) { throw std::invalid_argument("Amount of bits 'b' must be positive."); }
    std::fill(M_.begin(), M_.end(), r_min);
}

BaseLogExpSketch::BaseLogExpSketch(
    std::size_t sketch_size, 
    const std::vector<std::uint32_t>& seeds, 
    std::uint8_t amount_bits, 
    float logarithm_base,
    const std::vector<int>& registers
)
    : Sketch(sketch_size, seeds),
      amount_bits_(amount_bits),
      logarithm_base(logarithm_base),
      r_max((1 << (amount_bits - 1)) - 1),
      r_min(-(1 << (amount_bits - 1)) + 1),
      M_(amount_bits, sketch_size)
{
    if (amount_bits == 0) { throw std::invalid_argument("Amount of bits 'b' must be positive."); }
    if (registers.size() != sketch_size) { throw std::invalid_argument("Invalid state: registers vector size mismatch"); }
    for (std::size_t i = 0; i < size; ++i) {
        M_[i] = registers[i];
    }
}


size_t BaseLogExpSketch::memory_usage_total() const {
    size_t total_size = 0;
    total_size += sizeof(size); // 8
    total_size += seeds_.bytes(); // m * 4
    total_size += M_.bytes(); // mb/8
    total_size += sizeof(amount_bits_); // 1
    total_size += sizeof(r_max); // 4
    total_size += sizeof(r_min); // 4
    total_size += sizeof(logarithm_base); // 4
    return total_size; // m * 4 + mb/8 + 21
}

size_t BaseLogExpSketch::memory_usage_write() const {
    size_t write_size = 0;
    write_size += M_.bytes(); // mb/8
    return write_size; // mb/8
}

size_t BaseLogExpSketch::memory_usage_estimate() const {
    size_t estimate_size = M_.bytes(); // mb/8
    estimate_size += sizeof(logarithm_base); // 4
    return estimate_size; // mb/8 + 4
}

std::uint8_t BaseLogExpSketch::get_amount_bits() const { return amount_bits_; }
float BaseLogExpSketch::get_logarithm_base() const { return logarithm_base; }
std::vector<int> BaseLogExpSketch::get_registers() const {
    return std::vector<int>(M_.begin(), M_.end());
}

void BaseLogExpSketch::add(const std::string& elem, double weight){ 
    for (std::size_t i = 0; i < size; ++i) {
        std::uint64_t h = murmur64(elem, seeds_[i]);
        double u = to_unit_interval(h);   
        double g = -std::log(u) / weight;
        int q = static_cast<int>(std::floor(-std::log(g)/std::log(logarithm_base)));
        q = std::min(q, r_max);
        if (q > M_[i]){
            M_[i] = q;
        }
    }
} 


double BaseLogExpSketch::derivative_of_log_probability(int r, double lambda) const {
    const double k = this->logarithm_base;

    // // r = r_min
    // if (r == this->r_min) {
    //     // (ln P[M = r_min])' = k^{r_min+1} e^{-λ k^{r_min+1}} / (1 - e^{-λ k^{r_min+1}})
    //     const double k_rmin_plus_1 = std::pow(k, static_cast<double>(this->r_min + 1));
    //     const double exp_term      = std::exp(-lambda * k_rmin_plus_1);
    //     const double denom         = 1.0 - exp_term;
    //     // You may want to guard against denom ~ 0 in production code
    //     return k_rmin_plus_1 * exp_term / denom;
    // }

    // // r = r_max
    // if (r == this->r_max) {
    //     // (ln P[M = r_max])' = -k^{r_max}
    //     const double k_rmax = std::pow(k, static_cast<double>(this->r_max));
    //     return -k_rmax;
    // }

    // r_min < r < r_max
    // (ln P[M = r])' = -k^r * (1 - k e^{-λ k^r (k-1)}) / (1 - e^{-λ k^r (k-1)})

    const double k_r      = std::pow(k, static_cast<double>(-r-1));
    const double front_term = -k_r;
    const double a        = - lambda * k_r * (k-1);              // k^r (k-1)
    const double exp_term = std::exp(a);        // e^{-λ k^r (k-1)}
    const double denom    = 1.0 - exp_term;
    const double numer    = 1.0 - (k * exp_term);

    return front_term * (numer / denom);
}

double BaseLogExpSketch::second_derivative_of_log_probability(int r, double lambda) const {
    const double k = this->logarithm_base;

    // // r = r_min
    // if (r == this->r_min) {
    //     // (ln P[M = r_min])'' =
    //     // - k^{2 r_min + 2} * e^{-λ k^{r_min+1}} / (1 - e^{-λ k^{r_min+1}})^2
    //     const double k_rmin_plus_1 = std::pow(k, static_cast<double>(this->r_min + 1));
    //     const double exp_term      = std::exp(-lambda * k_rmin_plus_1);
    //     const double denom         = 1.0 - exp_term;
    //     const double k_pow         = k_rmin_plus_1 * k_rmin_plus_1; // k^{2 r_min + 2}
    //     // Again, denom could be very small if λ is tiny
    //     return -k_pow * exp_term / (denom * denom);
    // }

    // // r = r_max
    // if (r == this->r_max) {
    //     // (ln P[M = r_max])'' = 0
    //     return 0.0;
    // }

    // r_min < r < r_max
    // (ln P[M = r])'' =
    // - (k-1)^2 k^{2r} * e^{-λ k^r (k-1)} / (1 - e^{-λ k^r (k-1)})^2


    const double k_r      = std::pow(k, static_cast<double>(-r-1));
    const double a        = - lambda * k_r * (k-1);              // k^r (k-1)
    const double exp_term = std::exp(a);        // e^{-λ k^r (k-1)}
    const double numerator = exp_term;
    const double denom    = (1.0 - exp_term)*(1.0 - exp_term);

    const double front_term = - (k-1)*(k-1)* k_r * k_r;

    return front_term * numerator/denom;
}


/*
Here we could use this formula:
\begin{equation}
    \overline\lambda = \frac{(k-1)m }{\ln k\sum_{i=1}^m k^{-R_i}}
\end{equation}
now we use formula 
\begin{equation}
    \overline\lambda = \frac{(m-1}{\sum_{i=1}^m k^{-R_i}}
\end{equation}
*/
double BaseLogExpSketch::estimate_fast() const {
    double tmp_sum = 0.0;
    for(int r: M_) { 
        tmp_sum += std::pow(logarithm_base, -r);
    }
    return (double)(this->size-1) / tmp_sum;
}


/*
Here goal is to compute
\frac{\ell'(\lambda^{(t)})}{\ell''(\lambda^{(t)})},
\ell'(\lambda)  = \sum_{i=1}^{m} \left( \ln  \Pr[M=r_i] \right)'
\ell''(\lambda) = \sum_{i=1}^{m} \left( \ln  \Pr[M=r_i] \right)''
*/
double BaseLogExpSketch::ffunc_divided_by_dffunc(double lambda) const {
    double loss_derivative = 0;
    double loss_second_derivative = 0;
    for (int r: M_) {
        loss_derivative += derivative_of_log_probability(r, lambda);
        loss_second_derivative += second_derivative_of_log_probability(r, lambda); 
    }
    return loss_derivative / loss_second_derivative;
}

double BaseLogExpSketch::newton(double c0) const {
    double c1 = c0 - ffunc_divided_by_dffunc(c0);
    // std::cout << "c1=" << c1 << " c0=" << c0 << '\n';
    int it = 0;
    while (std::abs(c1 - c0) > NEWTON_MAX_ERROR) {
        c0 = c1;
        c1 = c0 - ffunc_divided_by_dffunc(c0);
        // std::cout << "c1=" << c1 << " c0=" << c0 << '\n';
        it += 1;
        if (it > NEWTON_MAX_ITERATIONS){ break; }
    }
    return c1;
}

double BaseLogExpSketch::estimate() const {
    return newton(estimate_fast());
}
