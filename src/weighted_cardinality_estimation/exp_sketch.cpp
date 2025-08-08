#include "exp_sketch.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
#include "hash_util.hpp"

ExpSketch::ExpSketch(std::size_t m, const std::vector<std::uint32_t>& seeds)
    : m_(m), seeds_(seeds),
      M_(m, std::numeric_limits<double>::infinity())
{
    if (seeds_.size() != m_)
        throw std::invalid_argument("Seeds vector must have length m");
}

void ExpSketch::add(const std::string& x, double weight)
{ 
    for (std::size_t i = 0; i < m_; ++i) {
        std::uint64_t h = murmur64(x, seeds_[i]);
        double u = to_unit_interval(h);   
        double g = -std::log(u) / weight;
        if (g < M_[i]) M_[i] = g;
    }
} 

double ExpSketch::estimate() const
{
    double sum = 0.0;
    for (double v : M_) sum += v;
    return (m_ - 1.0) / sum;
}

double ExpSketch::jaccard_struct(const ExpSketch& other) const
{
    if (other.m_ != m_) return 0.0;
    std::size_t equal = 0;
    for (std::size_t i = 0; i < m_; ++i)
        if (M_[i] == other.M_[i]) ++equal;
    return static_cast<double>(equal) / static_cast<double>(m_);
}
