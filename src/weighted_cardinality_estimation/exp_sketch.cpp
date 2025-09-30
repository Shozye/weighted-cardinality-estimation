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
    std::uint64_t hash_answer[2];
    for (std::size_t i = 0; i < m_; ++i) {
        std::uint64_t h = murmur64(x, seeds_[i], hash_answer);
        double u = to_unit_interval(h);   
        double g = -std::log(u) / weight;
        if (g < M_[i]) M_[i] = g;
    }
} 

size_t ExpSketch::memory_usage_total() const {
    size_t total_size = 0;
    total_size += sizeof(m_);
    total_size += seeds_.capacity() * sizeof(uint32_t);
    total_size += M_.capacity() * sizeof(double);
    return total_size;
}

size_t ExpSketch::memory_usage_write() const {
    return M_.capacity() * sizeof(double);
}

void ExpSketch::add_many(const std::vector<std::string>& elems,
                                  const std::vector<double>& weights) {
    if (elems.size() != weights.size()){
        throw std::invalid_argument("add_many: elems and weights size mismatch");
    }
    for (std::size_t i = 0; i < elems.size(); ++i) {
        this->add(elems[i], weights[i]);
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

ExpSketch::ExpSketch(std::size_t m, const std::vector<std::uint32_t>& seeds, const std::vector<double>& registers)
    : m_(m), seeds_(seeds), M_(registers)
{
}


std::size_t ExpSketch::get_m() const {
    return m_;
}

const std::vector<std::uint32_t>& ExpSketch::get_seeds() const {
    return seeds_;
}

const std::vector<double>& ExpSketch::get_registers() const {
    return M_;
}
