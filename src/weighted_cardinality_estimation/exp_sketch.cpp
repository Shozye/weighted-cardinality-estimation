#include "exp_sketch.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include "hash_util.hpp"

ExpSketch::ExpSketch(std::size_t sketch_size, const std::vector<std::uint32_t>& seeds)
    : size(sketch_size), seeds_(seeds),
      M_(sketch_size, std::numeric_limits<double>::infinity())
{
    if (sketch_size == 0) { throw std::invalid_argument("Sketch size 'm' must be positive."); }
    if ((!seeds.empty() && seeds.size() != size)) { 
        throw std::invalid_argument("Seeds must have length m or 0"); 
    }
}

void ExpSketch::add(const std::string& elem, double weight)
{ 
    for (std::size_t i = 0; i < size; ++i) {
        std::uint64_t h = murmur64(elem, seeds_[i], hash_answer);
        double u = to_unit_interval(h);   
        double g = -std::log(u) / weight;
        M_[i] = std::min(g, M_[i]);
    }
} 

size_t ExpSketch::memory_usage_total() const {
    size_t total_size = 0;
    total_size += sizeof(size);
    total_size += seeds_.bytes();
    total_size += M_.capacity() * sizeof(double);
    return total_size;
}

size_t ExpSketch::memory_usage_write() const {
    return M_.capacity() * sizeof(double);
}

size_t ExpSketch::memory_usage_estimate() const {
    size_t estimate_size = M_.capacity() * sizeof(double);
    return estimate_size;
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
    double total = 0.0;
    for (double value : M_) { total += value; }
    return ((double)this->size - 1.0) / total;
}

double ExpSketch::jaccard_struct(const ExpSketch& other) const
{
    if (other.size != size) { return 0.0; }
    std::size_t equal = 0;
    for (std::size_t i = 0; i < size; ++i) {
        if (M_[i] == other.M_[i]) { ++equal; }
    }
    return static_cast<double>(equal) / static_cast<double>(size);
}

ExpSketch::ExpSketch(std::size_t sketch_size, const std::vector<std::uint32_t>& seeds, const std::vector<double>& registers)
    : size(sketch_size), seeds_(seeds), M_(registers)
{
    if (sketch_size == 0) { throw std::invalid_argument("Sketch size 'm' must be positive."); }
    if ((!seeds.empty() && seeds.size() != size)) { 
        throw std::invalid_argument("Seeds must have length m or 0"); 
    }
}

std::size_t ExpSketch::get_sketch_size() const {
    return size;
}

std::vector<std::uint32_t> ExpSketch::get_seeds() const {
    return seeds_.toVector();
}

const std::vector<double>& ExpSketch::get_registers() const {
    return M_;
}
