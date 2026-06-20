#pragma once
#include <algorithm>
#include <limits>
#include <string>
#include <vector>
#include <cstdint>
#include "sketch.hpp"
#include "hash_util.hpp"

// Unweighted MinHash (k-mins sketch, Broder CPM 2000).
// Each register stores min_x{ -log(U(x,i)) } over all inserted elements.
// Estimator: (m-1) / sum(M[i])
// Supports merge (element-wise min) and Jaccard (fraction of equal registers).
class MinHash : public UnweightedSketch, public MergeableMixin, public JaccardMixin {
public:
    MinHash(std::size_t sketch_size, std::uint64_t master_seed)
        : UnweightedSketch(sketch_size, master_seed), M_(sketch_size, std::numeric_limits<double>::infinity()) {}

    MinHash(std::size_t sketch_size, std::uint64_t master_seed,
            const std::vector<double>& registers)
        : UnweightedSketch(sketch_size, master_seed), M_(registers) {}

    void add(const std::string& elem) override {
        for (std::size_t i = 0; i < size; ++i) {
            double g = -std::log(to_unit_interval(murmur64(elem, seeds_[i])));
            if (g < M_[i]) M_[i] = g;
        }
    }

    [[nodiscard]] double estimate() const override {
        double sum = 0.0;
        for (double v : M_) sum += v;
        return static_cast<double>(size - 1) / sum;
    }

    [[nodiscard]] double jaccard_struct(const MinHash& other) const {
        std::size_t equal = 0;
        for (std::size_t i = 0; i < size; ++i)
            if (M_[i] == other.M_[i]) ++equal;
        return static_cast<double>(equal) / static_cast<double>(size);
    }

    void merge(const MinHash& other) {
        if (other.size != size) throw std::invalid_argument("Cannot merge sketches of different sizes.");
        for (std::size_t i = 0; i < size; ++i)
            M_[i] = std::min(M_[i], other.M_[i]);
    }

    [[nodiscard]] const std::vector<double>& get_registers() const { return M_; }

    [[nodiscard]] size_t memory_usage(uint64_t flags) const override {
        uint64_t f = resolve_flags(flags);
        size_t s = 0;
        if (f & MemoryFlag::REGISTERS) s += M_.capacity() * sizeof(double);
        if (f & MemoryFlag::ALL_WRITE_NO_REGISTERS) s += sizeof(size);
        if (f & MemoryFlag::SEEDS) s += seeds_.bytes();
        return s;
    }

private:
    std::vector<double> M_;
};
