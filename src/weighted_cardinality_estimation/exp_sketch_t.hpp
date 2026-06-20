#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <string>
#include <cstdint>
#include "sketch.hpp"
#include "hash_util.hpp"

template <typename T>
class ExpSketchT : public Sketch, public MergeableMixin, public JaccardMixin {
public:
    ExpSketchT(std::size_t sketch_size, std::uint64_t master_seed)
        : Sketch(sketch_size, master_seed), M_(sketch_size, std::numeric_limits<T>::infinity()) {}

    ExpSketchT(std::size_t sketch_size, std::uint64_t master_seed, const std::vector<T>& registers)
        : Sketch(sketch_size, master_seed), M_(registers) {}

    void add(const std::string& elem, double weight = 1.0) override {
        validate_weight(weight);
        for (std::size_t i = 0; i < size; ++i) {
            std::uint64_t h = murmur64(elem, seeds_[i]);
            double u = to_unit_interval(h);
            double g = -std::log(u) / weight;
            M_[i] = std::min(static_cast<T>(g), M_[i]);
        }
    }

    [[nodiscard]] double estimate() const override {
        double total = 0.0;
        for (double value : M_) { total += value; }
        return (static_cast<double>(this->size) - 1.0) / total;
    }

    [[nodiscard]] double jaccard_struct(const ExpSketchT& other) const {
        if (other.size != size) { return 0.0; }
        std::size_t equal = 0;
        for (std::size_t i = 0; i < size; ++i) {
            if (M_[i] == other.M_[i]) { ++equal; }
        }
        return static_cast<double>(equal) / static_cast<double>(size);
    }

    const std::vector<T>& get_registers() const { return M_; }

    void merge(const ExpSketchT& other) {
        if (other.size != size) { throw std::invalid_argument("Cannot merge sketches of different sizes."); }
        for (std::size_t i = 0; i < size; ++i) {
            M_[i] = std::min(M_[i], other.M_[i]);
        }
    }

    [[nodiscard]] size_t memory_usage(uint64_t flags) const override {
        uint64_t f = resolve_flags(flags);
        size_t s = 0;
        if (f & MemoryFlag::REGISTERS) s += M_.capacity() * sizeof(T);
        if (f & MemoryFlag::ALL_WRITE_NO_REGISTERS) s += sizeof(size);
        if (f & MemoryFlag::SEEDS) s += seeds_.bytes();
        return s;
    }

private:
    std::vector<T> M_;
};
