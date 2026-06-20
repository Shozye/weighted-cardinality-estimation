#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <string>
#include <cstdint>
#include "fisher_yates.hpp"
#include "hash_util.hpp"
#include "rng_engine_type.hpp"
#include "sketch.hpp"

template <typename T>
class FastExpSketchT : public Sketch, public MergeableMixin, public JaccardMixin {
public:
    FastExpSketchT(std::size_t sketch_size, std::uint64_t master_seed, RngEngine engine = kDefaultRngEngine)
        : Sketch(sketch_size, master_seed),
          M_(sketch_size, std::numeric_limits<T>::infinity()),
          fisher_yates(sketch_size, engine),
          max(std::numeric_limits<T>::infinity())
    {}

    FastExpSketchT(std::size_t sketch_size, std::uint64_t master_seed, const std::vector<T>& registers, RngEngine engine = kDefaultRngEngine)
        : Sketch(sketch_size, master_seed),
          M_(registers),
          fisher_yates(sketch_size, engine)
    {
        max = *std::max_element(M_.begin(), M_.end());
    }

    void add(const std::string& elem, double weight = 1.0) override {
        validate_weight(weight);
        double S = 0;
        bool updateMax = false;

        fisher_yates.initialize(elem);
        for (std::size_t k = 0; k < size; ++k) {
            std::uint64_t hashed = murmur64(elem, seeds_[k]);
            double U = to_unit_interval(hashed);
            double E = -std::log(U) / weight;

            S += E / static_cast<double>(size - k);
            if (S >= static_cast<double>(max)) { break; }

            std::uint32_t j = fisher_yates.get_fisher_yates_element(k);

            if (M_[j] == max) { updateMax = true; }
            M_[j] = std::min(static_cast<T>(S), M_[j]);
        }

        if (updateMax) {
            max = *std::max_element(M_.begin(), M_.end());
        }
    }

    [[nodiscard]] double estimate() const override {
        double total = 0.0;
        for (T val : M_) { total += static_cast<double>(val); }
        return (static_cast<double>(size) - 1.0) / total;
    }

    [[nodiscard]] double jaccard_struct(const FastExpSketchT& other) const {
        if (other.size != size) { return 0.0; }
        std::size_t equal = 0;
        for (std::size_t i = 0; i < size; ++i) {
            if (M_[i] == other.M_[i]) { ++equal; }
        }
        return static_cast<double>(equal) / static_cast<double>(size);
    }

    const std::vector<T>& get_registers() const { return M_; }

    void merge(const FastExpSketchT& other) {
        if (other.size != size) { throw std::invalid_argument("Cannot merge sketches of different sizes."); }
        for (std::size_t i = 0; i < size; ++i) {
            M_[i] = std::min(M_[i], other.M_[i]);
        }
        max = *std::max_element(M_.begin(), M_.end());
    }

    [[nodiscard]] size_t memory_usage(uint64_t flags) const override {
        uint64_t f = resolve_flags(flags);
        size_t s = 0;
        if (f & MemoryFlag::REGISTERS) s += M_.capacity() * sizeof(T);
        if (f & MemoryFlag::ALL_WRITE_NO_REGISTERS) s += sizeof(size) + sizeof(max);
        s += fisher_yates.memory_usage(f);
        if (f & MemoryFlag::SEEDS) s += seeds_.bytes();
        return s;
    }

private:
    std::vector<T> M_;
    FisherYates fisher_yates;
    T max;
};
