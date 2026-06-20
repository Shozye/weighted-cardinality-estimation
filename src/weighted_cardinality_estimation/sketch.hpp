#pragma once

#include "seeds.hpp"
#include "memory_flag.hpp"
#include <cstddef>
#include <cmath>
#include <stdexcept>

// ─── Capability markers ───────────────────────────────────────────────────────

class CardinalitySketch {
public:
    virtual ~CardinalitySketch() = default;
};

class WeightedMixin {
public:
    virtual ~WeightedMixin() = default;
};

class MergeableMixin {
public:
    virtual ~MergeableMixin() = default;
};

class JaccardMixin {
public:
    virtual ~JaccardMixin() = default;
};

class NewtonMixin {
public:
    virtual ~NewtonMixin() = default;
};

// ─── Shared state base ────────────────────────────────────────────────────────
// Holds size + seeds; not exposed to Python.

class SketchBase : public CardinalitySketch {
public:
    SketchBase(std::size_t sketch_size, std::uint64_t master_seed)
        : size(sketch_size), seeds_(master_seed)
    {
        if (sketch_size == 0) throw std::invalid_argument("Sketch size 'm' must be positive.");
    }
    virtual ~SketchBase() = default;

    std::size_t get_sketch_size() const { return size; }
    std::uint64_t get_master_seed() const { return seeds_.get_master_seed(); }

    virtual double estimate() const = 0;
    [[nodiscard]] virtual size_t memory_usage(uint64_t flags) const = 0;

    virtual void add(const std::string& elem) = 0;

    void add_many(const std::vector<std::string>& elems) {
        for (const auto& e : elems) this->add(e);
    }

protected:
    std::size_t size;
    Seeds seeds_;

    static void validate_weight(double weight) {
        if (weight <= 0.0 || std::isnan(weight) || std::isinf(weight))
            throw std::invalid_argument("Weight must be a finite positive number.");
    }

    // Resolves the flag bitmask: if TOTAL bit is set, the other bits are EXCLUSIONS;
    // otherwise the bits specify which components to INCLUDE.
    static uint64_t resolve_flags(uint64_t flags) {
        constexpr uint64_t ALL_COMPONENTS =
            MemoryFlag::REGISTERS | MemoryFlag::ALL_WRITE_NO_REGISTERS |
            MemoryFlag::FISHER_YATES_PERM_INIT | MemoryFlag::FISHER_YATES_NON_PERM_INIT |
            MemoryFlag::SEEDS;
        if (flags & MemoryFlag::TOTAL) {
            return ALL_COMPONENTS & ~(flags & ~MemoryFlag::TOTAL);
        }
        return flags;
    }
};

// ─── UnweightedSketch ─────────────────────────────────────────────────────────
// Base for sketches that treat all elements as having implicit weight 1.

class UnweightedSketch : public SketchBase {
public:
    UnweightedSketch(std::size_t sketch_size, std::uint64_t master_seed)
        : SketchBase(sketch_size, master_seed) {}

    virtual void add(const std::string& elem) = 0;

    void add_many(const std::vector<std::string>& elems) {
        for (const auto& e : elems) this->add(e);
    }
};

// ─── Sketch (weighted) ────────────────────────────────────────────────────────
// Base for all weighted cardinality sketches.
// Inherits WeightedMixin so Python can use isinstance(sketch, WeightedMixin).

class Sketch : public SketchBase, public WeightedMixin {
public:
    Sketch(std::size_t sketch_size, std::uint64_t master_seed)
        : SketchBase(sketch_size, master_seed) {}

    virtual void add(const std::string& elem, double weight) = 0;

    // Satisfy SketchBase: unweighted add dispatches with weight=1
    void add(const std::string& elem) override { this->add(elem, 1.0); }

    void add_many(const std::vector<std::string>& elems,
                  const std::vector<double>& weights) {
        if (elems.size() != weights.size())
            throw std::invalid_argument("add_many: elems and weights size mismatch");
        for (std::size_t i = 0; i < elems.size(); ++i)
            this->add(elems[i], weights[i]);
    }

    // Convenience: add_many with implicit weight=1
    void add_many(const std::vector<std::string>& elems) {
        for (const auto& e : elems) this->add(e, 1.0);
    }
};
