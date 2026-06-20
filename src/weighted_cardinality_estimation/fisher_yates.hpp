#pragma once
#include <compact_vector.hpp>
#include <cstdint>
#include <random>
#include <variant>
#include "memory_flag.hpp"
#include "pcg_random.hpp"
#include "rng_engine_type.hpp"
#include "xoshiro.hpp"

class FisherYates {
public:
    explicit FisherYates(std::uint32_t sketch_size, RngEngine engine = kDefaultRngEngine);
    void initialize(const std::string& elem);
    uint32_t get_fisher_yates_element(uint32_t index);
    [[nodiscard]] size_t memory_usage(uint64_t flags) const;
    [[nodiscard]] RngEngine engine_type() const { return engine_type_; }
private:
    std::variant<pcg64, std::mt19937_64, xoshiro128pp, xoshiro256pp> rng_engine_;
    RngEngine engine_type_;
    compact::vector<uint32_t> permInit;
    compact::vector<uint32_t> permWork;
};
