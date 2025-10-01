#pragma once
#include <vector>
#include <string>
#include <cstdint>

class QSketchDyn {
public:
    QSketchDyn(std::size_t m, const std::vector<std::uint32_t>& seeds, std::uint8_t amount_bits, std::uint32_t g_seed);
    void add(const std::string& elem, double weight = 1.0);
    void add_many(const std::vector<std::string>& elems, const std::vector<double>& weights);
    [[nodiscard]] double estimate() const;

    QSketchDyn(
        std::size_t m,
        std::uint8_t amount_bits,
        std::uint32_t g_seed,
        const std::vector<std::uint32_t>& seeds,
        const std::vector<int>& registers,
        const std::vector<int>& t_histogram,
        double cardinality
    );
    std::size_t get_m() const;
    std::uint8_t get_amount_bits() const;
    std::uint32_t get_g_seed() const;
    const std::vector<std::uint32_t>& get_seeds() const;
    const std::vector<int>& get_registers() const;
    const std::vector<int>& get_t_histogram() const;
    double get_cardinality() const;

    [[nodiscard]] size_t memory_usage_total() const;
    [[nodiscard]] size_t memory_usage_write() const;
    [[nodiscard]] size_t memory_usage_estimate() const;

private:
    std::size_t m_;
    std::uint8_t amount_bits_;
    std::int32_t r_min;
    std::int32_t r_max;
    std::vector<std::uint32_t> seeds_;
    std::uint32_t g_seed_;
    std::vector<int> k_idx_;

    double cardinality_;
    double q_r_;
    std::vector<int> R_;
    std::vector<int> T_;

    uint64_t hash_answer[2];
};

