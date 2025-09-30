#pragma once
#include <vector>
#include <string>
#include <cstdint>

class ExpSketch {
public:
    ExpSketch(std::size_t m, const std::vector<std::uint32_t>& seeds);
    void add(const std::string& x, double weight = 1.0);
    void add_many(const std::vector<std::string>& elems, const std::vector<double>& weights);
    [[nodiscard]] double estimate() const;
    [[nodiscard]] double jaccard_struct(const ExpSketch& other) const;

    ExpSketch(std::size_t m, const std::vector<std::uint32_t>& seeds, const std::vector<double>& registers);

    std::size_t get_m() const;
    const std::vector<std::uint32_t>& get_seeds() const;
    const std::vector<double>& get_registers() const;
private:
    std::size_t m_;
    std::vector<std::uint32_t> seeds_;
    std::vector<double> M_;
};
