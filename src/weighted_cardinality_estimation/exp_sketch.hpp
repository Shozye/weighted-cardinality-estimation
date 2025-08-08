#pragma once
#include <vector>
#include <string>
#include <cstdint>

class ExpSketch {
public:
    ExpSketch(std::size_t m, const std::vector<std::uint32_t>& seeds);
    void add(const std::string& x, double weight = 1.0);
    [[nodiscard]] double estimate() const;
    [[nodiscard]] double jaccard_struct(const ExpSketch& other) const;

private:
    std::size_t m_;
    std::vector<std::uint32_t> seeds_;
    std::vector<double> M_;
};
