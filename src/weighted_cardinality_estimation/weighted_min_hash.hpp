#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include "sketch.hpp"

// Weighted MinHash based on Cohen, Katzir & Yehezkel (IPL 2015).
// Generalises the max-sketch cardinality estimator to weighted streams via
// the Beta(w,1) transform: h(x)^(1/w) ~ Beta(w,1).
// Estimator: w_hat = m*(m-1) / sum_k(1 - M[k])  (Method-1, ARE=1.00)
class WeightedMinHash : public Sketch {
public:
    WeightedMinHash(std::size_t sketch_size, const std::vector<std::uint32_t>& seeds);
    WeightedMinHash(std::size_t sketch_size, const std::vector<std::uint32_t>& seeds, const std::vector<double>& registers);

    void add(const std::string& elem, double weight = 1.0) override;
    [[nodiscard]] double estimate() const override;
    [[nodiscard]] const std::vector<double>& get_registers() const { return M_; }

    [[nodiscard]] size_t memory_usage_total() const override;
    [[nodiscard]] size_t memory_usage_write() const override;
    [[nodiscard]] size_t memory_usage_estimate() const override;

private:
    std::vector<double> M_; // max registers, initialised to 0
};
