#pragma once
#include <vector>
#include <string>
#include <cstdint>

class FastQSketch {
public:
    FastQSketch(std::size_t m, const std::vector<std::uint32_t>& seeds, std::uint8_t amount_bits);
    void add(const std::string& x, double weight = 1.0);
    [[nodiscard]] double estimate();

private:
    int rand(int min, int max);

    double initialValue();
    double ffunc(double w);
    double dffunc(double w);
    double Newton(double c0);

    void update_treshold();

    std::size_t m_;
    std::vector<std::uint32_t> seeds_; 
    std::int32_t r_max;
    std::int32_t r_min; 

    std::vector<int> M_; // sketch structure with elements between < r_min ... r_max >

    std::vector<uint32_t> permInit; // static structure, only used to fast copy to permWork
    std::vector<uint32_t> permWork; // used at the beginning of every update
    
    uint64_t rng_seed; 
    int min_sketch_value; 
    double min_value_to_change_sketch; // that's 2**{-min_sketch_value}
};
