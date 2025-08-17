#pragma once
#include <vector>
#include <string>
#include <cstdint>

class FastQSketch {
public:
    FastQSketch(std::size_t sketch_size, const std::vector<std::uint32_t>& seeds, std::uint8_t amount_bits);
    void add(const std::string& elem, double weight = 1.0);
    void add_many(const std::vector<std::string>& elems,
                         const std::vector<double>& weights);
    [[nodiscard]] double estimate();

private:
    uint32_t rand(uint32_t min, uint32_t max);

    double initialValue();
    double ffunc(double);
    double dffunc(double);
    double Newton(double);

    void update_treshold();

    std::size_t sketch_size_; // amount of registers used in sketch. Sketch uses linear memory to m and increases accuracy based on m
    std::vector<std::uint32_t> seeds_; // seeds used to hash 
    std::int32_t r_max; // maximum possible value in sketch due to amount of bits per register
    std::int32_t r_min; // minimum possible value in sketch due to amount of bits per register

    std::vector<int> M_; // sketch structure with elements between < r_min ... r_max >

    std::vector<uint32_t> permInit; // static structure, only used to fast copy to permWork
    std::vector<uint32_t> permWork; // used at the beginning of every update
    
    uint64_t rng_seed; 
    int min_sketch_value; 
    double min_value_to_change_sketch; // that's 2**{-min_sketch_value}
};
