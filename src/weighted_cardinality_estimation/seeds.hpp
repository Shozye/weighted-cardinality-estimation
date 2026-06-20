#pragma once
#include <compact_vector.hpp>
#include <cstdint>
#include <vector>

inline std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

class Seeds {
public:
    Seeds(std::uint64_t master_seed); // derive seeds lazily via splitmix64
    std::uint32_t get(uint32_t index) const; 
    std::uint32_t bytes() const;
    std::uint32_t operator[](uint32_t index) const;
    std::uint64_t get_master_seed() const { return master_seed_; }
private:
    std::uint64_t master_seed_;
};
