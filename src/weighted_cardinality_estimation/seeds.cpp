#include "seeds.hpp"

Seeds::Seeds(std::uint64_t master_seed)
    : master_seed_(master_seed) {}

std::uint32_t Seeds::get(uint32_t index) const {
    return static_cast<std::uint32_t>(splitmix64(splitmix64(master_seed_) ^ splitmix64(index)));
}

std::uint32_t Seeds::bytes() const {
    return sizeof(std::uint64_t);
}

std::uint32_t Seeds::operator[](uint32_t index) const {
    return get(index);
}
