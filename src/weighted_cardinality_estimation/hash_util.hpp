#pragma once
#include <cstdint>
#include <string>
#include "MurmurHash3.h"



inline std::uint64_t murmur64(const std::string& key, std::uint32_t seed)
{
    std::uint64_t hash128[2];                       // 16 B buffer
    MurmurHash3_x64_128(key.data(),
                        static_cast<int>(key.size()),
                        seed,
                        hash128);
    return hash128[0];                              // pierwsze 64 bity
}

inline double to_unit_interval(std::uint64_t h)
{
    return static_cast<double>(h + 1ULL) /
           static_cast<double>(std::numeric_limits<std::uint64_t>::max());
}
