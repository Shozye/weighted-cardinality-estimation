#pragma once
#include <cstdint>

enum class RngEngine : std::uint8_t {
    PCG64        = 0,
    MT19937      = 1,
    XOSHIRO128PP = 2,
    XOSHIRO256PP = 3,
};

inline constexpr RngEngine kDefaultRngEngine = RngEngine::XOSHIRO128PP;
