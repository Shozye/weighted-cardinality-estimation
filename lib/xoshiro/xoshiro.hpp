// xoshiro128++ and xoshiro256++ by David Blackman and Sebastiano Vigna (2019)
// Public domain (CC0). Reference: https://prng.di.unimi.it/
// C++ adaptation with seed() and UniformRandomBitGenerator interface.
#pragma once
#include <cstdint>
#include <limits>

// Splitmix64 for seeding (Vigna)
namespace xoshiro_detail {
inline uint64_t splitmix64(uint64_t& x) {
    uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
} // namespace xoshiro_detail

class xoshiro128pp {
    uint32_t s[4]{};
    static uint32_t rotl(uint32_t x, int k) { return (x << k) | (x >> (32 - k)); }
public:
    using result_type = uint32_t;
    xoshiro128pp() = default;
    explicit xoshiro128pp(uint64_t seed) { this->seed(seed); }
    void seed(uint64_t v) {
        uint64_t a = xoshiro_detail::splitmix64(v);
        uint64_t b = xoshiro_detail::splitmix64(v);
        s[0] = uint32_t(a); s[1] = uint32_t(a >> 32);
        s[2] = uint32_t(b); s[3] = uint32_t(b >> 32);
    }
    static constexpr result_type min() { return 0; }
    static constexpr result_type max() { return std::numeric_limits<uint32_t>::max(); }
    result_type operator()() {
        const uint32_t result = rotl(s[0] + s[3], 7) + s[0];
        const uint32_t t = s[1] << 9;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 11);
        return result;
    }
};

class xoshiro256pp {
    uint64_t s[4]{};
    static uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
public:
    using result_type = uint64_t;
    xoshiro256pp() = default;
    explicit xoshiro256pp(uint64_t seed) { this->seed(seed); }
    void seed(uint64_t v) {
        s[0] = xoshiro_detail::splitmix64(v);
        s[1] = xoshiro_detail::splitmix64(v);
        s[2] = xoshiro_detail::splitmix64(v);
        s[3] = xoshiro_detail::splitmix64(v);
    }
    static constexpr result_type min() { return 0; }
    static constexpr result_type max() { return std::numeric_limits<uint64_t>::max(); }
    result_type operator()() {
        const uint64_t result = rotl(s[0] + s[3], 23) + s[0];
        const uint64_t t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 45);
        return result;
    }
};
