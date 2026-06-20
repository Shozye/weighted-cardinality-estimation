#include "fisher_yates.hpp"
#include "hash_util.hpp"
#include <cmath>
#include <numeric>

static std::variant<pcg64, std::mt19937_64, xoshiro128pp, xoshiro256pp> make_engine(RngEngine engine) {
    switch (engine) {
        case RngEngine::MT19937:      return std::mt19937_64{std::random_device{}()};
        case RngEngine::XOSHIRO128PP: return xoshiro128pp{std::random_device{}()};
        case RngEngine::XOSHIRO256PP: return xoshiro256pp{std::random_device{}()};
        default:                      return pcg64{std::random_device{}()};
    }
}

FisherYates::FisherYates(std::uint32_t sketch_size, RngEngine engine)
    : rng_engine_(make_engine(engine)),
      engine_type_(engine),
      permInit(std::max(1U, static_cast<std::uint32_t>(std::ceil(std::log2(sketch_size)))), sketch_size),
      permWork(std::max(1U, static_cast<std::uint32_t>(std::ceil(std::log2(sketch_size)))), sketch_size) {
    std::iota(permInit.begin(), permInit.end(), 0);
}

void FisherYates::initialize(const std::string& elem) {
    std::uint64_t rng_seed = murmur64(elem, 1);
    std::visit([rng_seed](auto& rng) { rng.seed(rng_seed); }, rng_engine_);
    permWork = permInit;
}

std::uint32_t FisherYates::get_fisher_yates_element(uint32_t index) {
    std::uniform_int_distribution<uint32_t> dist(index, permInit.size() - 1);
    uint32_t r = std::visit([&dist](auto& rng) { return dist(rng); }, rng_engine_);

    std::uint32_t tmp = permWork[index];
    permWork[index] = permWork[r];
    permWork[r] = tmp;

    return permWork[index];
}

size_t FisherYates::memory_usage(uint64_t flags) const {
    size_t s = 0;
    if (flags & MemoryFlag::FISHER_YATES_PERM_INIT) s += permInit.bytes();
    if (flags & MemoryFlag::FISHER_YATES_NON_PERM_INIT)
        s += std::visit([](const auto& rng) { return sizeof(rng); }, rng_engine_) + permWork.bytes();
    return s;
}
