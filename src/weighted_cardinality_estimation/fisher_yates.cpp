#include "fisher_yates.hpp"
#include "hash_util.hpp"
#include <numeric>
#include<cmath>

FisherYates::FisherYates(std::uint32_t sketch_size):
      rng_engine(std::random_device{}()),
      permInit(static_cast<std::uint32_t>(std::ceil(std::log2(sketch_size))), sketch_size),
      permWork(static_cast<std::uint32_t>(std::ceil(std::log2(sketch_size))), sketch_size) {
    std::iota(permInit.begin(), permInit.end(), 0);
}

void FisherYates::initialize(const std::string& elem){
    std::uint64_t rng_seed = murmur64(elem, 1);
    this->rng_engine.seed(rng_seed);
    permWork = permInit;
}

std::uint32_t FisherYates::get_fisher_yates_element(uint32_t index){
    std::uniform_int_distribution<uint32_t> dist(index, this->permInit.size() - 1);
    uint32_t r = dist(this->rng_engine);

    std::uint32_t swap = permWork[index];
    permWork[index] = permWork[r];
    permWork[r] = swap;

    return permWork[index];
}

std::uint32_t FisherYates::bytes_total() const {
    size_t total_size = 0;
    total_size += permInit.bytes(); // m * ceil(log_2 m)
    total_size += permWork.bytes(); // m * ceil(log_2 m)
    return total_size; // 2m ceil(log_2 m) + 8
}

std::uint32_t FisherYates::bytes_write() const {
    size_t total_size = 0;
    total_size += permWork.bytes(); // m * ceil(log_2 m)
    return total_size; // m ceil(log_2 m) + 8
}