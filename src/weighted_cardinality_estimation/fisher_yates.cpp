#include "fisher_yates.hpp"
#include <numeric>
#include<cmath>

FisherYates::FisherYates(std::uint32_t sketch_size):
      permInit(static_cast<std::uint32_t>(std::ceil(std::log2(sketch_size))), sketch_size),
      permWork(static_cast<std::uint32_t>(std::ceil(std::log2(sketch_size))), sketch_size) {
    std::iota(permInit.begin(), permInit.end(), 1);
}

void FisherYates::initialize(std::uint64_t rng_seed){
    this->rng_seed = rng_seed;
    permWork = permInit;
}

std::uint32_t FisherYates::get_fisher_yates_element(uint32_t index){
    uint32_t r = rand(index, this->permInit.size());
    std::uint32_t swap = permWork[index];
    permWork[index] = permWork[r];
    permWork[r] = swap;
    std::uint32_t j = permWork[index] - 1;
    return j;
}

uint32_t FisherYates::rand(uint32_t min, uint32_t max){
    this->rng_seed = this->rng_seed * 1103515245 + 12345;
    auto temp = (unsigned)(this->rng_seed/65536) % 32768;
    return (temp % (max-min)) + min;
}

std::uint32_t FisherYates::bytes_total() const {
    size_t total_size = 0;
    total_size += sizeof(rng_seed); // 8
    total_size += permInit.bytes(); // m * ceil(log_2 m)
    total_size += permWork.bytes(); // m * ceil(log_2 m)
    return total_size; // 2m ceil(log_2 m) + 8
}

std::uint32_t FisherYates::bytes_write() const {
    size_t total_size = 0;
    total_size += sizeof(rng_seed); // 8
    total_size += permWork.bytes(); // m * ceil(log_2 m)
    return total_size; // m ceil(log_2 m) + 8
}