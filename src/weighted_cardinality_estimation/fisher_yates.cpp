#include "fisher_yates.hpp"
#include <numeric>

FisherYates::FisherYates(std::uint32_t sketch_size):
      permInit(sketch_size),
      permWork(sketch_size) {
    std::iota(permInit.begin(), permInit.end(), 1);
}

void FisherYates::initialize(std::uint64_t rng_seed){
    this->rng_seed = rng_seed;
    permWork = permInit;
}

std::uint32_t FisherYates::get_fisher_yates_element(uint32_t index){
    uint32_t r = rand(index, this->permInit.size());
    std::swap(permWork[index], permWork[r]);
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
    total_size += sizeof(rng_seed);
    total_size += permInit.capacity() * sizeof(uint32_t);
    total_size += permWork.capacity() * sizeof(uint32_t);
    return total_size;
}

std::uint32_t FisherYates::bytes_write() const {
    size_t total_size = 0;
    total_size += sizeof(rng_seed);
    total_size += permWork.capacity() * sizeof(uint32_t);
    return total_size;
}