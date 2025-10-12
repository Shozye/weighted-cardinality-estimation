#pragma once
#include <vector>
#include <cstdint>

#define NEWTON_MAX_ITERATIONS 5
#define NEWTON_MAX_ERROR 1e-5

void print_vector(std::vector<int> vec);
std::vector<uint32_t> range(uint32_t min, uint32_t max);
std::uint32_t argmax(std::vector<double> vec);
