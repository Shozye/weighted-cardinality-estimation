#pragma once
#include <compact_vector.hpp>
#include <vector>
#include <cstdint>

void print_vector(std::vector<int> vec);
void print_vector(std::vector<std::uint32_t> vec);
std::vector<uint32_t> range(uint32_t min, uint32_t max);
std::uint32_t argmax(std::vector<double> vec);
std::uint32_t argmin(compact::vector<int> vec);
std::uint32_t argmin(std::vector<int> vec);
