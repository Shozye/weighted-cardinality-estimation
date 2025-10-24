#include "h_fingerprint.hpp"
#include "hash_util.hpp"
#include <cmath>
#include <cstdint>

HFingerprint::HFingerprint(std::uint8_t jaccard_bits, std::size_t sketch_size)
    : H_(jaccard_bits, sketch_size), jaccard_bits(jaccard_bits) {
  std::fill(H_.begin(), H_.end(), 0);
}

HFingerprint::HFingerprint(
    std::uint8_t jaccard_bits,
    std::size_t sketch_size,
    const std::vector<std::uint32_t>& registers   
) : H_(jaccard_bits, sketch_size),
    jaccard_bits(jaccard_bits){
    if (registers.size() != sketch_size) { throw std::invalid_argument("Invalid state: registers vector size mismatch"); }
    for(size_t i = 0; i < H_.size(); i++){
        H_[i] = registers[i];
    }
}

void HFingerprint::set_elem(std::size_t index, const std::string& elem){
    const uint64_t mod = (uint64_t{1} << jaccard_bits) - 1; // =2^b - 1
    const uint64_t full_hash = murmur64(elem, 42);
    const uint32_t h = static_cast<uint32_t>((full_hash % mod) + 1U); // 1..(2^b - 1)

    this->H_[index] = h;
}

double HFingerprint::compute_jaccard(const HFingerprint& other) const {
    if (other.H_.size() != H_.size()) { return 0.0; }
    std::size_t equal = 0;
    for (std::size_t i = 0; i < H_.size(); ++i) {
        if (this->H_[i] == other.H_[i] && this->H_[i] != 0) { 
            ++equal; 
        }
    }
    double probability_h1_equals_h2 = static_cast<double>(equal) / static_cast<double>(H_.size());
    double g_max = std::pow(2, jaccard_bits) - 1;
    double jacc = (g_max * probability_h1_equals_h2 - 1)/(g_max - 1);
    if(jacc < 0){
        return 0;
    }
    return jacc;
}

std::vector<std::uint32_t> HFingerprint::get_h_registers() const {
    return std::vector<std::uint32_t>(H_.begin(), H_.end());
}

size_t HFingerprint::memory_usage_total() const {
    size_t total = H_.bytes();  // mb/8
    total += sizeof(jaccard_bits); // 1
    return total; // mb/8 + 1
}

size_t HFingerprint::memory_usage_write() const {
    size_t total = H_.bytes();  // mb/8
    return total; // mb/8
}
