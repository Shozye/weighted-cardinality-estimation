#include <compact_vector.hpp>
#include <vector>
#include <cstdint>
class HFingerprint {
public:
    HFingerprint(
        std::uint8_t jaccard_bits,
        std::size_t sketch_size
    );
    HFingerprint(
        std::uint8_t jaccard_bits,
        std::size_t sketch_size,
        const std::vector<std::uint32_t>& registers
    );

    void set_elem(std::size_t index, const std::string& elem);
    double compute_jaccard(const HFingerprint& other, std::vector<int> v1, std::vector<int> v2) const;
    std::vector<std::uint32_t> get_h_registers() const;

    [[nodiscard]] size_t memory_usage_total() const;
    [[nodiscard]] size_t memory_usage_write() const;
private:
    compact::vector<std::uint32_t> H_;
    std::uint8_t jaccard_bits;
};
