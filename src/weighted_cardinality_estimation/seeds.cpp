#include"seeds.hpp"
#include <cmath>

Seeds::Seeds(const std::vector<std::uint32_t>& seeds)
: seeds_(create_seeds_vector(seeds))
{
}

Seeds::Seeds()
: seeds_(1, 0) 
{

}

compact::vector<uint32_t> Seeds::create_seeds_vector(const std::vector<std::uint32_t>& seeds){
    if(!is_sequential_from_one(seeds)){
        const auto amount_bits = static_cast<std::uint32_t>(std::ceil(std::log2(seeds.size())));
        compact::vector<uint32_t> seeds_(amount_bits, seeds.size());
        for(size_t i = 0; i < seeds.size(); i++){
            seeds_[i] = seeds[i];
        }
        return seeds_;
    }
    return compact::vector<uint32_t> (1, 0);
   
}

std::uint32_t Seeds::get(uint32_t index) const {
    if(this->seeds_.empty()){
        return index + 1; // i don't like 0 being seed.
    }          
    return this->seeds_.at(index);
}

std::uint32_t Seeds::bytes() const {
    return this->seeds_.bytes();
}

std::uint32_t Seeds::operator[](uint32_t index) const {
    return this->get(index);
}

std::vector<std::uint32_t> Seeds::toVector() const {
    if (this->seeds_.empty()){
        return std::vector<std::uint32_t>();
    }
    return std::vector<std::uint32_t>(seeds_.begin(), seeds_.end());
}

bool Seeds::is_sequential_from_one(const std::vector<std::uint32_t>& vec){
    for(size_t i = 0; i < vec.size(); ++i){
        if (vec[i] != i+1){
            return false;
        }
    }
    return true;
};