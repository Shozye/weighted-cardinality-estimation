#include"seeds.hpp"
#include <cmath>

Seeds::Seeds(const std::vector<std::uint32_t>& seeds)
: seeds_(std::ceil(std::log2(seeds.size())), seeds.size())
{
    for(size_t i = 0; i < seeds.size(); i++){
        seeds_[i] = seeds[i];
    }
}

Seeds::Seeds()
: seeds_(1, 0) 
{

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
    return std::vector<std::uint32_t>(seeds_.begin(), seeds_.end());
}