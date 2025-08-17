#include <vector>
#include <iostream>
#include <cstdint>
#include<numeric>


int mypow(int a, int b){
    int res = 1;
    while(b > 0){
        res *= a;
        b--;
    }
    return res;
}
void print_vector(std::vector<int> vec){
    std::cout << "vec:[";
    for(size_t i = 0; i < vec.size();i++){
        std::cout << vec[i] << ", ";
    }
    std::cout << "]\n";
}

std::vector<uint32_t> range(uint32_t min, uint32_t max){
    std::vector<std::uint32_t> vec(max - min + 1);
    std::iota(vec.begin(), vec.end(), min); 
    return vec;
}