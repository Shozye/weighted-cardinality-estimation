#include <vector>
#include <iostream>


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
    std::cout << "]" << std::endl;
}
