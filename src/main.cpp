#include "demo.hpp"
#include "benchmark.hpp"
#include "benchmark_bitset.hpp"


int main(){
    using namespace efficient_DST;
    //demo();
    benchmark();
    //benchmark_bitset();
/*
    const size_t N = 8;
    std::string fod_definition[N] = {"e", "f", "g", "h", "i", "j", "l", "o"};
    //FOD fod(fod_definition);
    //FOD<N>::set b = 4;
    std::bitset<N> b = 4;
    //boost::dynamic_bitset<> b(5, 4);
    clock_t t;
    t = clock();
    b.to_ulong();
    t = clock() - t;
    std::cout << "time spent = " << ((float)t)/CLOCKS_PER_SEC << " sec" << std::endl;
    */
}
