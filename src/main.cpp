#include "demo.hpp"
#include "benchmark.hpp"
#include "benchmark_bitset.hpp"


int main(){
    using namespace efficient_DST;
    demo();
    //benchmark();
    //benchmark_bitset();

    /*
    //std::bitset<5> b = 4;
    boost::dynamic_bitset<> b(5, 4);
    clock_t t;
    t = clock();
    b.to_ulong();
    t = clock() - t;
    std::cout << "time spent = " << ((float)t)/CLOCKS_PER_SEC << " sec" << std::endl;
    */
}
