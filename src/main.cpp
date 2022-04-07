#include <bitset>
#include <iostream>
#include <string>
#include <time.h>
#include <random>

#include <powerset_btree.hpp>
#include "demo.hpp"
#include "demo_vector.hpp"
//#include "benchmark_bitset.hpp"
#include <benchmarking.hpp>



int main(){
    using namespace efficient_DST;
    demo();
//    demo_vector();
	typedef float T;
	//constexpr static size_t fod_sizes[] = {16, 17, 20, 22, 24, 25, 26};
	//constexpr static size_t fod_sizes[] = {50, 100, 200, 400, 800, 1600};
//		srand(time(NULL));

//		for (size_t n = 0; n < nb_N; ++n){
//		benchmark<9, T>().run();
//		benchmark<12, T>().run();
//		benchmark<15, T>().run();
//		benchmark<18, T>().run();
//	benchmarking<21, T>().run();
//		benchmark<24, T>().run();
//		benchmark<27, T>().run();
//		benchmark<30, T>().run();
	std::cout << "--- Done." << std::endl;
}
