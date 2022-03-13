#include <bitset>
#include <iostream>
#include <string>
#include <time.h>
#include <random>

#include <powerset_btree.hpp>
#include "demo.hpp"
#include "demo_vector.hpp"
//#include "benchmark_bitset.hpp"
//#include "benchmark.hpp"



int main(){
    using namespace efficient_DST;
    demo();
    demo_vector();
    //benchmark_main();

/*
    const size_t N = 5000;
	std::bitset<N> b = 0;
	std::bitset<N> b2 = 111233454655211;
    clock_t t;
    //const size_t N = 5000;
    //std::bitset<N> b = 0;

    b = benchmark<N>::random_bitset((double) 0.0001);
    //b.set(4000);
    //b = pow(2, 50)-1;
    //std::cout << b << std::endl;

    std::cout << b.count() << std::endl;

    size_t current_bit;
	size_t previous_bit;
    t = clock();
    current_bit = 0;
    previous_bit = 0;
    while (current_bit != N){
    	previous_bit = current_bit;
    	current_bit = b._Find_next(current_bit);
    }
    std::cout << previous_bit+1 << " ";
    std::cout << (((float) clock() - t)*100/CLOCKS_PER_SEC) << "\n";
    t = clock();
    current_bit = b._Find_first();
    previous_bit = 0;
    while (current_bit != N){
    	previous_bit = current_bit+1;
    	current_bit = b._Find_next(current_bit);
    }
    std::cout << previous_bit << " ";
    std::cout << (((float) clock() - t)*100/CLOCKS_PER_SEC) << "\n";
    t = clock();
    current_bit = N-1;
    //previous_bit = 0;
    while (!b[current_bit] && current_bit > 0){
    	--current_bit;
    }
    std::cout << current_bit+1 << " ";
    std::cout << (((float) clock() - t)*100/CLOCKS_PER_SEC) << "\n";
    t = clock();
    b & b2;
    std::cout << (((float) clock() - t)*100/CLOCKS_PER_SEC) << "\n";
    t = clock();
    b[1530];
    std::cout << (((float) clock() - t)*100/CLOCKS_PER_SEC) << "\n";

    //////////////////////////////////////////////////////////////

    const size_t N = 10;
    size_t size = 150000;
    size_t sample_size = 100000;

    std::vector<set_N_value<double, N> > elements(size);
    std::vector<set_N_value<double, N>* > pointers;
    pointers.reserve(sample_size);
    size_t n = 0;
	while (n < sample_size){
		size_t r = rand() % elements.size();
		pointers.emplace_back(&elements[r]);
		++n;
	}

	clock_t t = clock();
	for (size_t i = 0; i < pointers.size(); ++i){
		const std::bitset<N>& b = pointers[i]->set;
	}
	t = clock() - t;
	std::cout << ((double)t)/CLOCKS_PER_SEC << std::endl;

	std::vector<size_t> pointers_i;
	pointers_i.reserve(sample_size);
    n = 0;
	while (n < sample_size){
		size_t r = rand() % elements.size();
		pointers_i.emplace_back(r);
		++n;
	}

	t = clock();
	for (size_t i = 0; i < pointers_i.size(); ++i){
		const std::bitset<N>& b = elements[pointers_i[i]].set;
	}
	t = clock() - t;
	std::cout << ((double)t)/CLOCKS_PER_SEC << std::endl;
	*/
}
