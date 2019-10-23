#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <boost/dynamic_bitset.hpp>
#include <vector>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <iomanip>

const size_t s = 10000;

void b_and(std::bitset<s> b){
	b & b;
}


void benchmark_bitset(){

	srand (time(NULL));
	//srand (9);
	size_t n = rand() % (1000000000) + 1;
    clock_t t;
    t = clock();
    std::bitset<s> b = n;
    t = clock() - t;
    std::cout << "time spent = " << ((float)t)/CLOCKS_PER_SEC << " sec" << std::endl;
    t = clock();
    boost::dynamic_bitset<> b_dyn(s, n);
    t = clock() - t;
    std::cout << "time spent = " << ((float)t)/CLOCKS_PER_SEC << " sec" << std::endl;


    t = clock();
    for (size_t i = 0; i < b.size(); ++i){
    	std::cout << b[i];
    }
    t = clock() - t;
    std::cout << "\ntime spent = " << ((float)t)/CLOCKS_PER_SEC << " sec" << std::endl;


    t = clock();
    for (size_t i = 0; i < b.size(); ++i){
    	std::cout << b_dyn[i];
    }
    t = clock() - t;
    std::cout << "\ntime spent = " << ((float)t)/CLOCKS_PER_SEC << " sec" << std::endl;

    std::cout << "================== bitwise operations\n";

    t = clock();
    for (size_t i = 0; i < b.size(); ++i){
    	b_and(b);
    }
    t = clock() - t;
    std::cout << "time spent = " << ((float)t)/CLOCKS_PER_SEC << " sec" << std::endl;


    t = clock();
    for (size_t i = 0; i < b.size(); ++i){
    	b_dyn & b_dyn;
    }
    t = clock() - t;
    std::cout << "time spent = " << ((float)t)/CLOCKS_PER_SEC << " sec" << std::endl;
}
