/*
 * Copyright (C) 2019-2023  Maxime Chaveroche (maxime.chaveroche@gmail.com)
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the CeCILL License, either version 2.1 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * CeCILL License for more details.
 * 
 * You should have received a copy of the CeCILL License
 * along with this program. If not, see <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html>.
 */

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
