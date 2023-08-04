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

#include <bitset>
#include <iostream>
#include <string>
#include <time.h>
#include <random>

#include <powerset_btree.hpp>
#include "demo.hpp"
#include "demo_vector.hpp"
#include <benchmarking.hpp>


int main(){
    using namespace efficient_DST;
    demo();
//    demo_vector();
//	typedef float T;
//	benchmarking<9, T>().run();
//	benchmarking<10, T>().run();
//	benchmarking<11, T>().run();
//	benchmarking<12, T>().run();
//	benchmarking<13, T>().run();
//	benchmarking<14, T>().run();
//	benchmarking<15, T>().run();
//	benchmarking<16, T>().run();
//	benchmarking<17, T>().run();
//	benchmarking<18, T>().run();
//	benchmarking<19, T>().run();
//	benchmarking<20, T>().run();
//	benchmarking<21, T>().run();
//	benchmarking<22, T>().run();
//	benchmarking<23, T>().run();
//	benchmarking<24, T>().run();
//	benchmarking<25, T>().run();
//	benchmarking<26, T>().run();

	std::cout << "--- Done." << std::endl;
}
