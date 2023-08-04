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

#include <commonality_function.hpp>
#include <conjunctive_decomposition.hpp>
#include <disjunctive_decomposition.hpp>
#include <implicability_function.hpp>
#include <mass_function.hpp>
#include <plausibility_function.hpp>
#include <powerset_vector.hpp>
#include <benchmarking.hpp>


void demo(){
    using namespace efficient_DST;
    typedef float T;
    scheme_type_t scheme_type = scheme_type_t::semilattice;
    const bool adaptive_uncertainty = true;

    const size_t N = 20;
    std::string labels[N];
	for (size_t i = 0; i < N; ++i){
		labels[i] = std::to_string(i);
	}
    sample_space<N> outcomes(labels);

    conjunctive_decomposition<N, T> w0(outcomes, adaptive_uncertainty);

    std::cout << "\n============================================\n";

//    w0.assign_emptyset(0.82);
//    w0.assign({"1"}, 0.67);
//    w0.assign({"1", "2"}, 0.14286);
//    w0.assign({"1", "5"}, 0.02381);
//    w0.assign({"1", "6"}, 0.03571);
//    w0.assign({"0", "1", "4", "7"}, 0.125);

    mass_function<N, T> m000(outcomes);

    std::cout << "Mass generation\n";
    int seed = -1;
    benchmarking<N, T>().generate_random_mass_function(
		m000,
		mass_family_t::random,
		1,
		6,
		seed
	);
    const std::vector<set_N_value<N, T> const * >& elements = m000.get_definition().elements();
    for (size_t i = 0; i < elements.size(); ++i){
    	w0.assign(elements[i]->set, elements[i]->value);
    }

    std::cout << "\nConjunctive decomposition in w0" << std::endl;

    w0.print();

    std::cout << "\n============================================\n";

    commonality_function<N, T> q0(w0, scheme_type);

    std::cout << "\nCommonality values from w0" << std::endl;

    q0.print();

    std::cout << "\n============================================\n";

    mass_function<N, T> m0(q0);

    std::cout << "\nMass values from w0" << std::endl;

    m0.print();

    std::cout << "\n============================================\n";

    mass_function<N, T> m00(w0);

    std::cout << "\nMass values directly from w0" << std::endl;

    m00.print();

    std::cout << "\n============================================\n";

    disjunctive_decomposition<N, T> v0(w0, adaptive_uncertainty);

    std::cout << "\n============================================\n";

    std::cout << "\nDisjunctive decomposition in v0" << std::endl;

    v0.print();

    std::cout << "\n============================================\n";

    implicability_function<N, T> b0(v0, scheme_type);

    std::cout << "\nImplicability values from v0" << std::endl;

    b0.print();

    std::cout << "\n============================================\n";

    mass_function<N, T> mv0(b0);

    std::cout << "\nMass values from v0" << std::endl;

    mv0.print();

	std::cout << "\n============================================\n";

    mass_function<N, T> mv00(v0);

    std::cout << "\nMass values directly from v0" << std::endl;

    mv00.print();

    std::cout << "\n============================================\n";

	mass_function<N, T> m(outcomes);

    m.assign_emptyset(0.42);
//    m.assign({"f"}, 0.4);
    m.assign({"1", "2"}, 0.08);
    m.assign({"1", "6"}, 0.03);
    m.assign({"0", "1", "4", "7"}, 0.37);
    m.assign({"1", "4"}, 0.05);
    m.assign({"0", "1", "2", "4", "5", "6", "7"}, 0.05);

	std::cout << "\nMass values m" << std::endl;

	m.print();

    std::cout << "\n============================================\n";

    commonality_function<N, T> q(m, scheme_type);

    std::cout << "\nCommonality values from m" << std::endl;

    q.print();

    std::cout << "\n============================================\n";

    conjunctive_decomposition<N, T> w(q, adaptive_uncertainty);

    std::cout << "\nConjunctive decomposition in m" << std::endl;

    w.print();

	std::cout << "\n============================================\n";

    mass_function<N, T> mw(w);

    std::cout << "\nMass values directly from w" << std::endl;

    mw.print();

    std::cout << "\n============================================\n";

    implicability_function<N, T> b(m, scheme_type);

    std::cout << "\nImplicability values from m" << std::endl;

    b.print();

    std::cout << "\n============================================\n";

    disjunctive_decomposition<N, T> v(b, adaptive_uncertainty);

    std::cout << "\nDisjunctive decomposition in m" << std::endl;

    v.print();

	std::cout << "\n============================================\n";

    mass_function<N, T> mv(v);

    std::cout << "\nMass values directly from v" << std::endl;

    mv.print();
}
