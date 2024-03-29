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

#include <commonality_vector.hpp>
#include <conjunctive_decomposition_vector.hpp>
#include <disjunctive_decomposition_vector.hpp>
#include <implicability_vector.hpp>
#include <mass_vector.hpp>
#include <plausibility_vector.hpp>
#include <powerset_vector.hpp>


void demo_vector(){
    using namespace efficient_DST;
    typedef float T;
    bool core_reduced = true;
    bool adaptive_uncertainty = true;

    const size_t N = 8;
    std::string labels[] = {"e", "f", "g", "h", "i", "j", "l", "o"};
    sample_space<N> outcomes(labels);

    conjunctive_decomposition_vector<N, T> w0(outcomes, adaptive_uncertainty);

    std::cout << "\n============================================\n";

    w0.assign_emptyset(0.82);
    w0.assign({"f"}, 0.67);
    w0.assign({"f", "g"}, 0.14286);
    w0.assign({"f", "j"}, 0.02381);
    w0.assign({"f", "l"}, 0.03571);
    w0.assign({"e", "f", "i", "o"}, 0.125);

    std::cout << "\nConjunctive decomposition_vector in w0" << std::endl;

    w0.print();

    std::cout << "\n============================================\n";

    commonality_vector<N, T> q0(w0, core_reduced);

    std::cout << "\nCommonality values from w0" << std::endl;

    q0.print();

    std::cout << "\n============================================\n";

    mass_vector<N, T> m0(q0);

    std::cout << "\nMass values from w0" << std::endl;

    m0.print();

    std::cout << "\n============================================\n";

    mass_vector<N, T> m00(w0);

    std::cout << "\nMass values directly from w0" << std::endl;

    m00.print();

    std::cout << "\n============================================\n";

    disjunctive_decomposition_vector<N, T> v0(w0, adaptive_uncertainty);

    std::cout << "\n============================================\n";

    std::cout << "\nDisjunctive decomposition_vector in v0" << std::endl;

    v0.print();

    std::cout << "\n============================================\n";

    implicability_vector<N, T> b0(v0, core_reduced);

    std::cout << "\nImplicability values from v0" << std::endl;

    b0.print();

    std::cout << "\n============================================\n";

    mass_vector<N, T> mv0(b0);

    std::cout << "\nMass values from v0" << std::endl;

    mv0.print();

	std::cout << "\n============================================\n";

    mass_vector<N, T> mv00(v0);

    std::cout << "\nMass values directly from v0" << std::endl;

    mv00.print();

    std::cout << "\n============================================\n";

	mass_vector<N, T> m(outcomes);

    m.assign_emptyset(0.42);
//    m.assign({"f"}, 0.4);
    m.assign({"f", "g"}, 0.08);
    m.assign({"f", "l"}, 0.03);
    m.assign({"e", "f", "i", "o"}, 0.37);
    m.assign({"f", "i"}, 0.05);
    m.assign({"e", "f", "g", "i", "j", "l", "o"}, 0.05);

	std::cout << "\nMass values m" << std::endl;

	m.print();

    std::cout << "\n============================================\n";

    commonality_vector<N, T> q(m, core_reduced);

    std::cout << "\nCommonality values from m" << std::endl;

    q.print();

    std::cout << "\n============================================\n";

    conjunctive_decomposition_vector<N, T> w(q, adaptive_uncertainty);

    std::cout << "\nConjunctive decomposition_vector in m" << std::endl;

    w.print();

	std::cout << "\n============================================\n";

    mass_vector<N, T> mw(w);

    std::cout << "\nMass values directly from w" << std::endl;

    mw.print();

    std::cout << "\n============================================\n";

    implicability_vector<N, T> b(m, core_reduced);

    std::cout << "\nImplicability values from m" << std::endl;

    b.print();

    std::cout << "\n============================================\n";

    disjunctive_decomposition_vector<N, T> v(b, adaptive_uncertainty);

    std::cout << "\nDisjunctive decomposition_vector in m" << std::endl;

    v.print();

	std::cout << "\n============================================\n";

    mass_vector<N, T> mv(v);

    std::cout << "\nMass values directly from v" << std::endl;

    mv.print();
}
