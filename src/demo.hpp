#include <iostream>
#include <string>
#include <vector>

#include <commonality_function.hpp>
#include <conjunctive_decomposition.hpp>
#include <disjunctive_decomposition.hpp>
#include <implicability_function.hpp>
#include <mass_function.hpp>
#include <plausibility_function.hpp>
#include <rule_conjunctive.hpp>
//#include <rule_conjunctive_cautious.hpp>
//#include <rule_disjunctive_bold.hpp>
//#include <rule_disjunctive.hpp>
//#include <rule_dempster.hpp>
//#include <rule_dubois_prade.hpp>
#include <powerset_vector.hpp>


void demo(){
    using namespace efficient_DST;
    typedef float T;
    scheme_type_t scheme_type = scheme_type_t::lattice;
    const bool adaptive_uncertainty = true;

    const size_t N = 8;
    std::string labels[N];
	for (size_t i = 0; i < N; ++i){
		labels[i] = std::to_string(i);
	}
    sample_space<N> outcomes(labels);

    conjunctive_decomposition<N, T> w0(outcomes, adaptive_uncertainty);

    std::cout << "\n============================================\n";

    w0.assign_emptyset(0.82);
    w0.assign({"1"}, 0.67);
    w0.assign({"1", "2"}, 0.14286);
    w0.assign({"1", "5"}, 0.02381);
    w0.assign({"1", "6"}, 0.03571);
    w0.assign({"0", "1", "4", "7"}, 0.125);

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

//
//    std::string labels_vec[] = {"a", "b", "c", "d"};
//    FOD<4> fod_vec(labels_vec);
//
//    //q_vec corresponds to the mass function {0.01, 0.06, 0, 0, 0.38, 0, 0.05, 0, 0, 0, 0, 0.05, 0.27, 0, 0, 0.18}.
//    std::vector<double> q_vec = {1, 0.29, 0.28, 0.23, 0.88, 0.18, 0.23, 0.18, 0.5, 0.23, 0.23, 0.23, 0.45, 0.18, 0.18, 0.18};

//    commonality<double, 4> q_from_vec(q_vec, fod_vec);
//
//    std::cout << "\nCommonality values " << std::endl;
//
//    q_from_vec.get_definition().print(std::cout);
//
//	std::cout << "\n============================================\n";
//
//	//mass<> m_back_vec(q_from_vec);
//
//	std::cout << "\nMass values from commonality " << std::endl;
//
//	//const powerset_btree<double>& m_back_vec = computation_scheme<double>::EMT_with_lattice_Mobius_from_zeta_values(q_vec, fod_vec, transform_type_t::Mobius, order_relation_t::superset, operation_t::addition);
//	const powerset_btree<double, 4>& m_back_vec = computation_scheme<double, 4>::EMT_with_semilattice(q_from_vec.get_definition(), transform_type_t::Mobius, order_relation_t::superset, operation_t::addition);
//	m_back_vec.print(std::cout);
//	//print<>(std::cout, m_back_vec.get_definition());
//
//	std::cout << "\n============================================\n";
//
//    std::cout << "\nMass values from commonality by FMT " << std::endl;
//
//    //const std::vector<double>& m_values = computation_scheme<double, 4>::FMT(q_vec, transform_type_t::Mobius, order_relation_t::superset, operation_t::addition);
//	//for (size_t i = 0; i < m_values.size(); ++i){
//	//	std::cout << m_values[i] << "\t <- " << fod_vec.to_string(std::bitset<4>(i)) << std::endl;
//	//}
//    powerset_btree<double, 4> q_values(&fod_vec, q_vec.size());
//	for (size_t i = 0; i < q_vec.size(); ++i){
//		q_values.insert(std::bitset<4>(i), q_vec[i]);
//	}
//    const powerset_btree<double, 4>& m_values = computation_scheme<double, 4>::FMT(q_values, transform_type_t::Mobius, order_relation_t::superset, operation_t::addition);
//    m_values.print(std::cout);
//
//	std::cout << "\n============================================\n";
//
//	std::cout << "\nImplicability vector " << std::endl;
//	std::reverse(q_vec.begin(), q_vec.end());
//	for (size_t i = 0; i < q_vec.size(); ++i){
//		std::cout << q_vec[i] << "\t <- " << fod_vec.to_string(std::bitset<4>(i)) << std::endl;
//	}
//
//	implicability<double, 4> b_from_vec(q_vec, fod_vec);
//
//    std::cout << "\nImplicability values " << std::endl;
//
//    b_from_vec.get_definition().print(std::cout);
//
//	std::cout << "\n============================================\n";
//
//	//mass<> m_back_vec_2(b_from_vec);
//
//	std::cout << "\nMass values from implicability " << std::endl;
//
//	//const powerset_btree<double>& m_back_vec_2 = computation_scheme<double>::EMT_with_lattice_Mobius_from_zeta_values(q_vec, fod_vec, transform_type_t::Mobius, order_relation_t::subset, operation_t::addition);
//	const powerset_btree<double, 4>& m_back_vec_2 = computation_scheme<double, 4>::EMT_with_semilattice(b_from_vec.get_definition(), transform_type_t::Mobius, order_relation_t::subset, operation_t::addition);
//	m_back_vec_2.print(std::cout);
//	//print<>(std::cout, m_back_vec_2.get_definition());
//
//	std::cout << "\n============================================\n";
//
//    std::cout << "\nMass values from implicability by FMT " << std::endl;
//
//    //const std::vector<double>& m_values_2 = computation_scheme<double, 4>::FMT(q_vec, transform_type_t::Mobius, order_relation_t::subset, operation_t::addition);
//	//for (size_t i = 0; i < m_values_2.size(); ++i){
//	//	std::cout << m_values_2[i] << "\t <- " << fod_vec.to_string(std::bitset<4>(i)) << std::endl;
//	//}
//    powerset_btree<double, 4> b_values(&fod_vec, q_vec.size());
//	for (size_t i = 0; i < q_vec.size(); ++i){
//		b_values.insert(std::bitset<4>(i), q_vec[i]);
//	}
//    const powerset_btree<double, 4>& m_values_2 = computation_scheme<double, 4>::FMT(b_values, transform_type_t::Mobius, order_relation_t::subset, operation_t::addition);
//    m_values_2.print(std::cout);
//
//	std::cout << "\n============================================\n";
//
//    plausibility<double, N> pl(m);
//
//	std::cout << "\nPlausibility contour from mass " << std::endl;
//
//	std::vector<double> contour_pl = pl.get_contour();
//	for (size_t i = 0; i < contour_pl.size(); ++i){
//		std::cout << contour_pl[i] << "\t <- {" << fod.elements()[i].label << "}" << std::endl;
//	}
//
//	std::cout << "\n============================================\n";
//
//    pignistic_probability<double, N> bet_p(m);
//
//	std::cout << "\nPignistic probability contour from mass " << std::endl;
//
//	std::vector<double> contour_bet = bet_p.get_contour();
//	for (size_t i = 0; i < contour_bet.size(); ++i){
//		std::cout << contour_bet[i] << "\t <- {" << fod.elements()[i].label << "}" << std::endl;
//	}
//
//    std::cout << "\n============================================\n";
//
//    commonality<double, N> q(m);
//
//    std::cout << "\nCommonality values " << std::endl;
//
//    q.get_definition().print(std::cout);
//
//	std::cout << "\n============================================\n";
//
//	std::cout << q[{"f", "i", "l"}] << std::endl;
//	std::cout << q[{"f", "j"}] << std::endl;
//	std::cout << q[{"f", "i", "o"}] << std::endl;
//	std::cout << q[{"h", "o"}] << std::endl;
//
//	mass<double, N> m_back(q);
//
//	std::cout << "\nMass values from commonality " << std::endl;
//
//    m_back.get_definition().print(std::cout);
//
//	std::cout << "\n============================================\n";
//
//	//std::clog << q[{"j", "l"}];
//
//	conjunctive_weight<double, N> w(q);
//
//	std::cout << "\nConjunctive weight values " << std::endl;
//
//    w.get_definition().print(std::cout);
//
//	std::cout << "\n============================================\n";
//
//    commonality<double, N> q2(w);
//
//    std::cout << "\nCommonality values from weights " << std::endl;
//
//    q2.get_definition().print(std::cout);
//
//	std::cout << "\n============================================\n";
//
//	mass<double, N> m2(q2);
//
//	std::cout << "\nMass values from conjunctive weights " << std::endl;
//
//    m2.get_definition().print(std::cout);
//
//	std::cout << "\n============================================\n";
//
//
//	implicability<double, N> b(m);
//
//	std::cout << "\nImplicability values " << std::endl;
//
//	b.get_definition().print(std::cout);
//
//	std::cout << "\n============================================\n";
//
//	mass<double, N> m_back2(b);
//
//	std::cout << "\nMass values from implicability " << std::endl;
//
//    m_back2.get_definition().print(std::cout);
//
//	std::cout << "\n============================================\n";
//
//	disjunctive_weight<double, N> v(b);
//
//	std::cout << "\nDisjunctive weight values " << std::endl;
//
//    v.get_definition().print(std::cout);
//
//	std::cout << "\n============================================\n";
//
//    implicability<double, N> b2(v);
//
//    std::cout << "\nImplicability values from weights " << std::endl;
//
//    b2.get_definition().print(std::cout);
//
//	std::cout << "\n============================================\n";
//
//	mass<double, N> m3(b2);
//
//	std::cout << "\nMass values from disjunctive weights " << std::endl;
//
//    m3.get_definition().print(std::cout);
//
//	std::cout << "\n============================================\n";
//
//	conjunctive_weight<double, N> w4(w);
//
//    conjunctive_weight<double, N> w14 = w.apply<rule_Dempster<double, N> >(w4);
//
//	std::cout << "\nMass from conjunctive weight fusion values " << std::endl;
//
//    mass<double, N>(commonality<double, N>(w14)).get_definition().print(std::cout);
//
//	std::cout << "\n============================================\n";
//
//	commonality<double, N> q4(q);
//
//	commonality<double, N> q14 = q.apply<rule_conjunctive_cautious<double, N> >(q4);
//
//	std::cout << "\nMass from commonality fusion values " << std::endl;
//
//    mass<double, N>(q14).get_definition().print(std::cout);
//
//	std::cout << "\n============================================\n";
//
//	mass<double, N> m4(m);
//
//	mass<double, N> m14 = m.apply<rule_Dubois_Prade<double, N> >(m4);
//
//	std::cout << "\nMass fusion values " << std::endl;
//
//    m14.get_definition().print(std::cout);
//
//	std::cout << "\n============================================\n";
}
