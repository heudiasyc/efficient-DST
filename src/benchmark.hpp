#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <mass.hpp>
#include <commonality.hpp>
#include <belief.hpp>
#include <plausibility.hpp>
#include <pignistic_probability.hpp>
#include <conjunctive_weight.hpp>
#include <disjunctive_weight.hpp>
#include <implicability.hpp>
#include <fod.hpp>
#include <rule_conjunctive.hpp>
#include <rule_conjunctive_cautious.hpp>
#include <rule_disjunctive_bold.hpp>
#include <rule_disjunctive.hpp>
#include <rule_dempster.hpp>
#include <rule_dubois_prade.hpp>


namespace efficient_DST{

	enum class mass_family_t: bool { subnormal, non_dogmatic };

	static void generate_random_mass_function(
			const FOD& fod,
			std::vector<double>& m_vec,
			mass<double>& m,
			const mass_family_t& mass_family,
			const size_t& support_size
	){
		double mass = 1.0/support_size;

		if(mass_family == mass_family_t::non_dogmatic){
			m_vec.back() = mass;
			m.set_fod_value(mass);
		}else{
			m_vec[0] = mass;
			m.set_emptyset_value(mass);
		}
		size_t n = 1;
		while (n < support_size){
			size_t r = rand() % m_vec.size();
			if (m_vec[r] == 0) {
				m_vec[r] = mass;
				m.set_value_directly(boost::dynamic_bitset<>(fod.size(), r), mass);
				++n;
			}
		}
	}


	static void benchmark(){
		using namespace efficient_DST;

		clock_t t;
		//srand (time(NULL));
		srand (9);

		std::vector<std::string> alphabet({"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"});
		size_t subfod_size = rand() % 8 +1;
		subfod_size = 10;

		FOD fod(std::vector<std::string>(alphabet.begin(), alphabet.begin() + subfod_size));

		std::vector<double> m_vec(pow(2, fod.size()));
		mass<double> m(fod);

		size_t support_size = rand() % (m_vec.size()-1) + 1;
		support_size = 10;
		generate_random_mass_function(fod, m_vec, m, mass_family_t::non_dogmatic, support_size);

		std::cout << "\nOriginal mass values " << std::endl;
		print<>(std::cout, m.get_definition());

		std::cout << "\n============================================\n";
/*
		std::cout << "\nCommonality values from mass by FMT " << std::endl;

		t = clock();
		const std::vector<double>& q_vec = computation_scheme<double>::FMT(m_vec, fod.size(), transform_type_t::zeta, order_relation_t::superset, operation_t::addition);
		t = clock() - t;
		std::cout << "time spent = " << ((float)t)/CLOCKS_PER_SEC << " sec" << std::endl;
		std::cout << "size = " << q_vec.size() << std::endl;
		//for (size_t i = 0; i < q_vec.size(); ++i){
		//	std::cout << q_vec[i] << "\t <- " << fod.to_string(boost::dynamic_bitset<>(fod.size(), i)) << std::endl;
		//}
*/
		std::cout << "\nCommonality values from mass by EMT " << std::endl;

		t = clock();
		zeta_transform<double> q(m.get_definition(), order_relation_t::superset, operation_t::addition, scheme_type_t::lattice);
		const powerset_btree<double>& q_def = q.get_definition();
		//const powerset_btree<double>& q_def = computation_scheme<double>::EMT_with_lattice(m.get_definition(), transform_type_t::zeta, order_relation_t::superset, operation_t::addition);
		t = clock() - t;
		std::cout << "time spent = " << ((float)t)/CLOCKS_PER_SEC << " sec" << std::endl;
		std::cout << "size = " << q_def.size() << std::endl;
		//print<>(std::cout, q_def);
/*
		std::cout << "\n============================================\n";

		std::cout << "\nMass values from commonality by FMT " << std::endl;

		//const std::vector<double>& m_back_vec = computation_scheme<double>::FMT(q_vec, fod.size(), transform_type_t::Mobius, order_relation_t::superset, operation_t::addition);
		//for (size_t i = 0; i < m_back_vec.size(); ++i){
		//	std::cout << m_back_vec[i] << "\t <- " << fod.to_string(boost::dynamic_bitset<>(fod.size(), i)) << std::endl;
		//}

		std::cout << "\nMass values from commonality by EMT " << std::endl;

		t = clock();
		const powerset_btree<double>& m_back_def = q.inversion(operation_t::addition);
		//const powerset_btree<double>& m_back_def = computation_scheme<double>::EMT_with_lattice(q_def, transform_type_t::Mobius, order_relation_t::superset, operation_t::addition);
		t = clock() - t;
		std::cout << "time spent = " << ((float)t)/CLOCKS_PER_SEC << " sec" << std::endl;
		std::cout << "size = " << m_back_def.size() << std::endl;
		print<>(std::cout, m_back_def);*/
	}
}
