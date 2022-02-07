#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <boost/functional/hash.hpp>
#include <random>


#include <mass.hpp>
#include <commonality.hpp>
#include <belief.hpp>
#include <plausibility.hpp>
#include <pignistic_probability.hpp>
#include <conjunctive_weight.hpp>
#include <disjunctive_weight.hpp>
#include <implicability.hpp>
#include <rule_conjunctive.hpp>
#include <rule_conjunctive_cautious.hpp>
#include <rule_disjunctive_bold.hpp>
#include <rule_disjunctive.hpp>
#include <rule_dempster.hpp>
#include <rule_dubois_prade.hpp>


namespace efficient_DST{

	enum class mass_family_t: int8_t { random, almost_consonant, almost_bayesian };

	template <size_t N>
	class benchmark {
	protected:
		static const std::string alphabet[26];

	public:

		static inline std::bitset<N> random_bitset( double p = 1/N, size_t max_card = N) {

		std::bitset<N> bits;
		std::random_device rd;
		std::mt19937 gen( rd());
		std::bernoulli_distribution d( p);

		//size_t n = rand() % N;
		size_t n = 0;
		size_t count = 0;
		for(; n < N; ++n) {
			bool val = d(gen);
			bits[n] = val;
			if (val)
				++count;
			if (count >= max_card)
				break;
		}

		return bits;
		}

		static void generate_random_mass_function(
				std::vector<double>& m_vec,
				mass<double, N>& m,
				const mass_family_t& mass_family,
				const size_t& support_size
		){
			const size_t powerset_size = pow(2, N);
			size_t size = std::min(support_size, powerset_size);
			double mass;
			float consonant_proportion = 0.8;
			float bayesian_proportion = 0.8;
			float special_support_size_f = 0;
			size_t special_support_size = 0;
			std::bitset<N> set(0);
			size_t n = 0;
			switch (mass_family) {
				case mass_family_t::almost_consonant:

					special_support_size_f = std::min((float) consonant_proportion * size, (float) N+1);
					special_support_size = (size_t) special_support_size_f;
					size = (size_t) (special_support_size_f / consonant_proportion);
					mass = 1.0/size;

					while (n < special_support_size){
						m_vec[set.to_ulong()] = mass;
						m.set_value_directly(set, mass);
						if(n < N)
							set[n] = true;
						++n;
					}
					break;
				case mass_family_t::almost_bayesian:

					special_support_size_f = std::min((float) bayesian_proportion * size, (float) N);
					special_support_size = (size_t) special_support_size_f;
					size = (size_t) (special_support_size_f / bayesian_proportion);
					mass = 1.0/size;

					while (n < special_support_size){
						set = 0;
						set[n] = true;
						m_vec[set.to_ulong()] = mass;
						m.set_value_directly(set, mass);
						++n;
					}
					break;
				default:
					mass = 1.0/size;
					break;
			}
			if (mass == 0){
				std::cerr << "Number of elements in support exceeding machine precision\n";
				exit(1);
			}
			/*
			if(mass_family == mass_family_t::non_dogmatic){
				m_vec.back() = mass;
				m.set_fod_value(mass);
				++n;
			}else if(mass_family == mass_family_t::subnormal) {
				m_vec[0] = mass;
				m.set_emptyset_value(mass);
				++n;
			}*/
			while (n < size){
				size_t r = rand() % m_vec.size();
				if (m_vec[r] == 0) {
					m_vec[r] = mass;
					m.set_value_directly(std::bitset<N>(r), mass);
					++n;
				}
			}
		}

		static void generate_random_mass_function(
				mass<double, N>& m,
				const mass_family_t& mass_family,
				const size_t& support_size
		){
			size_t size = support_size;
			double mass;
			float consonant_proportion = 0.8;
			float bayesian_proportion = 0.8;
			float special_support_size_f = 0;
			size_t special_support_size = 0;
			std::bitset<N> set(0);
			size_t n = 0;
			switch (mass_family) {
				case mass_family_t::almost_consonant:

					special_support_size_f = std::min((float) consonant_proportion * size, (float) N+1);
					special_support_size = (size_t) special_support_size_f;
					size = (size_t) (special_support_size_f / consonant_proportion);
					mass = 1.0/size;

					while (n < special_support_size){
						m.set_value_directly(set, mass);
						if(n < N)
							set[n] = true;
						++n;
					}
					break;
				case mass_family_t::almost_bayesian:

					special_support_size_f = std::min((float) bayesian_proportion * size, (float) N);
					special_support_size = (size_t) special_support_size_f;
					size = (size_t) (special_support_size_f / bayesian_proportion);
					mass = 1.0/size;

					while (n < special_support_size){
						set = 0;
						set[n] = true;
						m.set_value_directly(set, mass);
						++n;
					}
					break;
				default:
					mass = 1.0/size;
					break;
			}
			if (mass == 0){
				std::cerr << "Number of elements in support exceeding machine precision\n";
				exit(1);
			}
/*
			std::unordered_set<std::pair<size_t, size_t>, boost::hash<std::pair<size_t, size_t> > > generated_sets;
			generated_sets.reserve(support_size);
			const size_t powerset_size = pow(2, std::min(N, (size_t) 32));
			while (n < size){
				size_t r = rand() % powerset_size;
				size_t offset = rand() % N;
				std::bitset<N> set(r);
				r <<= offset;
				std::pair<size_t, size_t> key(r, offset);
				bool insertion = generated_sets.emplace(key).second;
				if (insertion) {
					m.set_value_directly(set, mass);
					++n;
				}
			}*/
			std::unordered_set<std::bitset<N>> generated_sets;
			generated_sets.reserve(support_size);
			while (n < size){
				const std::bitset<N>& random_set = random_bitset();
				bool insertion = generated_sets.emplace(random_set).second;
				if (insertion) {
					m.set_value_directly(random_set, mass);
					++n;
				}
			}
		}

		static void one_step_EMT(
				order_relation_t order_relation,
				operation_t operation,
				scheme_type_t scheme_type,
				const mass<double, N>& m,
				bool build_persistence,
				clock_t& t,
				size_t& count_z,
				size_t& count_m,
				size_t& avg_semilattice_support_size
		){
			if (build_persistence){
				t = clock();
				zeta_transform<double, N> z(m.get_definition(), order_relation, operation, scheme_type);
				t = clock() - t;
				count_z += t;
				//z.get_definition().print(std::cout);
				avg_semilattice_support_size += z.get_definition().size();

				if (operation == operation_t::addition){
					t = clock();
					const mass<double, N>& m_back_def(z);
					t = clock() - t;
					count_m += t;
					//m_back_def.get_definition().print(std::cout);
				}else{
					if (order_relation == order_relation_t::subset){
						t = clock();
						const disjunctive_weight<double, N>& m_back_def(z);
						t = clock() - t;
						count_m += t;
					}else{
						t = clock();
						const conjunctive_weight<double, N>& m_back_def(z);
						t = clock() - t;
						count_m += t;
					}
				}
			}else{
				const powerset_btree<double, N>& m_def = m.get_definition();
				powerset_btree<double, N> z_def(m_def.get_FOD(), m_def.get_FOD_size() * m_def.size());
				std::vector<std::bitset<N> > iota_sequence;
				iota_sequence.reserve(N);
				t = clock();
				computation_scheme<double, N>::build_and_execute(
						m_def,
						z_def,
						transform_type_t::zeta,
						order_relation,
						operation,
						iota_sequence,
						scheme_type
				);
				t = clock() - t;
				count_z += t;
				avg_semilattice_support_size += z_def.size();

				iota_sequence.clear();
				powerset_btree<double, N> m_back_def(m_def.get_FOD(), m_def.get_FOD_size() * m_def.size());
				t = clock();
				computation_scheme<double, N>::build_and_execute(
						z_def,
						m_back_def,
						transform_type_t::Mobius,
						order_relation,
						operation,
						iota_sequence,
						scheme_type
				);
				t = clock() - t;
				count_m += t;
			}
		}

		static void run(order_relation_t order_relation, operation_t operation, scheme_type_t scheme_type, size_t support_size, mass_family_t mass_family, bool build_persistence, size_t nb_of_tests, std::vector<std::string>& out) {
			std::cout << "--- Benchmark with FOD of size " << N << std::endl;
			clock_t t;
			srand (time(NULL));

			std::string sub_fod_labels[N];
			for (size_t i = 0; i < N; ++i){
				sub_fod_labels[i] = alphabet[i];
			}
			FOD<N> fod(sub_fod_labels);

			std::string line = std::to_string(N);

			std::vector<double> m_vec(pow(2, N));
			mass<double, N> m(fod);

			generate_random_mass_function(m_vec, m, mass_family, support_size);

			//std::cout << "\nOriginal mass values " << std::endl;
			//m.get_definition().print(std::cout);

			t = clock();
			//const std::vector<double>& q_vec = computation_scheme<double, N>::FMT(m_vec, transform_type_t::zeta, order_relation, operation);
			const powerset_btree<double, N>& q = computation_scheme<double, N>::FMT(m.get_definition(), transform_type_t::zeta, order_relation, operation);
			t = clock() - t;
			line += "," + std::to_string(((double)t)/CLOCKS_PER_SEC);

			t = clock();
			//const std::vector<double>& m_back_vec = computation_scheme<double, N>::FMT(q_vec, transform_type_t::Mobius, order_relation, operation);
			const powerset_btree<double, N>& m_back = computation_scheme<double, N>::FMT(q, transform_type_t::Mobius, order_relation, operation);
			t = clock() - t;
			line += "," + std::to_string(((double)t)/CLOCKS_PER_SEC);

			size_t count_z = 0;
			size_t count_m = 0;
			size_t avg_semilattice_support_size = 0;
			for (size_t i = 0; i < nb_of_tests; ++i){
				one_step_EMT(
						order_relation,
						operation,
						scheme_type,
						m,
						build_persistence,
						t,
						count_z,
						count_m,
						avg_semilattice_support_size
				);
				std::fill(m_vec.begin(), m_vec.end(), 0);
				m.clear();
				generate_random_mass_function(m_vec, m, mass_family, support_size);
			}
			line += "," + std::to_string(((double)count_z / nb_of_tests)/CLOCKS_PER_SEC);
			line += "," + std::to_string(((double)count_m / nb_of_tests)/CLOCKS_PER_SEC);
			line += "," + std::to_string((double) avg_semilattice_support_size / nb_of_tests);
			out.emplace_back(line);
			std::cout << "--- Done." << std::endl;
		}

		static void run_without_FMT(order_relation_t order_relation, operation_t operation, scheme_type_t scheme_type, size_t support_size, mass_family_t mass_family, bool build_persistence, size_t nb_of_tests, std::vector<std::string>& out) {
			std::cout << "--- Benchmark with FOD of size " << N << std::endl;
			clock_t t;
			srand (time(NULL));

			std::string sub_fod_labels[N];
			for (size_t i = 0; i < N; ++i){
				sub_fod_labels[i] = std::to_string(i);
			}
			FOD<N> fod(sub_fod_labels);

			std::string line = std::to_string(N);
			line += ",-1,-1";

			size_t count_z = 0;
			size_t count_m = 0;
			size_t avg_semilattice_support_size = 0;
			mass<double, N> m(fod);
			for (size_t i = 0; i < nb_of_tests; ++i){
				generate_random_mass_function(m, mass_family, support_size);
				//m.get_definition().print(std::cout);
				one_step_EMT(
						order_relation,
						operation,
						scheme_type,
						m,
						build_persistence,
						t,
						count_z,
						count_m,
						avg_semilattice_support_size
				);
				m.clear();
			}
			line += "," + std::to_string(((double)count_z / nb_of_tests)/CLOCKS_PER_SEC);
			line += "," + std::to_string(((double)count_m / nb_of_tests)/CLOCKS_PER_SEC);
			line += "," + std::to_string((double) avg_semilattice_support_size / nb_of_tests);
			out.emplace_back(line);
			std::cout << "--- Done." << std::endl;
		}
	};

	template<size_t N>
	const std::string benchmark<N>::alphabet[] = {"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"};


	static void benchmark_main(){
		using namespace efficient_DST;

		constexpr static size_t fod_sizes[] = {15, 18, 20, 22, 23, 24, 25, 26};
		//constexpr static size_t fod_sizes[] = {16, 17, 20, 22, 24, 25, 26};
		//constexpr static size_t fod_sizes[] = {50, 100, 200, 400, 800, 1600};
		mass_family_t mass_family = mass_family_t::random;
		order_relation_t order_relation = order_relation_t::superset;
		bool test_FMT = true;
		size_t support_size = 50;
		bool build_persistence = true;
		size_t nb_of_tests = 1;
		std::vector<scheme_type_t> scheme_types({scheme_type_t::semilattice, scheme_type_t::direct});
		//scheme_type_t scheme_types[] = {scheme_type_t::semilattice, scheme_type_t::direct, scheme_type_t::lattice};

		std::string scheme_type_labels[3] = {"semilattice", "direct", "lattice"};
		std::string family_label;
		if (mass_family == mass_family_t::random){
			family_label = "random";
		}else if(mass_family == mass_family_t::almost_bayesian){
			family_label = "almost_bayesian";
		}else if(mass_family == mass_family_t::almost_consonant){
			family_label = "almost_consonant";
		}
		std::string order_label;
		if (order_relation == order_relation_t::subset){
			order_label = "subset";
		}else{
			order_label = "superset";
		}
		operation_t operation = operation_t::addition;

		for (size_t sc = 0; sc < scheme_types.size(); ++sc){
			std::string filename = "benchmark_scheme-"
									+ scheme_type_labels[sc]
									+ "_persistence-" + std::to_string(build_persistence)
									+ "_order-" + order_label
									+ "_support-" + std::to_string(support_size)
									+ "_family-" + family_label
									+ "_fodsize-" + std::to_string(fod_sizes[0]) + "_to_" + std::to_string(fod_sizes[5])
			+ ".csv";
			//if(std::remove(filename.c_str()) == 0){
			//	puts("Previous data file of this experiment deleted.");
			//}
			std::cout << "Benchmark with scheme type " << scheme_type_labels[sc] << std::endl;
			std::vector<std::string> out;
			out.reserve(6);
			if(test_FMT){
				benchmark<fod_sizes[0]>::run(order_relation, operation, scheme_types[sc], support_size, mass_family, build_persistence, nb_of_tests, out);
				//benchmark<fod_sizes[1]>::run(order_relation, operation, scheme_types[sc], support_size, mass_family, build_persistence, nb_of_tests, out);
				//benchmark<fod_sizes[2]>::run(order_relation, operation, scheme_types[sc], support_size, mass_family, build_persistence, nb_of_tests, out);
				//benchmark<fod_sizes[3]>::run(order_relation, operation, scheme_types[sc], support_size, mass_family, build_persistence, nb_of_tests, out);
				//benchmark<fod_sizes[4]>::run(order_relation, operation, scheme_types[sc], support_size, mass_family, build_persistence, nb_of_tests, out);
				//benchmark<fod_sizes[5]>::run(order_relation, operation, scheme_types[sc], support_size, mass_family, build_persistence, nb_of_tests, out);
			}else{
				benchmark<fod_sizes[0]>::run_without_FMT(order_relation, operation, scheme_types[sc], support_size, mass_family, build_persistence, nb_of_tests, out);
				benchmark<fod_sizes[1]>::run_without_FMT(order_relation, operation, scheme_types[sc], support_size, mass_family, build_persistence, nb_of_tests, out);
				benchmark<fod_sizes[2]>::run_without_FMT(order_relation, operation, scheme_types[sc], support_size, mass_family, build_persistence, nb_of_tests, out);
				benchmark<fod_sizes[3]>::run_without_FMT(order_relation, operation, scheme_types[sc], support_size, mass_family, build_persistence, nb_of_tests, out);
				benchmark<fod_sizes[4]>::run_without_FMT(order_relation, operation, scheme_types[sc], support_size, mass_family, build_persistence, nb_of_tests, out);
				//benchmark<fod_sizes[5]>::run_without_FMT(order_relation, operation, scheme_types[sc], support_size, mass_family, build_persistence, nb_of_tests, out);
			}

			std::ofstream myfile;
			myfile.open(filename);
			for (unsigned int i = 0; i < out.size(); ++i){
				myfile << out[i] << std::endl;
				std::cout << out[i] << std::endl;
			}
			myfile.close();
			std::cout << std::endl;
		}
	}
}
