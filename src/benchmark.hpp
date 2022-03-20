#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <random>


#include <commonality_function.hpp>
#include <commonality_vector.hpp>
#include <conjunctive_decomposition.hpp>
#include <disjunctive_decomposition.hpp>
#include <implicability_function.hpp>
#include <implicability_vector.hpp>
#include <mass_function.hpp>


namespace efficient_DST{

	enum class mass_family_t: int8_t { random, almost_consonant, almost_bayesian };
	enum class order_relation_t: bool { subset, superset };
	enum class transform_type_t: bool { zeta, Mobius };

	template <size_t N, typename T = float>
	class benchmark {
	public:

		static inline std::bitset<N> random_bitset( double p = 1/N, size_t max_card = N) {

			std::bitset<N> bits;
			std::random_device rd;
			std::mt19937 gen( rd());
			std::bernoulli_distribution d((double) 1/N);

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

		template<typename subset, class mass_func>
		static void generate_random_mass_function(
				mass_func& m,
//				mass_function<N, T>& m,
				const mass_family_t& mass_family,
				const float& proportion
		){
			size_t size;
			T mass;
			std::random_device rd;
			std::mt19937 gen(rd());
			std::bernoulli_distribution d(proportion);
			size_t n = 0;
			size_t special_size;
			switch (mass_family) {
				case mass_family_t::almost_consonant:
				{
					size = N+1;
					special_size = (size_t) (proportion * size);
					mass = 1.0/size;
					subset set(0);
//					std::bitset<N> set(0);
					while (n < special_size){
						if (m[0] == 0 && d(gen)){
							m.assign_emptyset(mass);
							++n;
						}
						subset singleton(1);
//						std::bitset<N> singleton(1);
						for (size_t i = 0; i < N; ++i){
							if (n == special_size)
								break;
							if ((set & singleton) == 0 && d(gen)){
								set |= singleton;
								m.assign(set, mass);
								++n;
							}
							singleton <<= 1;
						}
					}
				}break;
				case mass_family_t::almost_bayesian:
				{
					size = N;
					special_size = (size_t) (proportion * size);
					mass = 1.0/size;
					while (n < special_size){
						subset singleton(1);
//						std::bitset<N> singleton(1);
						for (size_t i = 0; i < N; ++i){
							if (n == special_size)
								break;
							if ((m[singleton] == 0) && d(gen)){
								m.assign(singleton, mass);
								++n;
							}
							singleton <<= 1;
						}
					}
				}break;
				default:
				{
					size = (size_t) (proportion * (1 << N));
					mass = 1.0/size;
				}break;
			}
			if (mass == 0){
				std::cerr << "Number of elements in support exceeding machine precision\n";
				exit(1);
			}
			while (n < size){
				const subset& random_set = random_bitset();
//				const std::bitset<N>& random_set = random_bitset();
				if (m[random_set] == 0) {
					m.assign(random_set, mass);
					++n;
				}
			}
		}

		static void run_transformations(
				const mass_function<N, T>& m,
				const scheme_type_t& scheme_type,
				const order_relation_t& order_relation,
				clock_t& t,
				size_t& t_zeta,
				size_t& t_mobius
		){
			if (order_relation == order_relation_t::superset){
				t = clock();
				commonality_function<N, T> q(m, scheme_type);
				t = clock() - t;
				t_zeta += t;
				t = clock();
				weight_function<N, T> w(q);
				t = clock() - t;
				t_mobius += t;
			}else{
				t = clock();
				implicability_function<N, T> b(m, scheme_type);
				t = clock() - t;
				t_zeta += t;
				t = clock();
				weight_function<N, T> v(b);
				t = clock() - t;
				t_mobius += t;
			}
		}

		static void run_vector_transformations(
				const mass_vector<N, T>& m,
				const bool& core_reduced,
				const order_relation_t& order_relation,
				clock_t& t,
				size_t& t_zeta,
				size_t& t_mobius
		){
			if (order_relation == order_relation_t::superset){
				t = clock();
				commonality_vector<N, T> q(m, core_reduced);
				t = clock() - t;
				t_zeta += t;
				t = clock();
				weight_vector<N, T> w(q);
				t = clock() - t;
				t_mobius += t;
			}else{
				t = clock();
				implicability_vector<N, T> b(m, core_reduced);
				t = clock() - t;
				t_zeta += t;
				t = clock();
				weight_vector<N, T> v(b);
				t = clock() - t;
				t_mobius += t;
			}
		}

		static void mass_function_to_vector(
				const mass_function<N, T>& m,
				mass_vector<N, T>& m_vec
		) {
			const std::vector<set_N_value<N, T>* >& elements = m.get_definition().elements();
			for (size_t i = 0; i < elements.size(); ++i){
				m_vec.assign(elements[i]->set.to_ullong(), elements[i]->value);
			}
		}

		static void run(){
			clock_t t;
			std::vector<mass_family_t> mass_families = {mass_family_t::almost_consonant, mass_family_t::almost_bayesian, mass_family_t::random};
			std::vector<float> proportions = {0.2, 0.4, 0.6, 0.8, 1};
			std::vector<order_relation_t> order_relations = {order_relation_t::superset, order_relation_t::subset};
			std::vector<scheme_type_t> scheme_types = {scheme_type_t::direct, scheme_type_t::semilattice, scheme_type_t::lattice};
			size_t nb_of_tests = 5;

			std::vector<std::string> scheme_type_labels = {"direct", "EMT_semilattice", "EMT_lattice", "reduced_FMT"};
			std::vector<std::string> family_labels = {"almost_consonant", "almost_bayesian", "random"};
			std::vector<std::string> order_labels = {"superset", "subset"};

			std::string outcome_labels[N];
			for (size_t i = 0; i < N; ++i){
				outcome_labels[i] = std::to_string(i);
			}
			sample_space<N> outcomes(outcome_labels);
			std::cout << "--- Benchmark with " << N << " possible outcomes.\n";
			for (size_t mf = 0; mf < mass_families.size(); ++mf){
				std::cout << "Mass family " << family_labels[mf] << std::endl;
				for (size_t p = 0; p < proportions.size(); ++p){
					std::cout << "Proportion " << proportions[p] << std::endl;
					for (size_t it = 0; it < nb_of_tests; ++it){
						std::cout << "Iteration " << it << std::endl;
						mass_function<N, T> m(outcomes);
						generate_random_mass_function<std::bitset<N>, mass_function<N, T>>(m, mass_families[mf], proportions[p]);
						for (size_t sc = 0; sc <= scheme_types.size(); ++sc){
							std::cout << "Scheme type " << scheme_type_labels[sc] << std::endl;
							if (sc < scheme_types.size()){
								for (size_t o = 0; o < order_relations.size(); ++o){
									size_t t_zeta = 0;
									size_t t_mobius = 0;
									run_transformations(
											m,
											scheme_types[sc],
											order_relations[o],
											t,
											t_zeta,
											t_mobius
									);
									std::string filename = "benchmark_N-" + std::to_string(N)
										+ "_family-" + family_labels[mf]
										+ "_prop-" + std::to_string(proportions[p])
										+ "_scheme-" + scheme_type_labels[sc]
										+ "_order-" + order_labels[o]
									+ ".csv";
									std::ofstream myfile;
									myfile.open(filename);
									//if(std::remove(filename.c_str()) == 0){
									//	puts("Previous data file of this experiment deleted.");
									//}
									std::string line = std::to_string((double) t_zeta / CLOCKS_PER_SEC);
									line += "," + std::to_string((double) t_mobius / CLOCKS_PER_SEC);
									myfile << line << std::endl;
									myfile.close();
								}
							}else{
								mass_vector<N, T> m_vec(outcomes);
								benchmark<N, T>::mass_function_to_vector(m, m_vec);
								for (size_t o = 0; o < order_relations.size(); ++o){
									size_t t_zeta = 0;
									size_t t_mobius = 0;
									run_vector_transformations(
											m_vec,
											true,
											order_relations[o],
											t,
											t_zeta,
											t_mobius
									);
									std::string filename = "benchmark/benchmark_N-" + std::to_string(N)
										+ "_family-" + family_labels[mf]
										+ "_prop-" + std::to_string(proportions[p])
										+ "_scheme-" + scheme_type_labels[sc]
										+ "_order-" + order_labels[o]
									+ ".csv";
									std::ofstream myfile;
									myfile.open(filename);
									//if(std::remove(filename.c_str()) == 0){
									//	puts("Previous data file of this experiment deleted.");
									//}
									std::string line = std::to_string((double) t_zeta / CLOCKS_PER_SEC);
									line += "," + std::to_string((double) t_mobius / CLOCKS_PER_SEC);
									myfile << line << std::endl;
									myfile.close();
								}
							}
						}
					}
				}
			}
		}
	};


	static void benchmark_main(){
		using namespace efficient_DST;

		typedef float T;
		//constexpr static size_t fod_sizes[] = {16, 17, 20, 22, 24, 25, 26};
		//constexpr static size_t fod_sizes[] = {50, 100, 200, 400, 800, 1600};
		srand(time(NULL));

//		for (size_t n = 0; n < nb_N; ++n){
		benchmark<15, T>::run();
		benchmark<18, T>::run();
		benchmark<21, T>::run();
		benchmark<24, T>::run();
		benchmark<27, T>::run();
		benchmark<30, T>::run();
		std::cout << "--- Done." << std::endl;
	}
}
