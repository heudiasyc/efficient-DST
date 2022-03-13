#ifndef EFFICIENT_DST_MASS_VECTOR_HPP
#define EFFICIENT_DST_MASS_VECTOR_HPP

#include <mobius_inversion.hpp>
#include <mobius_transform_vector.hpp>
#include <zeta_transform_vector.hpp>
#include <conjunctive_decomposition_vector.hpp>
#include <disjunctive_decomposition_vector.hpp>


namespace efficient_DST{

	enum class special_case_t: bool { degenerate, vacuous };

	template <size_t N, typename T = float>
	class mass_vector : public mobius_transform_vector<N, T>{
	public:

		mass_vector(const mass_vector<N, T>& m) : mobius_transform_vector<N, T>(m.outcomes, m.definition, 0)
		{}

		mass_vector(
			const sample_space<N>& outcomes,
			const std::vector<T>& focal_sets
		) : mobius_transform_vector<N, T>(outcomes, focal_sets, 0)
		{
			this->remove_negligible_values();
			this->normalize();
		}

		mass_vector(const sample_space<N>& outcomes) : mobius_transform_vector<N, T>(outcomes, 0)
		{}

		mass_vector(const sample_space<N>& outcomes, const special_case_t s_case) : mobius_transform_vector<N, T>(outcomes, 0)
		{
			switch(s_case){
				// create a mass function with all mass attributed to the empty set
			    case special_case_t::degenerate  : this->assign_emptyset(1);	break;
				// create a mass function with all mass attributed to the FOD set
			    case special_case_t::vacuous	:  this->assign_fullset(1);	break;
			}
		}

		mass_vector(const zeta_transform_vector<up_inclusion<N, T>, N, T >& q) : mass_vector<N, T>(q.get_sample_space(), q.inversion(operation_type_t::addition))
		{}

		mass_vector(const zeta_transform_vector<down_inclusion<N, T>, N, T >& b) : mass_vector<N, T>(b.get_sample_space(), b.inversion(operation_type_t::addition))
		{}

		mass_vector(const conjunctive_decomposition_vector<N, T >& w_dec) : mass_vector<N, T>(w_dec.get_sample_space())
		{
//			w_dec.get_definition().print(this->outcomes);
			fuse_decomposition<up_inclusion<N, T> >(w_dec, *this);
//			this->print(true);
			this->remove_negligible_values();
			normalize();
		}

		mass_vector(const disjunctive_decomposition_vector<N, T >& w_dec) : mass_vector<N, T>(w_dec.get_sample_space())
		{
//			w_dec.get_definition().print(this->outcomes);
			fuse_decomposition<down_inclusion<N, T> >(w_dec, *this);
//			this->print(true);
			this->remove_negligible_values();
			normalize();
		}


		template<class inclusion>
		mass_vector<N, T> fuse_decomposition(const decomposition_vector<inclusion, N, T>& w_dec){
			mass_vector<N, T> m1(this->outcomes);
			fuse_decomposition(w_dec, m1);
			this->remove_negligible_values();
			normalize();
			return m1;
		}

		template<class inclusion>
		void fuse_decomposition(const decomposition_vector<inclusion, N, T>& w_dec, mass_vector<N, T>& m1){
//			std::cout << "Decomposition to fuse:\n";
//			w_dec.print();
			const std::vector<T>& inverse_weights = w_dec.get_definition();
			mass_vector<N, T> m12(this->outcomes);
			m1.assign(w_dec.normalizing_set, 1);
			size_t j = 0;
			for (size_t i = 0; i < (1 << N); ++i){
				if (inclusion::set_operation(i, w_dec.normalizing_set) == i){
					if (i != w_dec.normalizing_set && inverse_weights[j] != 1){
						mass_vector<N, T> m2(this->outcomes);
						T weight = 1/inverse_weights[j];
						m2.assign(w_dec.normalizing_set, weight);
						m2.assign(i, 1-weight);
						m1.natural_fusion_with<inclusion>(m2, m12);
						m1.definition = m12.definition;
						m12.clear();
					}
					++j;
				}
			}
		}

		template<class inclusion>
		mass_vector<N, T> natural_fusion_with(const mass_vector<N, T>& m2) const {
			mass_vector<N, T> m12(this->outcomes);
			natural_fusion_with<inclusion>(m2, m12);
			m12.remove_negligible_values();
			m12.normalize();
			return m12;
		}

		template<class inclusion>
		void natural_fusion_with(const mass_vector<N, T>& m2, mass_vector<N, T>& m12) const {
//			std::cout << "Masses to fuse:\n";
//			this->print();
//			m2.print();
			for (size_t i1 = 0; i1 < this->definition.size(); ++i1){
				for (size_t i2 = 0; i2 < m2.definition.size(); ++i2){
					const size_t& set = inclusion::set_operation(i1, i2);
					m12.definition[set] += this->definition[i1] * m2[i2];
				}
			}
//			mass<N, T> m12(m1.get_sample_space(), focal_sets_12);
//			m12.remove_negligible_values();
//			m12.normalize();
		}


		template <class fusion_rule>
		mass_vector<N, T> fuse_with(const mass_vector<N, T>& m2) const {
			const fusion_rule fusion;
			return fusion(*this, m2);
		}

		void regularize() {
			this->definition[0] = 0;
			this->normalize();
		}

		void normalize() {
			normalize(this->definition);
		}

		static void normalize(std::vector<T>& definition) {
			T sum = 0;
			for (size_t i = 0; i < definition.size(); ++i) {
				sum += definition[i];
			}
			if(sum == 0){
				std::cerr << "\nSum of mass values equal to 0."
						<< "\nThis means that this mass function is either empty or contains as much positive values as negative values."
						<< "\nEither way, this mass function cannot be normalized into a valid mass function." << std::endl;;
				exit(1);
			}
			if(sum != 1){
				// normalize
				for (size_t i = 0; i < definition.size(); ++i) {
					definition[i] /= sum;
				}
			}
		}

//		void remove_negligible_values() {
//			remove_negligible_values(this->definition);
//		}
//
//		static void remove_negligible_values(powerset_btree<N, T>& definition) {
//			mobius_transform<N, T>::remove_negligible_values(definition, this->default_value);
//		}


		///////////////////////////////////////////////////////////////////


		/// Check if this mass function is valid (i.e. if the sum of all its images is 1)
		bool is_valid() const {
			T sum = 0;
			for (size_t i = 0; i < this->definition.size(); ++i) {
				if(this->definition[i] < 0 && !this->is_equivalent_to_zero(this->definition[i]))
					return false;
				sum += this->definition[i];
			}

			return this->is_equivalent_to_zero(1-sum);
		}

		/// Normal mass function has no mass at empty set (Conflict).
		bool is_normal() const {
			return this->is_equivalent_to_zero(this->at_emptyset());
		}

		/// Subnormal mass function has some mass at empty set (Conflict).
		bool is_subnormal() const {
			return !is_normal();
		}

		/// Dogmatic mass function has no mass at FOD.
		bool is_dogmatic() const {
			return this->is_equivalent_to_zero(this->at_fullset());
		}

		/// Non-dogmatic mass function has some mass at FOD.
		bool is_nondogmatic() const {
			return !is_dogmatic();
		}

		/// Vacuous mass function has all mass at FOD.
		bool is_vacuous() const {
			return this->is_equivalent_to_zero(1 - this->at_fullset());
		}

		/// Degenerate mass function has all mass at empty set (Conflict).
		bool is_degenerate() const {
			return this->is_equivalent_to_zero(1 - this->at_emptyset());
		}

		/// Categorical mass function has only one focal set.
		bool is_categorical() const {
			return this->definition.size() == 1;
		}

		/// Simple mass function has at most two focal sets, and if it has two,
		/// Omega is one of them.
		bool is_simple() const {
			if(this->definition.size() == 2){
				return this->at_fullset() != 0;
			}else if(this->definition.size() < 2){
				return true;
			}else{
				return false;
			}
		}

		/// In bayesian mass function, all focal sets are singletons.
		bool is_bayesian() const {
			if (this->definition.size() == 0)
				return false;
			T sum = 0;
			size_t singleton = 1;
			for (size_t i = 0; i < (size_t) log2(this->definition.size()); ++i) {
				sum += this->definition[singleton];
				singleton << 1;
			}
			return this->is_equivalent_to_zero(1-sum);
		}

		/// Mass function is consonant when all focal sets are nested.
		///
		/// Consonant evidence can be represented as a nested structure of subsets
		/// where the
		/// elements of the smallest set are included in the next larger set, all of
		/// whose elements are
		/// included in the next larger set and so on.
		bool is_consonant() const {
			return consonance_check<N, T>(this->definition);
		}

		/// Disjoint evidence implies that any two focal sets have no element in common.
		/// In other words, the intersection of each focal set with the union of other focal sets is empty.
		bool is_disjoint() const {
			size_t U = 0;
			for (size_t i = 0; i < this->definition.size(); ++i) {
				if ((i & U) != 0) {
					return false;
				}
				U |= i;
			}
			return true;
		}

		/// Focal sets of a partitioned mass function do not intersect and their
		/// union is equal to the frame of discernment (FOD).
		/// In other words, the union of all focal elements is the FOD,
		/// and the intersection of each focal set with the union of other focal sets is empty.
		bool is_partitioned() const {
			size_t U = 0;
			for (size_t i = 0; i < this->definition.size(); ++i) {
				if ((i & U) != 0) {
					return false;
				}
				U |= i;
			}
			return U == N-1;
		}

		/// Mass function without internal conflict is a one where all pairs of
		/// focal elements have a non-empty intersection.
		/// \f[
		/// \forall A, B \subseteq \Omega, m(A) > 0, m(B) > 0 : A \cap B \neq \emptyset
		/// \f]
		/// Therefore, there is internal conflict when there is at least one empty intersection between pairs of focal sets,
		/// i.e. when the intersection of all focal elements is empty.
		bool has_internal_conflict() const {
			size_t I = 0;
			I = ~I;
			for (size_t i = 0; i < this->definition.size(); ++i) {
				I &= i;
				if (I == 0){
					return true;
				}
			}
			return false;
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_MASS_VECTOR_HPP
