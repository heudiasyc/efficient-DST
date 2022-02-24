#ifndef EFFICIENT_DST_MASS_HPP
#define EFFICIENT_DST_MASS_HPP

#include <mobius_transform.hpp>
#include <zeta_transform.hpp>
#include <conjunctive_decomposition.hpp>


namespace efficient_DST{

	enum class special_case_t: bool { degenerate, vacuous };

	template <size_t N, typename T = float>
	class mass : public mobius_transform<N, T>{
	public:
		using typename powerset_function<N, T>::subset;
		using powerset_function<N, T>::emptyset;
		using powerset_function<N, T>::fullset;

		mass(const mass<N, T>& m) : mobius_transform<N, T>(m.outcomes, m.definition, 0)
		{}

		mass(
			const sample_space<N>& outcomes,
			const powerset_btree<N, T>& focal_sets
		) : mobius_transform<N, T>(outcomes, focal_sets, 0)
		{
			this->remove_negligible_values();
			this->normalize();
		}

		mass(const sample_space<N>& outcomes) : mobius_transform<N, T>(outcomes, 0)
		{}

		mass(const sample_space<N>& outcomes, const special_case_t s_case) : mobius_transform<N, T>(outcomes, 0)
		{
			switch(s_case){
				// create a mass function with all mass attributed to the empty set
			    case special_case_t::degenerate  : this->set_emptyset_value(1);	break;
				// create a mass function with all mass attributed to the FOD set
			    case special_case_t::vacuous	:  this->set_fullset_value(1);	break;
			}
		}

		mass(const zeta_transform<up_inclusion<N, T>, N, T >& q) : mass<N, T>(q.get_sample_space(), q.inversion(operation_type_t::addition))
		{}

		mass(const zeta_transform<down_inclusion<N, T>, N, T >& b) : mass<N, T>(b.get_sample_space(), b.inversion(operation_type_t::addition))
		{}

		mass(const conjunctive_decomposition<N, T >& w_dec) : mass<N, T>(w_dec.get_sample_space())
		{
//			w_dec.get_definition().print(this->outcomes);
			fuse_decomposition<up_inclusion<N, T> >(w_dec, *this);
//			this->print(true);
			this->remove_negligible_values();
			normalize();
		}


		template<class inclusion>
		mass<N, T> fuse_decomposition(const decomposition<inclusion, N, T>& w_dec){
			mass<N, T> m1(this->outcomes);
			fuse_decomposition(w_dec, m1);
			this->remove_negligible_values();
			normalize();
			return m1;
		}

		template<class inclusion>
		void fuse_decomposition(const decomposition<inclusion, N, T>& w_dec, mass<N, T>& m1){
//			std::cout << "Decomposition to fuse:\n";
//			w_dec.print();
			const powerset_btree<N, T>& inverse_weights = w_dec.get_definition();
			const std::vector<set_N_value<N, T>* >& elements = inverse_weights.elements();
			mass<N, T> m12(this->outcomes);
			size_t i = 0;
			if (elements[i]->set == w_dec.normalizing_set_assignment.set){
				++i;
			}
			T weight = 1/elements[i]->value;
			m1.assign(w_dec.normalizing_set_assignment.set, weight);
			m1.assign(elements[i]->set, 1-weight);
			for (++i; i < elements.size(); ++i){
				if (elements[i]->set != w_dec.normalizing_set_assignment.set){
					mass<N, T> m2(this->outcomes);
					weight = 1/elements[i]->value;
					m2.assign(w_dec.normalizing_set_assignment.set, weight);
					m2.assign(elements[i]->set, 1-weight);
					m1.natural_fusion_with<inclusion>(m2, m12);
					m1.clear();
					m1.definition.copy(m12.definition);
					m12.clear();
				}
			}
		}

		template<class inclusion>
		mass<N, T> natural_fusion_with(const mass<N, T>& m2) const {
			mass<N, T> m12(this->outcomes);
			natural_fusion_with<inclusion>(m2, m12);
			m12.remove_negligible_values();
			m12.normalize();
			return m12;
		}

		template<class inclusion>
		void natural_fusion_with(const mass<N, T>& m2, mass<N, T>& m12) const {
			const std::vector<set_N_value<N, T>* >& focal_sets_1 = this->definition.elements();
			const std::vector<set_N_value<N, T>* >& focal_sets_2 = m2.definition.elements();
			powerset_btree<N, T>& focal_sets_12 = m12.definition;

			for (size_t i1 = 0; i1 < focal_sets_1.size(); ++i1){
				for (size_t i2 = 0; i2 < focal_sets_2.size(); ++i2){
					const subset& set = inclusion::set_operation(focal_sets_1[i1]->set, focal_sets_2[i2]->set);
					set_N_value<N, T>* node = focal_sets_12[set];
					if (node){
						node->value += focal_sets_1[i1]->value * focal_sets_2[i2]->value;
					}else{
						focal_sets_12.insert(set, focal_sets_1[i1]->value * focal_sets_2[i2]->value);
					}
				}
			}
//			mass<N, T> m12(m1.get_sample_space(), focal_sets_12);
//			m12.remove_negligible_values();
//			m12.normalize();
		}


		template <class fusion_rule>
		mass<N, T> fuse_with(const mass<N, T>& m2) const {
			const fusion_rule fusion;
			return fusion(*this, m2);
		}

		void regularize() {
			this->definition.nullify(this->definition[emptyset]);
			this->normalize();
		}

		void normalize() {
			normalize(this->definition);
		}

		static void normalize(powerset_btree<N, T>& definition) {
			T sum = 0;
			const std::vector<set_N_value<N, T>* >& elements = definition.elements();
			for (size_t i = 0; i < elements.size(); ++i) {
				sum += elements[i]->value;
			}
			if(sum == 0){
				std::cerr << "\nSum of mass values equal to 0."
						<< "\nThis means that this mass function is either empty or contains as much positive values as negative values."
						<< "\nEither way, this mass function cannot be normalized into a valid mass function." << std::endl;;
				exit(1);
			}
			if(sum != 1){
				// normalize
				for (size_t i = 0; i < elements.size(); ++i) {
					elements[i]->value /= sum;
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
			const std::vector<set_N_value<N, T>* >& f_elements = this->definition.elements();
			T sum = 0;
			for (size_t i = 0; i < f_elements.size(); ++i) {
				if(f_elements[i]->value < 0 && !this->is_equivalent_to_zero(f_elements[i]->value))
					return false;
				sum += f_elements[i]->value;
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
			return this->definition.elements().size() == 1;
		}

		/// Simple mass function has at most two focal sets, and if it has two,
		/// Omega is one of them.
		bool is_simple() const {
			const std::vector<set_N_value<N, T>* >& elements = this->definition.elements();

			if(elements.size() == 2){
				return this->definition.full_set() != nullptr;
			}else if(elements.size() < 2){
				return true;
			}else{
				return false;
			}
		}

		/// In bayesian mass function, all focal sets are singletons.
		bool is_bayesian() const {
			T sum = 0;
			const std::vector<set_N_value<N, T>* >& singletons = this->definition.singletons();

			for (size_t i = 0; i < singletons.size(); ++i) {
				sum += singletons[i]->value;
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

			const std::vector<set_N_value<N, T>* >& elements = this->definition.elements();
			const subset& U = elements[0]->set;

			for (size_t i = 1; i < elements.size(); ++i) {
				const subset& set = elements[i]->set;
				if ((set & U) != 0) {
					return false;
				}
				U |= set;
			}
			return true;
		}

		/// Focal sets of a partitioned mass function do not intersect and their
		/// union is equal to the frame of discernment (FOD).
		/// In other words, the union of all focal elements is the FOD,
		/// and the intersection of each focal set with the union of other focal sets is empty.
		bool is_partitioned() const {
			const std::vector<set_N_value<N, T>* >& elements = this->definition.elements();
			const subset& U = elements[0]->set;

			for (size_t i = 1; i < elements.size(); ++i) {
				const subset& set = elements[i]->set;
				if ((set & U) != 0) {
					return false;
				}
				U |= set;
			}
			return U == fullset;
		}

		/// Mass function without internal conflict is a one where all pairs of
		/// focal elements have a non-empty intersection.
		/// \f[
		/// \forall A, B \subseteq \Omega, m(A) > 0, m(B) > 0 : A \cap B \neq \emptyset
		/// \f]
		/// Therefore, there is internal conflict when there is at least one empty intersection between pairs of focal sets,
		/// i.e. when the intersection of all focal elements is empty.
		bool has_internal_conflict() const {
			const std::vector<set_N_value<N, T>* >& elements = this->definition.elements();
			const subset& I = elements[0]->set;

			for (size_t i = 1; i < elements.size(); ++i) {
				const subset& set = elements[i]->set;
				I &= set;
				if (I == emptyset){
					return true;
				}
			}
			return false;
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_MASS_HPP
