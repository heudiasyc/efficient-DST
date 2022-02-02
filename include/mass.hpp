#ifndef EFFICIENT_DST_MASS_HPP
#define EFFICIENT_DST_MASS_HPP

#include <mobius_transform.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	enum class special_case_t: bool { degenerate, vacuous };

	template <typename T, size_t N>
	class mass : public mobius_transform<T, N>{
	public:

		mass(const mass<T, N>& m) : mobius_transform<T, N>(m.get_definition())
		{}

		mass(const powerset_btree<T, N>& support) : mobius_transform<T, N>(support)
		{
			this->remove_negligible_values();
			this->normalize();
		}

		mass(FOD<N>& fod) : mobius_transform<T, N>(fod)
		{}

		mass(FOD<N>& fod, const special_case_t s_case) : mobius_transform<T, N>(fod)
		{
			switch(s_case){
				// create a mass function with all mass attributed to the empty set
			    case special_case_t::degenerate  : this->set_emptyset_value(1);	break;
				// create a mass function with all mass attributed to the FOD set
			    case special_case_t::vacuous	:  this->set_fod_value(1);	break;
			}
		}

		mass(const zeta_transform<T, N, up_inclusion<T, N> >& q) : mass<T, N>(q.inversion(operation_type_t::addition))
		{}

		mass(const zeta_transform<T, N, down_inclusion<T, N> >& b) : mass<T, N>(b.inversion(operation_type_t::addition))
		{}


		template <class fusion_rule>
		mass<T, N> apply(const mass<T, N>& m2) const {
			const fusion_rule fusion;
			return fusion(*this, m2);
		}

		T at_emptyset() const {
			return powerset_function<T, N>::at_emptyset(0);
		}

		T at_fod() const {
			return powerset_function<T, N>::at_fod(0);
		}

		T operator[](const std::vector<std::string>& labels) const {
			return find(this->definition.get_FOD()->to_set(labels));
		}

		T find(const std::bitset<N>& set) const {
			return powerset_function<T, N>::find(set, 0);
		}

		void regularize() {
			this->nullify({});
			this->normalize();
		}

		void normalize() {
			normalize(this->definition);
		}

		static void normalize(powerset_btree<T, N>& definition) {
			T sum = 0;
			const std::vector<set_N_value<T, N>* >& elements = definition.elements();
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

		void remove_negligible_values() {
			remove_negligible_values(this->definition);
		}

		static void remove_negligible_values(powerset_btree<T, N>& definition) {
			mobius_transform<T, N>::remove_negligible_values(definition, 0);
		}


		///////////////////////////////////////////////////////////////////


		/// Check if this mass function is valid (i.e. if the sum of all its images is 1)
		bool is_valid() const {
			const std::vector<set_N_value<T, N>* >& f_elements = this->definition.elements();
			T sum = 0;
			for (size_t i = 0; i < f_elements.size(); ++i) {
				if(f_elements[i]->value < 0 && !this->is_equivalent_to_zero(f_elements[i]->value))
					return false;
				sum += f_elements[i]->value;
			}

			return this->is_equivalent_to_zero(1-sum);
		}

		/// A focal set has a non-zero mass attributed.
		bool is_focal(const std::bitset<N>& set) const {
			return !this->is_equivalent_to_zero(find(set));
		}

		bool is_focal(const std::vector<std::string>& labels) const {
			return !this->is_equivalent_to_zero(this->operator[](labels));
		}

		/// Normal mass function has no mass at empty set (Conflict).
		bool is_normal() const {
			return this->is_equivalent_to_zero(this->at_emptyset());
		}

		/// Regular mass function has no mass at empty set (Conflict).
		bool is_regular() const {
			return is_normal();
		}

		/// Subnormal mass function has some mass at empty set (Conflict).
		bool is_subnormal() const {
			return !is_normal();
		}

		/// Dogmatic mass function has no mass at FOD.
		bool is_dogmatic() const {
			return this->is_equivalent_to_zero(this->at_fod());
		}

		/// Non-dogmatic mass function has some mass at FOD.
		bool is_nondogmatic() const {
			return !is_dogmatic();
		}

		/// Vacuous mass function has all mass at FOD.
		bool is_vacuous() const {
			return this->is_equivalent_to_zero(1 - this->at_fod());
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
			const std::vector<set_N_value<T, N>* >& elements = this->definition.elements();

			if(elements.size() == 2){
				return at_fod() != nullptr;
			}else if(elements.size() < 2){
				return true;
			}else{
				return false;
			}
		}

		/// In bayesian mass function, all focal sets are singletons.
		bool is_bayesian() const {
			T belief = 0;
			const std::vector<set_N_value<T, N>* >& singletons = this->definition.singletons();

			for (size_t i = 0; i < singletons.size(); ++i) {
				belief += singletons[i]->value;
			}
			return this->is_equivalent_to_zero(1-belief);
		}

		/// Mass function is consonant when all focal sets are nested.
		///
		/// Consonant evidence can be represented as a nested structure of subsets
		/// where the
		/// elements of the smallest set are included in the next larger set, all of
		/// whose elements are
		/// included in the next larger set and so on.
		bool is_consonant() const {
			std::unordered_map<size_t, std::vector<set_N_value<T, N>* > > elements_by_cardinality = this->definition.elements_by_set_cardinality();

			for(auto kv : elements_by_cardinality) {
				if(kv.second.size() > 1)
					// if there is more than one element of same cardinality, the structure isn't nested
					return false;
			}
			// at this point, the structure has at most one element per cardinality

			// sort indices in ascending order to check if each focal set is contained in all bigger focal sets
			const std::vector<std::vector<set_N_value<T, N>* >* >& ordered_cardinalities = this->definition.get_sorted_cardinalities(elements_by_cardinality);
			size_t c = 0;
			while(c < ordered_cardinalities.size()-1){
				if (!FOD<N>::is_or_is_subset_of(
						elements_by_cardinality[ordered_cardinalities[c]][0]->set,
						elements_by_cardinality[ordered_cardinalities[c+1]][0]->set)) {
					return false;
				}
				++c;
			}

			return true;
		}

		/// Mass function is consistent when there is at least one focal set that is
		/// common to all focal sets.
		bool is_consistent() const {
			std::unordered_map<size_t, std::vector<set_N_value<T, N>* > > elements_by_cardinality = this->definition.elements_by_set_cardinality();
			// sort indices in ascending order to check if the smallest focal set is contained in all bigger focal sets
			const std::vector<std::vector<set_N_value<T, N>* >* >& ordered_cardinalities = this->definition.get_sorted_cardinalities(elements_by_cardinality);

			if(ordered_cardinalities[0]->size() > 1)
				// if there is more than one element of smallest cardinality, it is
				// inconsistent since two sets of same cardinality cannot contain each other
				return false;

			// at this point, there is only one set that has the smallest cardinality

			for (size_t c = 1; c < ordered_cardinalities.size(); ++c) {
				for (size_t i = 0; i < ordered_cardinalities.size(); ++i) {
					if (!FOD<N>::is_or_is_subset_of(
							elements_by_cardinality[ordered_cardinalities[0]][0]->set,
							elements_by_cardinality[ordered_cardinalities[c]][i]->set)) {
						return false;
					}
				}
			}
			return true;
		}

		/// Arbitrary evidence corresponds to the situation where there is no
		/// focal set common to all focal sets, though some focal sets may have elements in common.
		/// This is the opposite of a consistent structure.
		bool is_arbitrary() const {
			return !is_consistent();
		}

		/// Disjoint evidence implies that any two focal sets have no element in common.
		/// In other words, the intersection of each focal set with the union of other focal sets is empty.
		bool is_disjoint() const {

			const std::vector<set_N_value<T, N>* >& elements = this->definition.elements();
			const std::bitset<N>& U = elements[0]->set;

			for (size_t j = 1; j < elements.size(); ++j) {
				const std::bitset<N>& setj = elements[j]->set;
				if (!FOD<N>::are_disjoint(setj, U)) {
					return false;
				}
				U |= setj;
			}
			return true;
		}

		/// Focal sets of a partitioned mass function do not intersect and their
		/// union is equal to the frame of discernment (FOD).
		/// In other words, the union of all focal elements is the FOD,
		/// and the intersection of each focal set with the union of other focal sets is empty.
		bool is_partitioned() const {
			const std::vector<set_N_value<T, N>* >& elements = this->definition.elements();
			const std::bitset<N>& U = elements[0]->set;

			for (size_t j = 1; j < elements.size(); ++j) {
				const std::bitset<N>& setj = elements[j]->set;
				if (!FOD<N>::are_disjoint(setj, U)) {
					return false;
				}
				U |= setj;
			}
			return U.count() == N;
		}

		/// Mass function without internal conflict is a one where all pairs of
		/// focal elements have a non-empty intersection.
		/// \f[
		/// \forall A, B \subseteq \Omega, m(A) > 0, m(B) > 0 : A \cap B \neq \emptyset
		/// \f]
		/// Therefore, there is internal conflict when there is at least one empty intersection between pairs of focal sets,
		/// i.e. when the intersection of all focal elements is empty.
		bool has_internal_conflict() const {
			const std::vector<set_N_value<T, N>* >& elements = this->definition.elements();
			const std::bitset<N>& I = elements[0]->set;

			for (size_t j = 1; j < elements.size(); ++j) {
				const std::bitset<N>& setj = elements[j]->set;
				I &= setj;
				if (FOD<N>::is_emptyset(I)){
					return true;
				}
			}
			return false;
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_MASS_HPP
