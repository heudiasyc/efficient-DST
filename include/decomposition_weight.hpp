#ifndef EFFICIENT_DST_DECOMPOSITION_WEIGHT_HPP
#define EFFICIENT_DST_DECOMPOSITION_WEIGHT_HPP

#include <mobius_transform.hpp>
#include <mobius_aggregate.hpp>

namespace efficient_DST{

	template <typename T = double>
	class decomposition_weight : public mobius_transform<T>{

	public:

		decomposition_weight(const powerset_btree<T>& focal_log_sets_values) : mobius_transform<T>(focal_log_sets_values)
		{}

		decomposition_weight(const FOD& fod) : mobius_transform<T>(fod)
		{}

		decomposition_weight(const mobius_aggregate<T>& ma) : mobius_transform<T>(ma.inversion(mobius_transformation_form_t::multiplicative))
		{}


		virtual ~decomposition_weight(){}

		void nullify(const std::vector<std::string>& labels) {
			this->definition.nullify(this->definition[labels]);
		}

		void set_values(const std::unordered_map<std::vector<std::string>, T>& values) {
			for (std::pair<std::vector<std::string>, T> labels_U_value : values){
				set_value(labels_U_value.first, labels_U_value.second);
			}
		}

		void set_value(const std::vector<std::string>& labels, T value) {
			this->definition.insert(labels, value);
		}

		void set_emptyset_value(const T& value) {
			this->definition.set_value_of_sub_fod_of_size(0, value);
		}

		void set_fod_value(const T& value) {
			this->definition.set_value_of_sub_fod_of_size(this->definition.fod->size(), value);
		}

		T at_emptyset() const {
			set_N_value<T>* set_value = this->definition.sub_fod_of_size(0);
			if(set_value)
				return set_value->value;
			else
				return 1;
		}

		T at_fod() const {
			set_N_value<T>* set_value = this->definition.sub_fod_of_size(this->definition.fod->size());
			if(set_value)
				return set_value->value;
			else
				return 1;
		}

		T operator[](const std::vector<std::string>& labels) const {
			return find(this->definition.fod->to_set(labels));
		}

		T find(const boost::dynamic_bitset<>& key) const {
			set_N_value<T>* set_value = this->definition[key];
			if(set_value)
				return set_value->value;
			else
				return 1;
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_DECOMPOSITION_WEIGHT_HPP
