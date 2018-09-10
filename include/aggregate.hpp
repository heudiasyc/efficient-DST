#ifndef OW_BFT_AGGREGATE_HPP
#define OW_BFT_AGGREGATE_HPP

#include <bft_function.hpp>

namespace ow_bft{

	/*
	 * aggregate is a bft_function that results in the sum or multiplication of images of another function,
	 * i.e. :
	 * - commonality, belief, implicability, plausibility and pignistic probability as sums of images of a mass function
	 * - conjunctive decomposition as product of images of a commonality function
	 * - disjunctive decomposition as product of images of an implicability function
	 */
	template <typename T = double>
	class aggregate : public bft_function<T> {
	protected:

		virtual T compute_aggregation_at_emptyset() const = 0;

		virtual T compute_aggregation_at_fod() const = 0;

		virtual T compute_aggregation(const boost::dynamic_bitset<>& key) const = 0;

		virtual T compute_aggregation(const std::vector<fod_element*>& fod_elements) const = 0;

	public:

		virtual ~aggregate(){}

		virtual T at_emptyset() const = 0;

		virtual T at_fod() const = 0;

		virtual T operator[](const std::vector<std::string>& labels) const = 0;

		virtual T find(const boost::dynamic_bitset<>& key) const = 0;
	};
}		// namespace ow_bft

#endif // OW_BFT_AGGREGATE_HPP
