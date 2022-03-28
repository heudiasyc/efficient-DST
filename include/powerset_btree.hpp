#ifndef EFFICIENT_DST_POWERSET_BTREE_HPP
#define EFFICIENT_DST_POWERSET_BTREE_HPP

#include "macros.hpp"

#include <vector>
#include <string>
#include <vector>
#include <time.h>
#include <iomanip>
#include <map>
#include <sstream>

#include <sample_space.hpp>


namespace efficient_DST{

	// allow user to configure the floating-point tolerance
	static constexpr float precision = 5e-7;
	static constexpr int precision_digits = (int) -log10(precision);

	template <size_t N, class T = float>
	struct set_N_value{
		bool is_null;
		typename sample_space<N>::subset set;
		size_t cardinality;
		T value;

		set_N_value(const typename sample_space<N>::subset& _set, T _value) :
			is_null(false),
			set(_set),
			cardinality(_set.count()),
			value(_value)
		{}

		set_N_value(const typename sample_space<N>::subset& _set) :
			is_null(true),
			set(_set),
			cardinality(_set.count()),
			value(0)
		{}

		set_N_value() :
			is_null(true),
			set(0),
			cardinality(0),
			value(0)
		{}


		static inline const std::string to_string(const T& value){
			int displayed_precision = precision_digits;
			int magnitude;
			if (value < 0){
				--displayed_precision;
				magnitude = (int) log10(-value);
			}else{
				magnitude = (int) log10(value);
			}
			if (magnitude > 0)
				displayed_precision -= magnitude;
			std::ostringstream stm;
			stm << std::fixed;
			stm << std::setprecision(displayed_precision < 0 ? 0 : displayed_precision);
			stm << value;
			return stm.str() ;
		}

		const std::string to_string(const sample_space<N>& outcomes) const {
			if (this->is_null)
				return "NULL\t <- " + outcomes.to_string(this->set);
			else
				return to_string(this->value) + "\t <- " + outcomes.to_string(this->set);
		}

		const std::string to_string() const {
			if (this->is_null)
				return "NULL\t <- " + this->set.to_string();
			else
				return to_string(this->value) + "\t <- " + this->set.to_string();
		}
	};

	template <size_t N, class T = float>
	struct node : public set_N_value<N, T>{
		size_t depth;			// having a copy instead of a pointer ensures best performance as the depth
								// (that is heavily used in searches)
								// is at the same place in memory than all other information necessary to browse the tree
		size_t parent;
		size_t left;
		size_t right;

		node(const typename sample_space<N>::subset& _set, T _value) :
			set_N_value<N, T>(_set, _value),
			depth(0),
			parent(0),
			left(0),
			right(0)
		{}

		node(	const typename sample_space<N>::subset& _set,
				T _value,
				size_t _depth,
				size_t _parent_node,
				size_t _left_node,
				size_t _right_node) :
			set_N_value<N, T>(_set, _value),
			depth(_depth),
			parent(_parent_node),
			left(_left_node),
			right(_right_node)
		{}

		node(const typename sample_space<N>::subset& _set) :
			set_N_value<N, T>(_set),
			depth(0),
			parent(0),
			left(0),
			right(0)
		{}

		node(	const typename sample_space<N>::subset& _set,
				size_t _depth,
				size_t _parent_node,
				size_t _left_node,
				size_t _right_node) :
			set_N_value<N, T>(_set),
			depth(_depth),
			parent(_parent_node),
			left(_left_node),
			right(_right_node)
		{}
	};

	////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/*
	 * Spatial complexity O(P + D), where P = number of non null values manually assigned for subsets of sample space
	 * and D = number of disjunction nodes in tree
	 * However, a disjunction node is always created to separate two regular nodes (when needed), so there are at most P/2 disjunction nodes
	 * => O(P)
	 */
	template <size_t N, class T = float>
	class powerset_btree {
	protected:
		typedef typename sample_space<N>::subset subset;
		std::vector<node<N, T> > nodes;
		std::vector<size_t> free_slots;
		size_t number_of_non_null_values = 0;
		size_t cardinality_distribution_non_null[N+1] = {0};
		std::unordered_map<subset, size_t> manifest;

	public:
//		const size_t NULL_INDEX = 1 << N;
		const size_t INIT_SIZE = 32;

		powerset_btree(const size_t& init_size = 32) : INIT_SIZE(init_size)
		{
			init_tree();
		}

		powerset_btree(const powerset_btree<N, T>& p) :
			nodes(p.nodes),
			free_slots(p.free_slots),
			number_of_non_null_values(p.number_of_non_null_values),
			manifest(p.manifest)
		{
			std::copy(
				std::begin(p.cardinality_distribution_non_null),
				std::end(p.cardinality_distribution_non_null),
				std::begin(cardinality_distribution_non_null)
			);
		}


		/////////////////////////////////////////

		std::ostream& print(const sample_space<N>& outcomes, const bool& including_null = false) const {
			const std::vector<set_N_value<N, T> const * >& indices = this->elements(including_null);
			std::cout << std::endl;
			for (size_t i = 0; i < indices.size(); ++i) {
				std::cout << indices[i]->to_string(outcomes) << std::endl;
			}
//			this->print_layout();
			return std::cout;
		}

		std::ostream& print(const bool& including_null = false) const {
			const std::vector<set_N_value<N, T>* >& indices = this->elements(including_null);
			std::cout << std::endl;
			for (size_t i = 0; i < indices.size(); ++i) {
				std::cout << indices[i]->to_string() << std::endl;
			}

			return std::cout;
		}

		/////////////////////////////////////////

		void operator=(const powerset_btree<N, T>& p){
			nodes = p.nodes;
			free_slots = p.free_slots;
			number_of_non_null_values = p.number_of_non_null_values;
			manifest = p.manifest;
			std::copy(
				std::begin(p.cardinality_distribution_non_null),
				std::end(p.cardinality_distribution_non_null),
				std::begin(cardinality_distribution_non_null)
			);
		}
//
//		void copy_sets(const powerset_btree<N, T>& p, T default_value){
//			std::vector<set_N_value<N, T>* > elem = p.elements();
//			for (size_t i = 0; i < elem.size(); ++i) {
//				set_N_value<N, T>* node = (*this)[elem[i]->set];
//				if(!node)
//					insert(elem[i]->set, default_value);
//			}
//		}

		void reserve(const size_t& n){
			nodes.reserve(n);
			manifest.reserve(n);
		}

		size_t number_of_nodes() const {
			return this->nodes.size();
		}

		const size_t& size() const {
			return this->number_of_non_null_values;
		}

		const size_t& get_nb_sets_of_cardinality(const size_t card) const {
			return this->cardinality_distribution_non_null[card];
		}

		void nullify(size_t index){

			if(index >= nodes.size())
				return;
			node<N, T>& sentenced_node = nodes[index];

			if(sentenced_node.is_null)
				// if n is already NULL, then it is simply a disjunction node => nothing to do
				return;

			if(index == 0){
				nullify_without_repercussion(sentenced_node);
				return;
			}

			if(sentenced_node.left && sentenced_node.right){
				// if n has 2 children, transform it into a disjunction node
				nullify_without_repercussion(sentenced_node);
				//insert_disjunction_node(n);
			}else if(sentenced_node.left || sentenced_node.right){
				// if n has exactly one child, then erase n and link its child to its parent
				size_t child_index;
				if(sentenced_node.right)
					child_index = sentenced_node.right;
				else
					child_index = sentenced_node.left;

				replace_node_with(index, child_index);
			}else{
				// if n has no child, just erase it
				erase_node_without_child(index);
			}
		}

		void nullify(){
			this->manifest.clear();
			this->nodes.clear();
			this->free_slots.clear();
			this->number_of_non_null_values = 0;
			for (size_t c = 0; c <= N; ++c){
				this->cardinality_distribution_non_null[c] = 0;
			}
			init_tree();
		}

		/////////////////////////////////////////

		size_t insert(const subset& set, const T& value) {
			node<N, T>& emptyset_node = nodes[0];
			if (set == 0){
				if(emptyset_node.is_null){
					emptyset_node.is_null = false;
					++this->number_of_non_null_values;
					this->cardinality_distribution_non_null[0] = 1;
				}
				emptyset_node.value = value;
				return 0;
			}

			const size_t& final_depth = get_final_element_number(set);
			size_t depth = 0;
			subset cursor = 1;
			subset mask = 1;
			size_t index = 1;
			size_t inserted_index;
			bool is_left_child;

			for(;;){
				node<N, T>& leaf = nodes[index];
				if(depth < leaf.depth){
					if (leaf.set == set){
						if(leaf.is_null){
							leaf.is_null = false;
							++this->number_of_non_null_values;
							++this->cardinality_distribution_non_null[leaf.cardinality];
						}
						leaf.value = value;
						return index;
					}else{
						// if leaf->set is not equal to set
						// take skipped depths into account
						const subset& div_bits = set ^ leaf.set;
						size_t min_depth = std::min(final_depth, leaf.depth);
						while(depth < min_depth && (div_bits & cursor) == 0){
							cursor <<= 1;
							mask |= cursor;
							++depth;
						}
						if(depth < leaf.depth){
							// Since depth < leaf->depth, there is a disjunction before leaf->depth between set and leaf->set.
							// We cannot have depth after final_depth, because that would make it a subset of leaf->set (condition already checked earlier)
							// That being said, we can have depth <= final_depth
							// In any case, create a regular node at final_depth corresponding to set.
							// if depth == final_depth, then put leaf at the left of this new node.
							// otherwise, create a disjunction node at depth
							// and check which of set or leaf->set has depth. The one that has it must be at the right of this node, the other at its left.
							size_t parent_index = leaf.parent;
							if (depth == final_depth){
								if((leaf.set & cursor) == 0){
									// if leaf->set is a superset of set if we ignore the element at final_depth
									inserted_index = create_node(set, value, final_depth, parent_index, index, 0);
								}else{
									// if leaf->set is a superset of set (and is not set itself)
									inserted_index = create_node(set, value, final_depth, parent_index, 0, index);
								}
								leaf.parent = inserted_index;
								if (is_left_child){
									nodes[parent_index].left = inserted_index;
								}else{
									nodes[parent_index].right = inserted_index;
								}
							}else{
								// if depth is both < final_depth and < leaf->depth
								inserted_index = create_node(set, value, final_depth, 0, 0, 0);
								const subset& disjunction_set = mask & set;
								size_t disjunction_index;
								if ((set & cursor) == 0){
									// if leaf->set has the element at depth
									disjunction_index = create_disjunction_node(disjunction_set | cursor, depth, parent_index,
											inserted_index,
											index
									);
								}else{
									// if set has the element at depth
									disjunction_index = create_disjunction_node(disjunction_set, depth, parent_index,
											index,
											inserted_index
									);
								}
								nodes[inserted_index].parent = disjunction_index;
								leaf.parent = disjunction_index;
								if (is_left_child){
									nodes[parent_index].left = disjunction_index;
								}else{
									nodes[parent_index].right = disjunction_index;
								}
							}
							return inserted_index;
						}
					}
				}
				if((set & cursor) != 0){
					if(depth != final_depth){
						if(leaf.right){
							++depth;
							cursor <<= 1;
							mask |= cursor;
							is_left_child = false;
							index = leaf.right;
						}else{
							inserted_index = create_node(set, value, final_depth, index, 0, 0);
							leaf.right = inserted_index;
							return inserted_index;
						}
					}else{
						if(leaf.is_null){
							leaf.is_null = false;
							++this->number_of_non_null_values;
							++this->cardinality_distribution_non_null[leaf.cardinality];
						}
						leaf.value = value;
						return index;
					}
				}else {
					if(leaf.left){
						++depth;
						cursor <<= 1;
						mask |= cursor;
						is_left_child = true;
						index = leaf.left;
					}else{
						inserted_index = create_node(set, value, final_depth, index, 0, 0);
						leaf.left = inserted_index;
						return inserted_index;
					}
				}
			}
			return nodes.size();
		}

		void update(size_t index, const T& value){
			if (index < nodes.size()){
				node<N, T>& inserted_node = nodes[index];
//				std::cout << "Updating set " << inserted_node->set << " with value " << value << "\n";
				if (inserted_node.is_null){
					inserted_node.is_null = false;
					++this->number_of_non_null_values;
					++this->cardinality_distribution_non_null[inserted_node.cardinality];
				}
				inserted_node.value = value;
			}else{
				std::cerr << "No pointer given. Ignoring update.\n";
			}
		}

		size_t insert_set_if_smaller_than(const subset& set, const T& value, const size_t& card){
			size_t index = find(set);
			if (index < nodes.size()){
				node<N, T>& inserted_node = nodes[index];
				if(inserted_node.cardinality > card){
					return nodes.size();
				}else if(inserted_node.is_null){
					inserted_node.is_null = false;
					++this->number_of_non_null_values;
					++this->cardinality_distribution_non_null[inserted_node.cardinality];
					inserted_node.value = value;
				}
			}else{
				if(set.count() <= card){
					index = insert(set, value);
				}else{
					return nodes.size();
				}
			}
			return index;
		}

		size_t insert_set(const subset& set, const T& value){
			size_t index = find(set);
			if (index < nodes.size()){
				set_N_value<N, T>& inserted_node = nodes[index];
				if (inserted_node.is_null){
					inserted_node.is_null = false;
					++this->number_of_non_null_values;
					++this->cardinality_distribution_non_null[inserted_node.cardinality];
					inserted_node.value = value;
				}else{
					return nodes.size();
				}
			}else{
				index = insert(set, value);
			}
			return index;
		}

		size_t update_or_insert(const subset& set, const T& value){
			size_t index = find(set);
			if (index < nodes.size()){
				set_N_value<N, T>& inserted_node = nodes[index];
//				std::cout << "Updating set " << set << " with value " << value << "\n";
				if (inserted_node.is_null){
					inserted_node.is_null = false;
					++this->number_of_non_null_values;
					++this->cardinality_distribution_non_null[inserted_node.cardinality];
				}
				inserted_node.value = value;
			}else{
				index = insert(set, value);
			}
			return index;
		}

		/////////////////////////////////////////

		const node<N, T>& get_node(const size_t& index) const {
			return nodes[index];
		}

		node<N, T>& _node(const size_t& index) {
			return nodes[index];
		}

		size_t operator[](const subset& set) const {
			auto occurrence = this->manifest.find(set);
			if (occurrence == this->manifest.end() || nodes[occurrence->second].is_null){
				return nodes.size();
			} else {
				return occurrence->second;
			}
		}

		size_t find(const subset& set) const {
			auto occurrence = this->manifest.find(set);
			if (occurrence == this->manifest.end()){
				return nodes.size();
			} else {
				return occurrence->second;
			}
		}

		/////////////////////////////////////////

		std::vector<set_N_value<N, T>* > singletons() const {
			size_t index = 1;
			std::vector<set_N_value<N, T>* > singletons;
			singletons.reserve(std::min(N, this->size()));
			subset singleton = 1;
			size_t depth = 0;
			for(;;){
				node<N, T>& cursor = nodes[index];
				singleton <<= cursor.depth - depth;
				depth = cursor.depth;
				if(cursor.set != singleton){
					break;
				}
				if(!cursor.is_null){
					singletons.emplace_back(&cursor);
				}
				if(cursor.left){
					index = cursor.left;
				}else{
					break;
				}
			}
			return singletons;
		}

		size_t empty_set() const {
			if (nodes[0].is_null)
				return nodes.size();
			else
				return 0;
		}

		std::map<size_t, std::vector<set_N_value<N, T> const * >, std::less<size_t> > elements_by_ascending_cardinality() const {
			std::map<size_t, std::vector<set_N_value<N, T> const * >, std::less<size_t> > all_values;
			if(!nodes[0].is_null){
//				all_values.emplace(0);
				all_values[0].emplace_back(&nodes[0]);
			}
			elements_by_ascending_cardinality(1, all_values);
			return all_values;
		}

		std::map<size_t, std::vector<set_N_value<N, T>* >, std::less<size_t> > _elements_by_ascending_cardinality() {
			std::map<size_t, std::vector<set_N_value<N, T>* >, std::less<size_t> > all_values;
			if(!nodes[0].is_null){
//				all_values.emplace(0);
				all_values[0].emplace_back(&nodes[0]);
			}
			_elements_by_ascending_cardinality(1, all_values);
			return all_values;
		}

		std::map<size_t, std::vector<set_N_value<N, T>* >, std::greater<size_t> > _elements_by_descending_cardinality() {
			std::map<size_t, std::vector<set_N_value<N, T>* >, std::greater<size_t> > all_values;
			if(!nodes[0].is_null){
//				all_values.emplace(0, (std::vector<set_N_value<N, T>* >) {&nodes[0]});
				all_values[0].emplace_back(&nodes[0]);
			}
			_elements_by_descending_cardinality(1, all_values);
			return all_values;
		}

		void print_layout() const {
			std::cout << "\nAll elements in tree:\nEMPTYSET : ";
			if(nodes[0].is_null){
				std::cout << "null";
			}else{
				std::cout << nodes[0].value;
			}
			std::cout << "\n[0]\t. : ";
			print_layout(nodes[1]);
			std::cout << std::endl;
		}

		std::vector<set_N_value<N, T>* > _elements(const bool& including_null = false) {
			std::vector<set_N_value<N, T>* > all_values;
			all_values.reserve(including_null ? nodes.size() - free_slots.size() : this->size());
			if(!including_null && nodes[0].is_null){
			}else{
				all_values.emplace_back(&nodes[0]);
			}
			_elements(1, all_values, including_null);
			return all_values;
		}

		std::vector<set_N_value<N, T> const * > elements(const bool& including_null = false) const {
			std::vector<set_N_value<N, T> const * > all_values;
			all_values.reserve(including_null ? nodes.size() - free_slots.size() : this->size());
			if(!including_null && nodes[0].is_null){
			}else{
				all_values.emplace_back(&nodes[0]);
			}
			elements(1, all_values, including_null);
			return all_values;
		}

		std::vector<size_t> elements_indices(const bool& including_null = false) const {
			std::vector<size_t> all_values;
			all_values.reserve(including_null ? nodes.size() - free_slots.size() : this->size());
			if(!including_null && nodes[0].is_null){
			}else{
				all_values.emplace_back(0);
			}
			elements_indices(1, all_values, including_null);
			return all_values;
		}

		std::vector<T> values() const {
			std::vector<T> all_values;
			all_values.reserve(this->size());

			if(!nodes[0].is_null){
				all_values.emplace_back(nodes[0].value);
			}
			values(nodes[1], all_values);
			return all_values;
		}

		/////////////////////////////////////////

		void subsets_of(const subset& set, std::vector<size_t>& subset_values) const {
			if(!nodes[0].is_null){
				subset_values.emplace_back(0);
			}
			if(set != 0){
				subsets_of(set, get_final_element_number(set), subset_values, 0, (subset) 1, 1);
			}
		}

		void subsets_of(const subset& set, std::vector<set_N_value<N, T> const * >& subset_values) const {
			if(!nodes[0].is_null){
				subset_values.emplace_back(&nodes[0]);
			}
			if(set != 0){
				subsets_of(set, get_final_element_number(set), subset_values, 0, (subset) 1, 1);
			}
		}

		std::vector<set_N_value<N, T> const * > subsets_of(const subset& set) const {
			std::vector<set_N_value<N, T> const * > subset_values;
			subset_values.reserve(this->size());
			subsets_of(set, subset_values);
			return subset_values;
		}

		/////////////////////////////////////////

		void supersets_of(const subset& set, std::vector<size_t>& superset_values) const {
			if(set == 0){
				if(!nodes[0].is_null){
					superset_values.emplace_back(0);
				}
				elements_indices(1, superset_values);
			}else{
				supersets_of(set, get_final_element_number(set), superset_values, 0, (subset) 1, (subset) 1, 1);
			}
		}

		void supersets_of(const subset& set, std::vector<set_N_value<N, T> const * >& superset_values) const {
			if(set == 0){
				if(!nodes[0].is_null){
					superset_values.emplace_back(&nodes[0]);
				}
				elements(1, superset_values);
			}else{
				supersets_of(set, get_final_element_number(set), superset_values, 0, (subset) 1, (subset) 1, 1);
			}
		}

		std::vector<set_N_value<N, T> const * > supersets_of(const subset& set) const {
			std::vector<set_N_value<N, T> const * > superset_values;
			superset_values.reserve(this->size());
			supersets_of(set, superset_values);
			return superset_values;
		}

	protected:

		void init_tree(){
			this->reserve(INIT_SIZE);
			create_disjunction_node(subset(0), 0, 0, 0, 0);
			create_disjunction_node(subset(1), 0, 0, 0, 0);
		}

		void destroy_tree(const size_t& index){
			if(index){
				node<N, T>& leaf = nodes[index];
				destroy_tree(leaf.left);
				destroy_tree(leaf.right);

				if(!leaf.is_null){
					nullify_without_repercussion(leaf);
				}
				erase_node(index);
			}
		}

		void erase_node(const size_t& index){
			manifest.erase(nodes[index].set);
			nodes[index].parent = 0;
			free_slots.emplace_back(index);
		}

		void nullify_without_repercussion(node<N, T>& ghost){
			ghost.is_null = true;
			--this->number_of_non_null_values;
			--this->cardinality_distribution_non_null[ghost.cardinality];
		}

		void replace_node_with(const size_t& index, const size_t& chosen_one){
			node<N, T>& sentenced_node = nodes[index];
			size_t parent_index = sentenced_node.parent;

			if(!sentenced_node.is_null){
				nullify_without_repercussion(sentenced_node);
			}

			if(parent_index){
//				if(!keep_root){
//					// if sentenced_node is the tree root and we don't want to keep it, replace it by chosen one
//					this->root = chosen_one;
//					this->root->parent = nullptr;
//					erase_node(sentenced_node);
//				}
//			}else{
				node<N, T>& parent = nodes[parent_index];
				if(chosen_one)
					nodes[chosen_one].parent = parent_index;

				if(parent.left == index)
					parent.left = chosen_one;
				else
					parent.right = chosen_one;

				erase_node(index);
			}
		}

		void erase_node_without_child(const size_t& index){
			node<N, T>& sentenced_node = nodes[index];
			size_t parent_index = sentenced_node.parent;

			if(!parent_index){
				// if sentenced_node is the tree root, just nullify it without erasing it
				if(!sentenced_node.is_null){
					nullify_without_repercussion(sentenced_node);
				}
				//sentenced_node->set = subset(this->fod.size(), 1);
				return;
			}

			node<N, T>& parent = nodes[parent_index];
			if(parent.is_null && parent.parent){
				// if parent is a disjunction node other than root, erase also parent
				size_t parent_other_child;
				if(parent.left == index)
					parent_other_child = parent.right;
				else
					parent_other_child = parent.left;

				replace_node_with(parent_index, parent_other_child);

				if(!sentenced_node.is_null){
					nullify_without_repercussion(sentenced_node);
				}

				erase_node(index);
			}else{
				// if parent is a regular node
				replace_node_with(index, 0);
			}
		}

		/////////////////////////////////////////

		static inline size_t get_final_element_number(const subset& set){
//		    size_t current_bit = set._Find_first();
//		    size_t previous_bit = N;
//		    while (current_bit != N){
//		    	previous_bit = current_bit;
//		    	current_bit = set._Find_next(current_bit);
//		    }
//		    return previous_bit;

		    subset last_bit_cursor = 1;
		    last_bit_cursor <<= (N-1);
		    size_t last_bit_set_reverse_index = 1;
		    while ((set & last_bit_cursor) == 0){
		    	++last_bit_set_reverse_index;
		    	last_bit_cursor >>= 1;
		    }
		    if (last_bit_set_reverse_index > N){
		    	return N;
		    } else {
		    	return N-last_bit_set_reverse_index;
		    }
		}

		/////////////////////////////////////////

		size_t insert_node(const subset& set, const T& value, const size_t& depth, const size_t& parent_node, const size_t& left_node, const size_t right_node){
			size_t index;
			if (free_slots.size() > 0){
				index = free_slots.back();
				free_slots.pop_back();
				nodes[index] = node<N, T>(set, value, depth, parent_node, left_node, right_node);
			}else{
				index = nodes.size();
				nodes.emplace_back(set, value, depth, parent_node, left_node, right_node);
			}

			manifest.emplace(set, index);
			node<N, T>& new_node = nodes[index];
			if(left_node)
				nodes[new_node.left].parent = index;
			if(right_node)
				nodes[new_node.right].parent = index;

			return index;
		}

		size_t create_node(const subset& set, const T& value, const size_t& depth, const size_t& parent_node, const size_t& left_node, const size_t right_node){
			size_t index = insert_node(set, value, depth, parent_node, left_node, right_node);
			++this->number_of_non_null_values;
			++this->cardinality_distribution_non_null[nodes[index].cardinality];
			return index;
		}

		size_t create_disjunction_node(const subset& set, const size_t& depth, const size_t& parent_node, const size_t& left_node, const size_t right_node){
			size_t index = insert_node(set, 0, depth, parent_node, left_node, right_node);
			nodes[index].is_null = true;
			return index;
		}

		/////////////////////////////////////////

//		/*
//		 * 'from' must be the position of a bit evaluating to 1
//		 */
//		static void tile_set_bit(subset& mask, size_t from, size_t to){
//			size_t shift_diff = 1;
//			size_t diff = 1 + to - from;
//			while ((shift_diff << 1) <= diff){
//				mask |= mask << shift_diff;
//				shift_diff <<= 1;
//			}
//			mask |= mask << (diff-shift_diff);
//		}

		/////////////////////////////////////////


		void elements_by_ascending_cardinality(const size_t& index, std::map<size_t, std::vector<set_N_value<N, T> const * >, std::less<size_t> >& all_values) const {
			const node<N, T>& leaf = nodes[index];
			if(!leaf.is_null){
				all_values[leaf.cardinality].emplace_back(&leaf);
			}

			if(leaf.left){
				elements_by_ascending_cardinality(leaf.left, all_values);
			}

			if(leaf.right){
				elements_by_ascending_cardinality(leaf.right, all_values);
			}
		}

		void _elements_by_ascending_cardinality(const size_t& index, std::map<size_t, std::vector<set_N_value<N, T>* >, std::less<size_t> >& all_values) {
			node<N, T>& leaf = nodes[index];
			if(!leaf.is_null){
				all_values[leaf.cardinality].emplace_back(&leaf);
			}

			if(leaf.left){
				_elements_by_ascending_cardinality(leaf.left, all_values);
			}

			if(leaf.right){
				_elements_by_ascending_cardinality(leaf.right, all_values);
			}
		}

		void _elements_by_descending_cardinality(const size_t& index, std::map<size_t, std::vector<set_N_value<N, T>* >, std::greater<size_t> >& all_values) {
			node<N, T>& leaf = nodes[index];
			if(!leaf.is_null){
				all_values[leaf.cardinality].emplace_back(&leaf);
			}

			if(leaf.left){
				_elements_by_descending_cardinality(leaf.left, all_values);
			}

			if(leaf.right){
				_elements_by_descending_cardinality(leaf.right, all_values);
			}
		}


		void print_layout(const node<N, T>& leaf) const {

			if(leaf.is_null)
				std::cout << "null";
			else
				std::cout << leaf.value;


			if(leaf.left){
				const node<N, T>& left_leaf = nodes[leaf.left];
				std::string tab = "";
				for (size_t i = 0; i < leaf.depth; ++i) {
					tab += "|";
				}
				if(leaf.depth > 0)
					tab += "";
				tab += "-";
				for (size_t i = leaf.depth; i < left_leaf.depth; ++i) {
					tab += ".";
				}
				std::cout << "\n" << "[" << left_leaf.depth << "]\t" << tab << " : ";
				print_layout(left_leaf);
			}

			if(leaf.right){
				const node<N, T>& right_leaf = nodes[leaf.right];
				std::string tab = "";
				for (size_t i = 0; i < leaf.depth; ++i) {
					tab += "|";
				}
				if(leaf.depth > 0)
					tab += "";
				tab += "+";
				for (size_t i = leaf.depth; i < right_leaf.depth; ++i) {
					tab += ".";
				}
				std::cout << "\n" << "[" << right_leaf.depth << "]\t" << tab << " : ";
				print_layout(right_leaf);
			}
		}

		void _elements(const size_t& index, std::vector<set_N_value<N, T>* >& all_values, const bool& including_null = false) {
			node<N, T>& leaf = nodes[index];
			if(including_null || !leaf.is_null){
				all_values.emplace_back(&leaf);
			}

			if(leaf.left){
				_elements(leaf.left, all_values, including_null);
			}

			if(leaf.right){
				_elements(leaf.right, all_values, including_null);
			}
		}

		void elements(const size_t& index, std::vector<set_N_value<N, T> const * >& all_values, const bool& including_null = false) const {
			const node<N, T>& leaf = nodes[index];
			if(including_null || !leaf.is_null){
				all_values.emplace_back(&leaf);
			}

			if(leaf.left){
				elements(leaf.left, all_values, including_null);
			}

			if(leaf.right){
				elements(leaf.right, all_values, including_null);
			}
		}


		void elements_indices(const size_t& index, std::vector<size_t>& all_values, const bool& including_null = false) const {
			const node<N, T>& leaf = nodes[index];
			if(including_null || !leaf.is_null){
				all_values.emplace_back(index);
			}

			if(leaf.left){
				elements_indices(leaf.left, all_values, including_null);
			}

			if(leaf.right){
				elements_indices(leaf.right, all_values, including_null);
			}
		}


		void values(const size_t& index, std::vector<T>& all_values) const {
			const node<N, T>& leaf = nodes[index];
			if(!leaf.is_null){
				all_values.emplace_back(leaf.value);
			}
			if(leaf.left){
				values(leaf.left, all_values);
			}
			if(leaf.right){
				values(leaf.right, all_values);
			}
		}


		/*
		 * Complexity O(S), where S is the number of subsets assigned by user of the set defined by key
		 * All subsets of key are of depth <= final_depth and don't contain any other elements than the ones of key
		 */
		void subsets_of(
				const subset& set,
				const size_t& final_depth,
				std::vector<size_t>& subset_values,
				size_t depth,
				subset cursor,
				const size_t& index) const {

			const node<N, T>& leaf = nodes[index];
			if(depth < leaf.depth){
				// take skipped depths into account
				cursor <<= leaf.depth - depth;
				depth = leaf.depth;
				if ((leaf.set & (set | cursor)) != leaf.set){
					// if leaf.set (ignoring current depth) is not a subset of set, then we already know that our search stops here
					return;
				}
			}
			if((set & cursor) != 0){
			//if(set[depth]){
				if(!leaf.is_null){
					subset_values.emplace_back(index);
				}
				if(depth != final_depth){
					++depth;
					cursor <<= 1;
					if(leaf.left){
						subsets_of(set, final_depth, subset_values, depth, cursor, leaf.left);
					}
					if(leaf.right){
						subsets_of(set, final_depth, subset_values, depth, cursor, leaf.right);
					}
				}
			}else {
				if(leaf.left){
					++depth;
					cursor <<= 1;
					subsets_of(set, final_depth, subset_values, depth, cursor, leaf.left);
				}
			}
		}


		void subsets_of(
				const subset& set,
				const size_t& final_depth,
				std::vector<set_N_value<N, T> const * >& subset_values,
				size_t depth,
				subset cursor,
				const size_t& index) const {

			const node<N, T>& leaf = nodes[index];
			if(depth < leaf.depth){
				// take skipped depths into account
				cursor <<= leaf.depth - depth;
				depth = leaf.depth;
				if ((leaf.set & (set | cursor)) != leaf.set){
					// if leaf.set (ignoring current depth) is not a subset of set, then we already know that our search stops here
					return;
				}
			}
			if((set & cursor) != 0){
			//if(set[depth]){
				if(!leaf.is_null){
					subset_values.emplace_back(&leaf);
				}
				if(depth != final_depth){
					++depth;
					cursor <<= 1;
					if(leaf.left){
						subsets_of(set, final_depth, subset_values, depth, cursor, leaf.left);
					}
					if(leaf.right){
						subsets_of(set, final_depth, subset_values, depth, cursor, leaf.right);
					}
				}
			}else {
				if(leaf.left){
					++depth;
					cursor <<= 1;
					subsets_of(set, final_depth, subset_values, depth, cursor, leaf.left);
				}
			}
		}

//		void subsets_of(
//				const subset& set,
//				const size_t& final_depth,
//				std::unordered_map<size_t, std::vector<set_N_value<N, T>* > >& subset_values,
//				size_t depth,
//				subset cursor,
//				const size_t& index) const {
//
//			const node<N, T>& leaf = nodes[index];
//			if(depth < leaf.depth){
//				// take skipped depths into account
//				cursor <<= leaf.depth - depth;
//				depth = leaf.depth;
//				if ((leaf.set & (set | cursor)) != leaf.set){
//					// if leaf.set (ignoring current depth) is not a subset of set, then we already know that our search stops here
//					return;
//				}
//			}
//			if((set & cursor) != 0){
//			//if(set[depth]){
//				if(!leaf.is_null){
//					subset_values[leaf.cardinality].emplace_back(&leaf);
//				}
//				if(depth != final_depth){
//					++depth;
//					cursor <<= 1;
//					if(leaf.left){
//						subsets_of(set, final_depth, subset_values, depth, cursor, leaf.left);
//					}
//					if(leaf.right){
//						subsets_of(set, final_depth, subset_values, depth, cursor, leaf.right);
//					}
//				}
//			}else {
//				if(leaf.left){
//					++depth;
//					cursor <<= 1;
//					subsets_of(set, final_depth, subset_values, depth, cursor, leaf.left);
//				}
//			}
//		}

		/////////////////////////////////////////

		/*
		 * Complexity O(S), where S is the number of subsets assigned by user of the set defined by key
		 * All subsets of key are of depth <= final_depth and don't contain any other elements than the ones of key
		 */
		void supersets_of(
				const subset& set,
				const size_t& final_depth,
				std::vector<size_t>& superset_values,
				size_t depth,
				subset cursor,
				subset mask,
				const size_t& index) const {

			const node<N, T>& leaf = nodes[index];
			if(depth < leaf.depth){
				////////////////////////////
				// take skipped depths into account
				size_t mask_length = depth + 1;
				size_t diff = leaf.depth - depth;
				cursor <<= diff;
				while (mask_length <= diff){
					mask |= mask << mask_length;
					diff -= mask_length;
					mask_length *= 2;
				}
				mask |= mask << diff;

//				subset mask = (const subset) cursor;
//				tile_set_bit(mask, depth, leaf.depth);
//				cursor <<= leaf.depth - depth;
				if ((mask & set & ~(leaf.set | cursor)) != 0){
					return;
				}
				////////////////////////////
				// cursor <<= leaf.depth - depth;
				///////////////////////////
				depth = leaf.depth;
			}
			if(final_depth <= depth){
				////////////////////////////
				// Search for elements except the one at leaf.depth that are in set but not in leaf.set
				//if (set & ~leaf.set != 0){
				//	return;
				//}
				////////////////////////////
				if(!leaf.is_null){
					superset_values.emplace_back(index);
				}
				if (final_depth != depth || ((set & cursor) == 0)){
					++depth;
					cursor <<= 1;
					mask |= cursor;
					if(leaf.left){
						supersets_of(set, final_depth, superset_values, depth, cursor, mask, leaf.left);
					}
					if(leaf.right){
						supersets_of(set, final_depth, superset_values, depth, cursor, mask, leaf.right);
					}
				}else{
					if(leaf.right){
						++depth;
						cursor <<= 1;
						mask |= cursor;
						supersets_of(set, final_depth, superset_values, depth, cursor, mask, leaf.right);
					}
				}
			}else{
				if ((set & cursor) == 0){
					++depth;
					cursor <<= 1;
					mask |= cursor;
					if(leaf.left){
						supersets_of(set, final_depth, superset_values, depth, cursor, mask, leaf.left);
					}
					if(leaf.right){
						supersets_of(set, final_depth, superset_values, depth, cursor, mask, leaf.right);
					}
				}else{
					if(leaf.right){
						++depth;
						cursor <<= 1;
						mask |= cursor;
						supersets_of(set, final_depth, superset_values, depth, cursor, mask, leaf.right);
					}
				}
			}
		}


		void supersets_of(
				const subset& set,
				const size_t& final_depth,
				std::vector<set_N_value<N, T> const * >& superset_values,
				size_t depth,
				subset cursor,
				subset mask,
				const size_t& index) const {

			const node<N, T>& leaf = nodes[index];
			if(depth < leaf.depth){
				////////////////////////////
				// take skipped depths into account
				size_t mask_length = depth + 1;
				size_t diff = leaf.depth - depth;
				cursor <<= diff;
				while (mask_length <= diff){
					mask |= mask << mask_length;
					diff -= mask_length;
					mask_length *= 2;
				}
				mask |= mask << diff;

//				subset mask = (const subset) cursor;
//				tile_set_bit(mask, depth, leaf.depth);
//				cursor <<= leaf.depth - depth;
				if ((mask & set & ~(leaf.set | cursor)) != 0){
					return;
				}
				////////////////////////////
				// cursor <<= leaf.depth - depth;
				///////////////////////////
				depth = leaf.depth;
			}
			if(final_depth <= depth){
				////////////////////////////
				// Search for elements except the one at leaf.depth that are in set but not in leaf.set
				//if (set & ~leaf.set != 0){
				//	return;
				//}
				////////////////////////////
				if(!leaf.is_null){
					superset_values.emplace_back(&leaf);
				}
				if (final_depth != depth || ((set & cursor) == 0)){
					++depth;
					cursor <<= 1;
					mask |= cursor;
					if(leaf.left){
						supersets_of(set, final_depth, superset_values, depth, cursor, mask, leaf.left);
					}
					if(leaf.right){
						supersets_of(set, final_depth, superset_values, depth, cursor, mask, leaf.right);
					}
				}else{
					if(leaf.right){
						++depth;
						cursor <<= 1;
						mask |= cursor;
						supersets_of(set, final_depth, superset_values, depth, cursor, mask, leaf.right);
					}
				}
			}else{
				if ((set & cursor) == 0){
					++depth;
					cursor <<= 1;
					mask |= cursor;
					if(leaf.left){
						supersets_of(set, final_depth, superset_values, depth, cursor, mask, leaf.left);
					}
					if(leaf.right){
						supersets_of(set, final_depth, superset_values, depth, cursor, mask, leaf.right);
					}
				}else{
					if(leaf.right){
						++depth;
						cursor <<= 1;
						mask |= cursor;
						supersets_of(set, final_depth, superset_values, depth, cursor, mask, leaf.right);
					}
				}
			}
		}


//		/*
//		 * Complexity O(S), where S is the number of subsets assigned by user of the set defined by key
//		 * All subsets of key are of depth <= final_depth and don't contain any other elements than the ones of key
//		 */
//		void supersets_of(
//				const subset& set,
//				const size_t& final_depth,
//				std::unordered_map<size_t, std::vector<set_N_value<N, T>* > >& superset_values,
//				size_t depth,
//				subset cursor,
//				subset mask,
//				const size_t& index) const {
//
//			const node<N, T>& leaf = nodes[index];
//			if(depth < leaf.depth){
//				////////////////////////////
//				// take skipped depths into account
//				size_t mask_length = depth + 1;
//				size_t diff = leaf.depth - depth;
//				cursor <<= diff;
//				while (mask_length <= diff){
//					mask |= mask << mask_length;
//					diff -= mask_length;
//					mask_length *= 2;
//				}
//				mask |= mask << diff;
//
////				subset mask = (const subset) cursor;
////				tile_set_bit(mask, depth, leaf.depth);
////				cursor <<= leaf.depth - depth;
//				if ((mask & set & ~(leaf.set | cursor)) != 0){
//					return;
//				}
//				////////////////////////////
//				// cursor <<= leaf.depth - depth;
//				///////////////////////////
//				depth = leaf.depth;
//			}
//			if(final_depth <= depth){
//				////////////////////////////
//				// Search for elements except the one at leaf.depth that are in set but not in leaf.set
//				//if (set & ~leaf.set != 0){
//				//	return;
//				//}
//				////////////////////////////
//				if(!leaf.is_null){
//					superset_values[leaf.cardinality].emplace_back(&leaf);
//				}
//				if (final_depth != depth || ((set & cursor) == 0)){
//					++depth;
//					cursor <<= 1;
//					mask |= cursor;
//					if(leaf.left){
//						supersets_of(set, final_depth, superset_values, depth, cursor, mask, leaf.left);
//					}
//					if(leaf.right){
//						supersets_of(set, final_depth, superset_values, depth, cursor, mask, leaf.right);
//					}
//				}else{
//					if(leaf.right){
//						++depth;
//						cursor <<= 1;
//						mask |= cursor;
//						supersets_of(set, final_depth, superset_values, depth, cursor, mask, leaf.right);
//					}
//				}
//			}else{
//				if ((set & cursor) == 0){
//					++depth;
//					cursor <<= 1;
//					mask |= cursor;
//					if(leaf.left){
//						supersets_of(set, final_depth, superset_values, depth, cursor, mask, leaf.left);
//					}
//					if(leaf.right){
//						supersets_of(set, final_depth, superset_values, depth, cursor, mask, leaf.right);
//					}
//				}else{
//					if(leaf.right){
//						++depth;
//						cursor <<= 1;
//						mask |= cursor;
//						supersets_of(set, final_depth, superset_values, depth, cursor, mask, leaf.right);
//					}
//				}
//			}
//		}
	};
}	// namespace efficient_DST

#endif // EFFICIENT_DST_POWERSET_BTREE_HPP
