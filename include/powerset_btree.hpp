#ifndef EFFICIENT_DST_POWERSET_BTREE_HPP
#define EFFICIENT_DST_POWERSET_BTREE_HPP

#include "macros.hpp"

#include <algorithm>
#include <math.h>
#include <vector>

#include <memory_pool.hpp>
#include <fod.hpp>

namespace efficient_DST{

	template <class T, size_t N>
	struct set_N_value{
		bool is_null;
		std::bitset<N> set;
		size_t cardinality;
		T value;

		set_N_value(std::bitset<N> _set, T _value) :
			is_null(false),
			set(_set),
			cardinality(_set.count()),
			value(_value)
		{}

		set_N_value(std::bitset<N> _set) :
			is_null(true),
			set(_set),
			cardinality(_set.count()),
			value(0)
		{}


		static inline const std::string to_string(const T& n){
			std::ostringstream stm ;
			stm << std::fixed;
			stm << std::setprecision(5);
			stm << n ;
			return stm.str() ;
		}

		const std::string to_string(const FOD<N>& fod) const {
			return to_string(this->value) + "\t <- " + fod.to_string(this->set);
		}
	};

	template <class T, size_t N>
	struct node : public set_N_value<T, N>{
		size_t depth;			// having a copy instead of a pointer ensures best performance as the depth
								// (that is heavily used in searches)
								// is at the same place in memory than all other information necessary to browse the tree
		node<T, N>* parent;
		node<T, N>* left;
		node<T, N>* right;

		node(std::bitset<N>& _set, T _value) :
			set_N_value<T, N>(_set, _value),
			depth(0),
			parent(nullptr),
			left(nullptr),
			right(nullptr)
		{}

		node(	std::bitset<N>& _set,
				T _value,
				size_t _depth,
				node<T, N>* _parent_node,
				node<T, N>* _left_node,
				node<T, N>* _right_node) :
			set_N_value<T, N>(_set, _value),
			depth(_depth),
			parent(_parent_node),
			left(_left_node),
			right(_right_node)
		{}

		node(std::bitset<N>& _set) :
			set_N_value<T, N>(_set),
			depth(0),
			parent(nullptr),
			left(nullptr),
			right(nullptr)
		{}

		node(	std::bitset<N>& _set,
				size_t _depth,
				node<T, N>* _parent_node,
				node<T, N>* _left_node,
				node<T, N>* _right_node) :
			set_N_value<T, N>(_set),
			depth(_depth),
			parent(_parent_node),
			left(_left_node),
			right(_right_node)
		{}
	};

	////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/*
	 * Spatial complexity O(P + D), where P = number of non null values manually assigned for subsets of fod
	 * and D = number of disjunction nodes in tree
	 * However, a disjunction node is always created to separate two regular nodes (when needed), so there are at most P/2 disjunction nodes
	 * => O(P)
	 */
	template <class T, size_t N>
	class powerset_btree {
	protected:
		memory_pool<node<T, N> > node_pool;
		node<T, N> *root;
		set_N_value<T, N> *emptyset;
		size_t number_of_non_null_values;
		FOD<N> *fod;
		size_t block_size;
		size_t cardinality_distribution_non_null[N+1] = {0};

	public:

		/*
		 * block_size is the number of predicted values to store in tree
		 * i.e. the number of nodes that should be reserved (for contiguous objects in memory)
		 * if the number of nodes exceeds block_size, this class will reserve another block of block_size nodes without reallocating the previous one
		 * in order to preserve pointers and references
		 */
		/*powerset_btree(FOD* _fod, const size_t _block_size, std::bitset<N> depths_to_mark) :
			node_pool(_block_size),
			number_of_non_null_values(0),
			fod(_fod),
			block_size(_block_size),
			depths_to_mark(depths_to_mark)
		{
			init_tree();
		}*/

		powerset_btree(FOD<N>* fod, const size_t block_size) :
			node_pool(block_size),
			number_of_non_null_values(0),
			fod(fod),
			block_size(block_size)
		{
			init_tree();
		}

		powerset_btree(const powerset_btree<T, N>& p) :
			powerset_btree(p.fod, p.block_size)//, p.depths_to_mark)
		{
			copy(p);
		}

		// constructor for default initialization (required in unordered_map)
		powerset_btree() :
			node_pool(0),
			root(nullptr),
			emptyset(nullptr),
			number_of_non_null_values(0),
			fod(),
			block_size(0)
		{}

		/////////////////////////////////////////

		~powerset_btree()
		{
			DEBUG(std::clog << "powerset_btree destroyed." << std::endl;);
		}

		std::ostream& print(std::ostream& os) const {
			std::vector<set_N_value<T, N>* > values = this->elements();
			os << std::endl;
			for (size_t i = 0; i < values.size(); ++i) {
				os << values[i]->to_string(*this->fod) << std::endl;
			}

			return os;
		}

		/////////////////////////////////////////

		void init_tree(){
			this->root = create_disjunction_node(std::bitset<N>(1), 0, nullptr, nullptr, nullptr);
			this->emptyset = this->node_pool.emplace(std::bitset<N>(0));
		}

		void copy(const powerset_btree<T, N>& p){

			if(p.fod != this->fod){
				std::vector<set_N_value<T, N>* > elem = p.elements();
				for (size_t i = 0; i < elem.size(); ++i) {
					insert(this->fod->to_labels(elem[i]->set), elem[i]->value);
				}
			}else{
				std::vector<set_N_value<T, N>* > elem = p.elements();
				for (size_t i = 0; i < elem.size(); ++i) {
					insert(elem[i]->set, elem[i]->value);
				}
			}
		}

		void copy_sets(const powerset_btree<T, N>& p, T default_value){

			if(p.fod != this->fod){
				std::vector<set_N_value<T, N>* > elem = p.elements();
				for (size_t i = 0; i < elem.size(); ++i) {
					const std::vector<std::string>& set = this->fod->to_labels(elem[i]->set);
					set_N_value<T, N>* node = (*this)[set];
					if(!node)
						insert(set, default_value);
				}
			}else{
				std::vector<set_N_value<T, N>* > elem = p.elements();
				for (size_t i = 0; i < elem.size(); ++i) {
					set_N_value<T, N>* node = (*this)[elem[i]->set];
					if(!node)
						insert(elem[i]->set, default_value);
				}
			}
		}

		FOD<N>* get_FOD() const {
			return this->fod;
		}

		static const size_t get_FOD_size() {
			return N;
		}

		const size_t& get_block_size() const {
			return this->block_size;
		}

		size_t size() const {
			return this->number_of_non_null_values;
		}

		void nullify(set_N_value<T, N>* s_N_v){

			if(!s_N_v)
				return;

			if(s_N_v->is_null)
				// if n is already NULL, then it is simply a disjunction node => nothing to do
				return;

			if(s_N_v == this->emptyset){
				nullify_without_repercussion(s_N_v);
				return;
			}

			node<T, N>* n = (node<T, N>*) s_N_v;
			if(n->left && n->right){
				// if n has 2 children, transform it into a disjunction node
				nullify_without_repercussion(n);
				//insert_disjunction_node(n);
			}else if(n->left || n->right){
				// if n has exactly one child, then erase n and link its child to its parent
				node<T, N>* child;
				if(n->right)
					child = n->right;
				else
					child = n->left;

				replace_node_with(n, child);
			}else{
				// if n has no child, just erase it
				erase_node_without_child(n);
			}
		}

		void nullify(){
			destroy_tree(this->root->left);
			destroy_tree(this->root->right);
			this->root->is_null = true;
			this->root->left = nullptr;
			this->root->right = nullptr;
			this->emptyset->is_null = true;
			this->number_of_non_null_values = 0;
			for (size_t c = 0; c <= N; ++c){
				this->cardinality_distribution_non_null[c] = 0;
			}
		}

		/////////////////////////////////////////

		set_N_value<T, N>* insert(const std::vector<std::string>& labels, T value){
			return insert(this->fod->to_elements(labels), value);
		}

		set_N_value<T, N>* insert(const std::vector<fod_element*>& fod_elements, T value){
			return insert(this->fod->to_set(fod_elements), value);
		}

		set_N_value<T, N>* insert(const std::bitset<N>& set, const T& value){
			set_N_value<T, N>* target = nullptr;
			size_t depth = 0;
			size_t final_depth = get_final_element_number(set);
			if(final_depth == 0){
				if(this->emptyset->is_null){
					this->emptyset->is_null = false;
					++this->number_of_non_null_values;
					this->cardinality_distribution_non_null[0] = 1;
				}
				this->emptyset->value = value;
				target = this->emptyset;
			}else{
				--final_depth;
				//return insert(set, value, final_depth, nullptr, this->root, depth);

				node<T, N>* parent_node = nullptr;
				bool left_child;
				node<T, N>* cursor = this->root;
				while(cursor){
					// From here, we know that cursor->set[depth] is true since cursor->depth == depth
					if(set[depth]){
						if(depth == final_depth){
							if(cursor->is_null){
								cursor->is_null = false;
								++this->number_of_non_null_values;
								++this->cardinality_distribution_non_null[cursor->cardinality];
							}
							cursor->value = value;
							//mark_subsequent_depths(cursor);
							target = cursor;
							cursor = nullptr;
						}else{
							++depth;
							if(cursor->right){
								parent_node = cursor;
								cursor = cursor->right;
								left_child = false;
							}else{
								/*std::bitset<N> disjunction_set = (const std::bitset<N>&) cursor->set;
								bool mark_depth = false;
								while(depth < final_depth){
									if(this->depths_to_mark[depth]){
										mark_depth = true;
										break;
									}else{
										disjunction_set[depth] = set[depth];
									}
									++depth;
								}

								if(mark_depth){
									disjunction_set[depth] = true;
									cursor->right = create_disjunction_node(disjunction_set, depth, cursor, nullptr, nullptr);
									parent_node = cursor;
									cursor = cursor->right;
									left_child = false;
								}else{*/
								// link this node to the position of the fod element corresponding to *final_depth
								cursor->right = create_node(set, value, final_depth, cursor, nullptr, nullptr);
								target = cursor->right;
								cursor = nullptr;
								//}
							}
						}
					}else {
						++depth;
						if(cursor->left){
							parent_node = cursor;
							cursor = cursor->left;
							left_child = true;
						}else{/*
							std::bitset<N> disjunction_set = (const std::bitset<N>&) cursor->set;
							disjunction_set[depth-1] = false;
							bool mark_depth = false;
							while(depth < final_depth){
								if(this->depths_to_mark[depth]){
									mark_depth = true;
									break;
								}else{
									disjunction_set[depth] = set[depth];
								}
								++depth;
							}

							if(mark_depth){
								disjunction_set[depth] = true;
								cursor->left = create_disjunction_node(disjunction_set, depth, cursor, nullptr, nullptr);
								parent_node = cursor;
								cursor = cursor->left;
								left_child = true;
							}else{*/
							// link this node to the position of the fod element corresponding to *final_depth
							cursor->left = create_node(set, value, final_depth, cursor, nullptr, nullptr);
							target = cursor->left;
							cursor = nullptr;
							//}
						}
					}
					if(cursor){
						// take skipped depths into account
						while(cursor->depth > depth){
							if(set[depth] != std::bitset<N>(cursor->set)[depth]){
								// if there is a disjunction between set and cursor for the element at depth in fod
								if(set[depth]){
									// if key has the first different bit
									if(depth == final_depth){
										// if (cursor U first different bit) & set == set
										if(left_child){
											parent_node->left = create_node(set, value, final_depth, parent_node, cursor, nullptr);
											target = parent_node->left;
										}else{
											parent_node->right = create_node(set, value, final_depth, parent_node, cursor, nullptr);
											target = parent_node->right;
										}
										cursor = nullptr;
										break;
									}else{
										node<T, N>* right_node = create_node(set, value, final_depth, nullptr, nullptr, nullptr);
										std::bitset<N> disjunction_set = set;
										for (size_t k = depth + 1; k < disjunction_set.size(); ++k) {
											disjunction_set[k] = false;
										}
										if(left_child)
											parent_node->left = create_disjunction_node(disjunction_set, depth, parent_node, cursor, right_node);
										else
											parent_node->right = create_disjunction_node(disjunction_set, depth, parent_node, cursor, right_node);
										target = right_node;
										cursor = nullptr;
										break;
									}
								}else{
									// if cursor & set != set and cursor has the first different bit
									// Given that depth < leaf->depth, leaf can't be a parent of the node described by set
									// Hence this only possible case :
									node<T, N>* left_node = create_node(set, value, final_depth, nullptr, nullptr, nullptr);
									std::bitset<N> disjunction_set = cursor->set;
									for (size_t k = depth + 1; k < disjunction_set.size(); ++k) {
										disjunction_set[k] = false;
									}
									if(left_child)
										parent_node->left = create_disjunction_node(disjunction_set, depth, parent_node, left_node, cursor);
									else
										parent_node->right = create_disjunction_node(disjunction_set, depth, parent_node, left_node, cursor);
									target = left_node;
									cursor = nullptr;
									break;
								}
							}else if(depth == final_depth){
									// if cursor & set == set
									if(left_child){
										parent_node->left = create_node(set, value, final_depth, parent_node, nullptr, cursor);
										target = parent_node->left;
									}else{
										parent_node->right = create_node(set, value, final_depth, parent_node, nullptr, cursor);
										target = parent_node->right;
									}
									cursor = nullptr;
									break;
							}
							++depth;
						}
					}
				}
			}
			return target;
		}

		/////////////////////////////////////////

		set_N_value<T, N>* operator[](const std::vector<std::string>& labels) const {
			return operator[](this->fod->to_elements(labels));
		}

		set_N_value<T, N>* operator[](const std::vector<fod_element*>& fod_elements) const {
			return operator[](this->fod->to_set(fod_elements));
		}

		set_N_value<T, N>* operator[](const std::bitset<N>& set) const {
			return find(set, get_final_element_number(set));
		}

		/////////////////////////////////////////

		set_N_value<T, N>* singleton_of_fod_element(const size_t final_depth) const {
			node<T, N> *cursor = this->root;
			size_t depth = 0;
			while (true) {
				if(cursor->depth > final_depth){
					//return this->default_value_copy();
					return nullptr;
				}
				// take skipped depths into account
				while(depth < cursor->depth){
					if(cursor->set[depth]){
						// if set has the element at depth (which isn't final_depth)
						// then there is no non null value for the singleton containing the element at final_deph in fod
						return nullptr;
						//return this->default_value_copy();
					}
					++depth;
				}
				if(depth != final_depth){
					if(cursor->left){
						cursor = cursor->left;
						++depth;
					}else{
						return nullptr;
						//return this->default_value_copy();
					}
				}else{
					if(cursor->is_null)
						return nullptr;
					else
						return cursor;
				}
			}
		}

		set_N_value<T, N>* sub_fod_of_size(const size_t size) const {
			if(size == 0){
				if(this->emptyset->is_null)
					return nullptr;
				else
					return this->emptyset;
			}

			node<T, N> *cursor = this->root;
			size_t depth = 0;
			const size_t final_depth = size-1;
			while (true) {
				if(cursor->depth > final_depth){
					return nullptr;
					//return this->default_value_copy();
				}
				// take skipped depths into account
				while(depth < cursor->depth){
					if(!std::bitset<N>(cursor->set)[depth]){
						// if set hasn't the element at depth
						// then there is no non null value for a fod of size superior to depth
						return nullptr;
						//return this->default_value_copy();
					}
					++depth;
				}
				if(depth != final_depth){
					if(cursor->right){
						cursor = cursor->right;
						++depth;
					}else{
						return nullptr;
						//return this->default_value_copy();
					}
				}else{
					if(cursor->is_null)
						return nullptr;
					else
						return cursor;
				}
			}
		}

		/////////////////////////////////////////

		std::vector<set_N_value<T, N>* > singletons() const {
			node<T, N> *cursor = this->root;
			std::vector<set_N_value<T, N>* > singletons;
			singletons.reserve(N);

			size_t depth = 0;
			while (true) {
				// take skipped depths into account
				while(depth < cursor->depth){
					if(cursor->set[depth]){
						// if set has the element at depth (which isn't its final_depth)
						// then there is no non null value for the singleton containing the element at final_deph in fod
						return singletons;
					}
					++depth;
				}
				if(!cursor->is_null)
					singletons.emplace_back(cursor);
				if(cursor->left){
					cursor = cursor->left;
					++depth;
				}else{
					return singletons;
				}
			}
			return singletons;
		}

		std::vector<set_N_value<T, N>* > sub_fods() const {
			node<T, N> *cursor = this->root;
			std::vector<set_N_value<T, N>* > fods;
			fods.reserve(N+1);

			if(!this->emptyset->is_null)
				fods.emplace_back(this->emptyset);

			size_t depth = 0;
			while (true) {
				// take skipped depths into account
				while(depth < cursor->depth){
					if(!cursor->set[depth]){
						// if set hasn't the element at depth
						// then there is no non null value for a fod of size superior to depth
						return fods;
					}
					++depth;
				}
				if(!cursor->is_null)
					fods.emplace_back(cursor);
				if(cursor->right){
					cursor = cursor->right;
					++depth;
				}else{
					return fods;
				}
			}
			return fods;
		}

		/////////////////////////////////////////

		void fill_with_union_of_powersets(
				const powerset_btree<T, N>& powerset1,
				const powerset_btree<T, N>& powerset2,
				std::function<T(const T&, const T&)> operation,
				const T& default_value
		){
			if(powerset1.fod != powerset2.fod){
				std::cerr << "FOD in powerset1 is not the same as FOD in powerset2. Their union has been aborted." << std::endl;
				return;
			}

			T val1 = default_value, val2 = default_value;
			if(!powerset1.emptyset->is_null)
				val1 = powerset1.emptyset->value;
			if(!powerset2.emptyset->is_null)
				val2 = powerset2.emptyset->value;

			if(this->emptyset->is_null){
				this->emptyset->is_null = false;
				++this->number_of_non_null_values;
				this->cardinality_distribution_non_null[0] = 1;
			}
			this->emptyset->value = operation(val1, val2);

			fill_with_union_of_powersets(
				powerset1.root,
				powerset2.root,
				operation,
				0,
				default_value
			);
		}

		std::unordered_map<size_t, std::vector<set_N_value<T, N>* > > elements_by_set_cardinality() const {
			std::unordered_map<size_t, std::vector<set_N_value<T, N>* > > all_values;
			if(!this->emptyset->is_null){
				all_values.emplace(0, (std::vector<set_N_value<T, N>* >) {this->emptyset});
			}
			elements_by_cardinality(this->root, all_values);
			return all_values;
		}

		std::vector<set_N_value<T, N>* > elements() const {
			std::vector<set_N_value<T, N>* > all_values;
			all_values.reserve(this->size());
			if(!this->emptyset->is_null){
				all_values.emplace_back(this->emptyset);
			}
			elements(this->root, all_values);
			return all_values;
		}

		std::vector<T> values() const {
			std::vector<T> all_values;
			all_values.reserve(this->size());

			if(!this->emptyset->is_null){
				all_values.emplace_back(this->emptyset->value);
			}
			values(this->root, all_values);
			return all_values;
		}

		/////////////////////////////////////////

		std::vector<set_N_value<T, N>* > subsets_of(const std::bitset<N>& set) const {
			std::vector<set_N_value<T, N>* > subset_values;
			subsets_of(set, get_final_element_number(set), false, subset_values);
			return subset_values;

		}

		std::vector<set_N_value<T, N>* > strict_subsets_of(const std::bitset<N>& set) const {
			std::vector<set_N_value<T, N>* > subset_values;
			subsets_of(set, get_final_element_number(set), true, subset_values);
			return subset_values;
		}

		std::unordered_map<size_t, std::vector<set_N_value<T, N>* > > strict_subsets_of_by_cardinality(const std::bitset<N>& set) const {
			std::unordered_map<size_t, std::vector<set_N_value<T, N>* > > subset_values;
			subsets_of_by_cardinality(set, get_final_element_number(set), true, subset_values);
			return subset_values;
		}

		/////////////////////////////////////////

		std::vector<set_N_value<T, N>* > supersets_of(const std::bitset<N>& set) const {
			std::vector<set_N_value<T, N>* > superset_values;
			supersets_of(set, get_final_element_number(set), false, superset_values);
			return superset_values;
		}

		std::unordered_map<size_t, std::vector<set_N_value<T, N>* > > supersets_of_by_cardinality(const std::bitset<N>& set) const {
			std::unordered_map<size_t, std::vector<set_N_value<T, N>* > > superset_values;
			supersets_of_by_cardinality(set, get_final_element_number(set), false, superset_values);
			return superset_values;
		}

		std::vector<set_N_value<T, N>* > strict_supersets_of(const std::bitset<N>& set) const {
			std::vector<set_N_value<T, N>* > superset_values;
			supersets_of(set, get_final_element_number(set), true, superset_values);
			return superset_values;
		}

		std::unordered_map<size_t, std::vector<set_N_value<T, N>* > > strict_supersets_of_by_cardinality(const std::bitset<N>& set) const {
			std::unordered_map<size_t, std::vector<set_N_value<T, N>* > > superset_values;
			supersets_of_by_cardinality(set, get_final_element_number(set), true, false, superset_values);
			return superset_values;
		}

	protected:

		//Complexity O(P), where P = number of nodes assigned by user
		static void fill_with_first_powerset(node<T, N>* leaf, std::function<T(const T&, const T&)> operation, const T& default_value){
			if(!leaf)
				return;

			if(!leaf->is_null)
				insert(leaf->set, operation(leaf->value, default_value));
			fill_with_first_powerset(leaf->left, operation, default_value);
			fill_with_first_powerset(leaf->right, operation, default_value);
		}

		static void fill_with_second_powerset(node<T, N>* leaf, std::function<T(const T&, const T&)> operation, const T& default_value){
			if(!leaf)
				return;

			if(!leaf->is_null)
				insert(leaf->set, operation(default_value, leaf->value));
			fill_with_second_powerset(leaf->left, operation, default_value);
			fill_with_second_powerset(leaf->right, operation, default_value);
		}

		//Complexity O(P1 + P2), where P1 = number of nodes in powerset1, P2 = number of nodes in powerset2
		static void fill_with_union_of_powersets(
				node<T, N>* leaf1,
				node<T, N>* leaf2,
				std::function<T(const T&, const T&)> operation,
				size_t depth,
				const T& default_value){

			// take skipped depths into account
			while(depth < leaf1->depth && depth < leaf2->depth){
				if(leaf1->set[depth] != leaf2->set[depth]){
					fill_with_first_powerset(leaf1, operation, default_value);
					fill_with_second_powerset(leaf2, operation, default_value);
					return;

				}
				++depth;
			}

			if(leaf1->depth < leaf2->depth){
				if(!leaf1->is_null)
					insert(leaf1->set, operation(leaf1->value, default_value));
				if(leaf2->set[depth]){
					if(leaf1->right){
						fill_with_union_of_powersets(leaf1->right, leaf2, operation, depth+1, default_value);
					}
					fill_with_first_powerset(leaf1->left, operation, default_value);
				}else{
					if(leaf1->left){
						fill_with_union_of_powersets(leaf1->left, leaf2, operation, depth+1, default_value);
					}
					fill_with_first_powerset(leaf1->right, operation, default_value);
				}
				return;
			}

			if(leaf1->depth > leaf2->depth){
				if(!leaf2->is_null)
					insert(leaf2->set, operation(default_value, leaf2->value));
				if(leaf1->set[depth]){
					if(leaf2->right){
						fill_with_union_of_powersets(leaf1, leaf2->right, operation, depth+1, default_value);
					}
					fill_with_second_powerset(leaf2->left, operation, default_value);
				}else{
					if(leaf2->left){
						fill_with_union_of_powersets(leaf1, leaf2->left, operation, depth+1, default_value);
					}
					fill_with_second_powerset(leaf2->right, operation, default_value);
				}
				return;
			}

			if(!leaf1->is_null || !leaf2->is_null){
				T val1 = default_value, val2 = default_value;
				if(!leaf1->is_null)
					val1 = leaf1->value;
				if(!leaf2->is_null)
					val2 = leaf2->value;
				insert(leaf1->set, operation(val1, val2));
			}

			++depth;
			if(leaf1->right){
				if(leaf2->right){
					fill_with_union_of_powersets(leaf1->right, leaf2->right, operation, depth, default_value);
				}else{
					fill_with_first_powerset(leaf1->right, operation, default_value);
				}
			}else{
				if(leaf2->right){
					fill_with_second_powerset(leaf2->right, operation, default_value);
				}
			}
			if(leaf1->left){
				if(leaf2->left){
					fill_with_union_of_powersets(leaf1->left, leaf2->left, operation, depth, default_value);
				}else{
					fill_with_first_powerset(leaf1->left, operation, default_value);
				}
			}else{
				if(leaf2->left){
					fill_with_second_powerset(leaf2->left, operation, default_value);
				}
			}
		}


		void destroy_tree(node<T, N> *leaf){
			if(leaf){
				destroy_tree(leaf->left);
				destroy_tree(leaf->right);

				if(!leaf->is_null){
					nullify_without_repercussion(leaf);
				}
				this->node_pool.erase(leaf);
			}
		}

		void nullify_without_repercussion(set_N_value<T, N>* s_N_v){
			s_N_v->is_null = true;
			--this->number_of_non_null_values;
			--this->cardinality_distribution_non_null[s_N_v->cardinality];
		}

		void replace_node_with(node<T, N>* sentenced_node, node<T, N>* chosen_one, bool keep_root=true){
			node<T, N>* parent = sentenced_node->parent;

			if(!sentenced_node->is_null){
				nullify_without_repercussion(sentenced_node);
			}

			if(!parent){
				if(!keep_root){
					// if sentenced_node is the tree root and we don't want to keep it, replace it by chosen one
					this->root = chosen_one;
					this->root->parent = nullptr;
					this->node_pool.erase(sentenced_node);
				}
			}else{
				if(chosen_one)
					chosen_one->parent = parent;

				if(parent->left == sentenced_node)
					parent->left = chosen_one;
				else
					parent->right = chosen_one;

				this->node_pool.erase(sentenced_node);
			}
		}

		void erase_node_without_child(node<T, N> *sentenced_node){
			node<T, N>* parent = sentenced_node->parent;

			if(!parent){
				// if sentenced_node is the tree root, just nullify it without erasing it
				if(!sentenced_node->is_null){
					nullify_without_repercussion(sentenced_node);
				}
				//sentenced_node->set = std::bitset<N>(this->fod.size(), 1);
				return;
			}

			if(parent->is_null && parent->parent){
				// if parent is a disjunction node other than root, erase also parent
				node<T, N>* parent_other_child;
				if(parent->left == sentenced_node)
					parent_other_child = parent->right;
				else
					parent_other_child = parent->left;

				replace_node_with(parent, parent_other_child);

				if(!sentenced_node->is_null){
					nullify_without_repercussion(sentenced_node);
				}

				this->node_pool.erase(sentenced_node);
			}else{
				// if parent is a regular node
				replace_node_with(sentenced_node, nullptr);
			}
		}

		/////////////////////////////////////////

		static inline size_t get_final_element_number(const std::bitset<N>& set){
			size_t final_element_number = 0;
			size_t i = set.size()-1;
			for(; i > 0; --i){
				if(set[i]){
					final_element_number = i + 1;
					break;
				}
			}
			// out of loop because -1 converted to size_t (which is unsigned) gives a positive number which causes an infinite loop
			if(i == 0 && set[i]){
				final_element_number = 1;
			}
			return final_element_number;
		}

		static inline size_t get_final_element_number(const std::vector<fod_element*>& fod_elements){
			size_t final_element_number = 0;
			for(size_t i = 0; i < fod_elements.size(); ++i){
				if(fod_elements[i]->position_in_fod >= final_element_number){
					final_element_number = fod_elements[i]->position_in_fod + 1;
				}
			}
			return final_element_number;
		}

		/////////////////////////////////////////

		node<T, N>* create_node(const std::bitset<N>& set, T value, size_t depth, node<T, N>* parent_node, node<T, N>* left_node, node<T, N>* right_node){
			node<T, N>* new_node = this->node_pool.emplace(set, value, depth, parent_node, left_node, right_node);

			++this->number_of_non_null_values;
			++this->cardinality_distribution_non_null[new_node->cardinality];

			if(left_node)
				new_node->left->parent = new_node;
			if(right_node)
				new_node->right->parent = new_node;
			return new_node;
		}

		node<T, N>* create_disjunction_node(const std::bitset<N>& set, size_t depth, node<T, N>* parent_node, node<T, N>* left_node, node<T, N>* right_node){
			node<T, N>* new_node = this->node_pool.emplace(set, depth, parent_node, left_node, right_node);


			if(left_node)
				new_node->left->parent = new_node;
			if(right_node)
				new_node->right->parent = new_node;
			return new_node;
		}

	/*	void mark_intermediate_depths(node<T, N>* leaf, node<T, N>* deeper_leaf){
			if(leaf->depth+1 >= deeper_leaf->depth){
				return;
			}
			std::bitset<N> disjunction_set = (const std::bitset<N>&) leaf->set;
			node<T, N>* branch_cursor = leaf;

			for(size_t depth = leaf->depth+1; depth < deeper_leaf->depth; ++depth){
				disjunction_set[depth-1] = deeper_leaf->set[depth-1];
				disjunction_set[depth] = true;

				if(deeper_leaf->set[depth-1]){
					branch_cursor->right = create_disjunction_node(disjunction_set, depth, branch_cursor, nullptr, nullptr);
					branch_cursor = branch_cursor->right;
				}else{
					branch_cursor->left = create_disjunction_node(disjunction_set, depth, branch_cursor, nullptr, nullptr);
					branch_cursor = branch_cursor->left;
				}
			}
			if(branch_cursor->depth+1 == deeper_leaf->depth){
				if(deeper_leaf->set[branch_cursor->depth]){
					branch_cursor->right = deeper_leaf;
				}else{
					branch_cursor->left = deeper_leaf;
				}
				deeper_leaf->parent = branch_cursor;
			}
		}*/
/*
		void mark_subsequent_depths(node<T, N>* leaf){
			if(!leaf->right && leaf->depth < this->fod->size()-1){
				std::bitset<N> disjunction_set = (const std::bitset<N>&) leaf->set;
				size_t depth = leaf->depth + 1;
				for(; depth < this->fod->size(); ++depth){
					if(this->depths_to_mark[depth]){
						disjunction_set[depth] = true;
						leaf->right = create_disjunction_node(disjunction_set, depth, leaf, nullptr, nullptr);
						disjunction_set[depth] = false;
						break;
					}
				}
				++depth;
				node<T, N>* branch_cursor = leaf->right;
				for(; depth < this->fod->size(); ++depth){
					if(this->depths_to_mark[depth]){
						disjunction_set[depth] = true;
						branch_cursor->left = create_disjunction_node(disjunction_set, depth, branch_cursor, nullptr, nullptr);
						disjunction_set[depth] = false;
						branch_cursor = branch_cursor->left;
					}
				}
			}
		}
*/
		/////////////////////////////////////////

		/*
		 * Complexity O(N), where N = fod size
		 *
		set_N_value<T, N>* insert(const std::bitset<N>& set, const T& value, const size_t& final_depth,
				node<T, N>* parent_node, node<T, N>*& leaf, size_t& depth){

			if(depth < leaf->depth){
				// take skipped depths into account
				const std::bitset<N>& next_set = leaf->set;
				while(depth < leaf->depth){
					if(set[depth] != next_set[depth]){
						// if there is a disjunction between set and next_set for the element at *depth in fod
						if(set[depth]){
							// if key has the first different bit
							if(depth != final_depth){
								node<T, N> *right_node = create_node(set, value, final_depth, nullptr, nullptr, nullptr);
								std::bitset<N> disjunction_set = set;
								for (size_t k = depth + 1; k < disjunction_set.size(); ++k) {
									disjunction_set[k] = false;
								}
								leaf = create_disjunction_node(disjunction_set, depth, parent_node, leaf, right_node);
								return right_node;
							}else{
								// if (next_set U first different bit) & key == key
								leaf = create_node(set, value, final_depth, parent_node, leaf, nullptr);
								return leaf;
							}
						}else{
							// if next_set & set != set and next_set has the first different bit
							// Given that *depth < leaf->depth, leaf can't be a parent of the node described by set
							// Hence this only possible case :
							node<T, N> *left_node = create_node(set, value, final_depth, nullptr, nullptr, nullptr);
							std::bitset<N> disjunction_set = next_set;
							for (size_t k = depth + 1; k < disjunction_set.size(); ++k) {
								disjunction_set[k] = false;
							}
							leaf = create_disjunction_node(disjunction_set, depth, parent_node, left_node, leaf);
							return left_node;
						}
					}
					if(depth == final_depth){
						// if next_set & set == set
						leaf = create_node(set, value, final_depth, parent_node, nullptr, leaf);
						return leaf;
					}
					++depth;
				}
			}
			//depth = leaf->depth;
			if(set[depth]){
				if(depth != final_depth){
					++depth;
					if(leaf->right){
						return insert(set, value, final_depth, leaf, leaf->right, depth);
					}else{
						//if(this->depths_to_mark[depth]){
						//	std::bitset<N> disjunction_set = (const std::bitset<N>&) leaf->set;
						//	disjunction_set[depth] = true;
						//	leaf->right = create_disjunction_node(disjunction_set, depth, leaf, nullptr, nullptr);
						//	return insert(set, value, final_depth, leaf, leaf->right, depth);
						//}else{
						// link this node to the position of the fod element corresponding to *final_depth
						leaf->right = create_node(set, value, final_depth, leaf, nullptr, nullptr);
						return leaf->right;
						//}
					}
				}else{
					if(leaf->is_null){
						// disjunction node
						//move_disjunction_to_regular_node(leaf);
						leaf->is_null = false;
						++this->number_of_non_null_values;
						//this->non_null_nodes_by_depth[leaf->depth + 1].push_back(leaf);
					}
					leaf->value = value;
					mark_subsequent_depths(leaf);
					return leaf;
				}
			}else {
				++depth;
				if(leaf->left){
					return insert(set, value, final_depth, leaf, leaf->left, depth);
				}else{
					//if(this->depths_to_mark[depth]){
					//	std::bitset<N> disjunction_set = (const std::bitset<N>&) leaf->set;
					//	disjunction_set[depth-1] = false;
					//	disjunction_set[depth] = true;
					//	leaf->left = create_disjunction_node(disjunction_set, depth, leaf, nullptr, nullptr);
					//	return insert(set, value, final_depth, leaf, leaf->left, depth);
					//}else{
					// link this node to the position of the fod element corresponding to *final_depth
					leaf->left = create_node(set, value, final_depth, leaf, nullptr, nullptr);
					return leaf->left;
					//}
				}
			}
		}
*/
		/////////////////////////////////////////
/*
		set_N_value<T, N>* find(const std::bitset<N>& set, const size_t& final_depth, node<T, N> *leaf, size_t& depth, const bool& smallest_superset_search) const {

			if(smallest_superset_search){
				if(depth > final_depth){
					while(leaf->is_null){
						if(leaf->left){
							leaf = leaf->left;
						}else if(leaf->right){
							leaf = leaf->right;
						}else{
							return nullptr;
						}
					}
					return leaf;
				}
			}else{
				// if leaf->depth > final_depth, then leaf has other elements than the ones of key
				if(leaf->depth > final_depth){
					return nullptr;
				}
			}

			// take skipped depths into account
			while(depth < leaf->depth && depth < final_depth){
				if(smallest_superset_search){
					if(set[depth] && !leaf->set[depth])
						return nullptr;
				}else{
					if(set[depth] != leaf->set[depth]){
						// if there is a disjunction between key and next_set for the element at depth in fod
						return nullptr;
					}
				}
				++depth;
			}

			if(set[depth]){
				if(depth == final_depth){
					if(leaf->is_null){
						if(smallest_superset_search && leaf->right){
							++depth;
							return find(set, final_depth, leaf->right, depth, smallest_superset_search);
						}else{
							return nullptr;
						}
					}else{
						return leaf;
					}
				}else{
					if(!leaf->right){
						return nullptr;
					}
					++depth;
					return find(set, final_depth, leaf->right, depth, smallest_superset_search);
				}
			}else {
				if(!leaf->left){
					if(smallest_superset_search && leaf->right){
						++depth;
						return find(set, final_depth, leaf->right, depth, smallest_superset_search);
					}else{
						return nullptr;
					}
				}
				++depth;
				return find(set, final_depth, leaf->left, depth, smallest_superset_search);
			}
		}
*/
		set_N_value<T, N>* find(const std::bitset<N>& set, size_t final_depth) const {

			set_N_value<T, N>* target = nullptr;
			size_t depth = 0;
			if(final_depth == 0){
				if(!this->emptyset->is_null)
					target = this->emptyset;
			}else{
				--final_depth;
				//return find(set, final_depth, this->root, depth, smallest_superset_search);

				node<T, N>* cursor = this->root;
				while(cursor){
					// if cursor->depth > final_depth, then cursor has other elements than the ones of set
					if(cursor->depth > final_depth){
						cursor = nullptr;
						break;
					}
					// take skipped depths into account
					while(cursor->depth > depth){
						if(set[depth] != cursor->set[depth]){
							// if there is a disjunction between set and cursor for the element at depth in fod
							cursor = nullptr;
							break;
						}
						++depth;
					}

					if(cursor){
						// From here, we know that cursor->set[depth] is true since cursor->depth == depth
						if(set[depth]){
							if(depth == final_depth){
								if(!cursor->is_null){
									target = cursor;
								}
								cursor = nullptr;
							}else{
								++depth;
								cursor = cursor->right;
							}
						}else {
							if(cursor->left){
								++depth;
								cursor = cursor->left;
							}else{
								cursor = nullptr;
							}
						}
					}
				}
			}
			return target;
		}

		/////////////////////////////////////////


		static void elements_by_cardinality(node<T, N>* leaf, std::unordered_map<size_t, std::vector<set_N_value<T, N>* > >& all_values) {

			if(!leaf->is_null){
				all_values[leaf->cardinality].emplace_back(leaf);
			}
			DEBUG_TREE({
				if(leaf->is_null)
					std::clog << "null";
				else
					std::clog << leaf->value;
			});


			if(leaf->left){
				DEBUG_TREE({
					std::string tab = "";
					for (size_t i = 0; i < leaf->depth; ++i) {
						tab += "|";
					}
					if(leaf->depth > 0)
						tab += "";
					tab += "-";
					for (size_t i = leaf->depth; i < leaf->left->depth; ++i) {
						tab += ".";
					}
					std::clog << "\n" << "[" << leaf->left->depth << "]\t" << tab << " : ";
				});
				elements_by_cardinality(leaf->left, all_values);
			}

			if(leaf->right){
				DEBUG_TREE({
					std::string tab = "";
					for (size_t i = 0; i < leaf->depth; ++i) {
						tab += "|";
					}
					if(leaf->depth > 0)
						tab += "";
					tab += "+";
					for (size_t i = leaf->depth; i < leaf->right->depth; ++i) {
						tab += ".";
					}
					std::clog << "\n" << "[" << leaf->right->depth << "]\t" << tab << " : ";
				});
				elements_by_cardinality(leaf->right, all_values);
			}
		}


		static void elements(node<T, N>* leaf, std::vector<set_N_value<T, N>* >& all_values) {

			if(!leaf->is_null){
				all_values.emplace_back(leaf);
			}
			DEBUG_TREE({
				if(leaf->is_null)
					std::clog << "null";
				else
					std::clog << leaf->value;
			});


			if(leaf->left){
				DEBUG_TREE({
					std::string tab = "";
					for (size_t i = 0; i < leaf->depth; ++i) {
						tab += "|";
					}
					if(leaf->depth > 0)
						tab += "";
					tab += "-";
					for (size_t i = leaf->depth; i < leaf->left->depth; ++i) {
						tab += ".";
					}
					std::clog << "\n" << "[" << leaf->left->depth << "]\t" << tab << " : ";
				});
				elements(leaf->left, all_values);
			}

			if(leaf->right){
				DEBUG_TREE({
					std::string tab = "";
					for (size_t i = 0; i < leaf->depth; ++i) {
						tab += "|";
					}
					if(leaf->depth > 0)
						tab += "";
					tab += "+";
					for (size_t i = leaf->depth; i < leaf->right->depth; ++i) {
						tab += ".";
					}
					std::clog << "\n" << "[" << leaf->right->depth << "]\t" << tab << " : ";
				});
				elements(leaf->right, all_values);
			}
		}


		static void values(node<T, N>* leaf, std::vector<T>& all_values) {
			if(!leaf->is_null){
				all_values.emplace_back(leaf->value);
			}
			if(leaf->left){
				values(leaf->left, all_values);
			}
			if(leaf->right){
				values(leaf->right, all_values);
			}
		}


		/*
		 * Complexity O(S), where S is the number of subsets assigned by user of the set defined by key
		 * All subsets of key are of depth <= final_depth and don't contain any other elements than the ones of key
		 */
		static void subsets_of(
				const std::bitset<N>& set,
				const size_t& final_depth,
				std::vector<set_N_value<T, N>* >& subset_values,
				size_t depth,
				node<T, N> *leaf,
				bool is_set,
				const bool& strict) {

			// if leaf->depth > *final_depth, then leaf has other elements than the ones of key
			if(leaf->depth > final_depth){
				return;
			}
			// take skipped depths into account
			while(depth < leaf->depth){
				if(!set[depth] && std::bitset<N>(leaf->set)[depth]){
					// if next_set has an element that key doesn't
					return;
				}else if(set[depth] && !std::bitset<N>(leaf->set)[depth]){
					// if key has an element that next_set doesn't
					is_set = false;
				}
				++depth;
			}
			if(set[depth]){
				if(depth != final_depth){
					// get value only if leaf doesn't correspond to key (only strict subsets)
					if(!leaf->is_null){
						subset_values.emplace_back(leaf);
					}
					++depth;
					if(leaf->left){
						subsets_of(set, final_depth, subset_values, depth, leaf->left, false, strict);
					}
					if(leaf->right){
						subsets_of(set, final_depth, subset_values, depth, leaf->right, is_set, strict);
					}
				}else if(!strict || !is_set){//this->fod->is_subset_of(this->fod->to_set(leaf->fod_elements), key)){
					// get value if you don't only want strict subsets
					// OR if this set is a *strict* subset of key
					if(!leaf->is_null){
						subset_values.emplace_back(leaf);
					}
				}
			}else {
				if(leaf->left){
					++depth;
					subsets_of(set, final_depth, subset_values, depth, leaf->left, is_set, strict);
				}
			}
		}

		static void subsets_of_by_cardinality(
				const std::bitset<N>& set,
				const size_t& final_depth,
				std::unordered_map<size_t, std::vector<set_N_value<T, N>* > >& subset_values,
				size_t depth,
				node<T, N> *leaf,
				bool is_set,
				const bool& strict) {

			// if leaf->depth > *final_depth, then leaf has other elements than the ones of key
			if(leaf->depth > final_depth){
				return;
			}
			// take skipped depths into account
			while(depth < leaf->depth){
				if(!set[depth] && leaf->set[depth]){
					// if next_set has an element that key doesn't
					return;
				}else if(set[depth] && !leaf->set[depth]){
					// if key has an element that next_set doesn't
					is_set = false;
				}
				++depth;
			}
			if(set[depth]){
				if(depth != final_depth){
					// get value only if leaf doesn't correspond to key (only strict subsets)
					if(!leaf->is_null){
						subset_values[leaf->cardinality].emplace_back(leaf);
					}
					++depth;
					if(leaf->left){
						subsets_of_by_cardinality(set, final_depth, subset_values, depth, leaf->left, false, strict);
					}
					if(leaf->right){
						subsets_of_by_cardinality(set, final_depth, subset_values, depth, leaf->right, is_set, strict);
					}
				}else if(!strict || !is_set){//this->fod->is_subset_of(this->fod->to_set(leaf->fod_elements), key)){
					// get value if you don't only want strict subsets
					// OR if this set is a *strict* subset of key
					if(!leaf->is_null){
						subset_values[leaf->cardinality].emplace_back(leaf);
					}
				}
			}else {
				if(leaf->left){
					++depth;
					subsets_of_by_cardinality(set, final_depth, subset_values, depth, leaf->left, is_set, strict);
				}
			}
		}

		void subsets_of(const std::bitset<N>& set, size_t final_depth, const bool& strict, std::vector<set_N_value<T, N>* >& subset_values) const {

			if(final_depth == 0){
				if(!strict && !this->emptyset->is_null){
					subset_values.emplace_back(this->emptyset);
				}
			}else{
				subset_values.reserve(this->size());
				if(!this->emptyset->is_null){
					subset_values.emplace_back(this->emptyset);
				}
				--final_depth;
				size_t depth = 0;
				subsets_of(set, final_depth, subset_values, depth, this->root, true, strict);
			}
		}

		void subsets_of_by_cardinality(
				const std::bitset<N>& set, size_t final_depth, const bool& strict, std::unordered_map<size_t, std::vector<set_N_value<T, N>* > >& subset_values) const {

			if(final_depth == 0){
				if(!strict && !this->emptyset->is_null){
					subset_values[0].emplace_back(this->emptyset);
				}
			}else{
				subset_values.reserve(N);

				for (size_t i = 0; i < N; ++i) {
					subset_values[i].reserve(this->cardinality_distribution_non_null[i]);
				}
				--final_depth;
				size_t depth = 0;
				subsets_of_by_cardinality(set, final_depth, subset_values, depth, this->root, true, strict);
			}
		}

		/////////////////////////////////////////

		/*
		 * Complexity O(S), where S is the number of supersets assigned by user of the set defined by key
		 * All supersets of key are of depth >= final_depth and contain every element of key
		 */
		static void supersets_of_by_cardinality(
				const std::bitset<N>& set,
				const size_t& final_depth,
				std::unordered_map<size_t, std::vector<set_N_value<T, N>* > >& superset_values,
				size_t depth,
				node<T, N> *leaf,
				bool is_set,
				const bool& strict
		) {

			// take skipped depths into account
			while(depth < leaf->depth){
				if(set[depth]){
					if(!leaf->set[depth]){
						// if key has an element that next_set doesn't have
						return;
					}
				}else if(leaf->set[depth]){
					// if next_set has an element that key doesn't have
					is_set = false;
				}
				++depth;
			}
			if(depth < final_depth){
				if(set[depth++]){
					// if key has an element at depth, so its supersets
					if(leaf->right){
						supersets_of_by_cardinality(set, final_depth, superset_values, depth, leaf->right, is_set, strict);
					}
				}else{
					// if key has no element at depth, its supersets have one or not
					if(leaf->left){
						supersets_of_by_cardinality(set, final_depth, superset_values, depth, leaf->left, is_set, strict);
					}
					if(leaf->right){
						supersets_of_by_cardinality(set, final_depth, superset_values, depth, leaf->right, false, strict);
					}
				}
			}else if(depth == final_depth){

				if(!strict || !is_set){//!strict || this->fod->is_superset_of(this->fod->to_set(leaf->fod_elements), key)){
					// get value if you don't only want strict supersets
					// OR if this set is a *strict* superset of key
					if(!leaf->is_null){
						superset_values[leaf->cardinality].emplace_back(leaf);
					}
				}
				++depth;
				if(leaf->right){
					supersets_of_by_cardinality(set, final_depth, superset_values, depth, leaf->right, false, strict);
				}
			}else{
				// if depth > final_depth, then we are visiting strict supersets of key
				if(!leaf->is_null){
					superset_values[leaf->cardinality].emplace_back(leaf);
				}

				++depth;
				// Now, no matter what additional elements they have, they are all supersets of key
				if(leaf->left){
					supersets_of_by_cardinality(set, final_depth, superset_values, depth, leaf->left, is_set, strict);
				}
				if(leaf->right){
					supersets_of_by_cardinality(set, final_depth, superset_values, depth, leaf->right, is_set, strict);
				}
			}
		}

		static void supersets_of(
				const std::bitset<N>& set,
				const size_t& final_depth,
				std::vector<set_N_value<T, N>* >& superset_values,
				size_t depth,
				node<T, N> *leaf,
				bool is_set,
				const bool& strict
		) {

			// take skipped depths into account
			while(depth < leaf->depth){
				if(set[depth]){
					if(!std::bitset<N>(leaf->set)[depth]){
						// if key has an element that next_set doesn't
						return;
					}
				}else if(std::bitset<N>(leaf->set)[depth]){
					// if next_set has an element that key doesn't
					is_set = false;
				}
				++depth;
			}
			if(depth < final_depth){
				if(set[depth++]){
					// if key has an element at depth, so its supersets
					if(leaf->right){
						supersets_of(set, final_depth, superset_values, depth, leaf->right, is_set, strict);
					}
				}else{
					// if key has no element at depth, its supersets have one or not
					if(leaf->left){
						supersets_of(set, final_depth, superset_values, depth, leaf->left, is_set, strict);
					}
					if(leaf->right){
						supersets_of(set, final_depth, superset_values, depth, leaf->right, false, strict);
					}
				}
			}else if(depth == final_depth){

				if(!strict || !is_set){//!strict || this->fod->is_superset_of(this->fod->to_set(leaf->fod_elements), key)){
					// get value if you don't only want strict supersets
					// OR if this set is a *strict* superset of key
					if(!leaf->is_null){
						superset_values.emplace_back(leaf);
					}
				}
				++depth;
				if(leaf->right){
					supersets_of(set, final_depth, superset_values, depth, leaf->right, false, strict);
				}
			}else{
				// if depth > final_depth, then we are visiting strict supersets of key
				if(!leaf->is_null){
					superset_values.emplace_back(leaf);
				}

				++depth;
				// Now, no matter what additional elements they have, they are all supersets of key
				if(leaf->left){
					supersets_of(set, final_depth, superset_values, depth, leaf->left, is_set, strict);
				}
				if(leaf->right){
					supersets_of(set, final_depth, superset_values, depth, leaf->right, is_set, strict);
				}
			}
		}

		void supersets_of(const std::bitset<N>& set, size_t final_depth, const bool& strict, std::vector<set_N_value<T, N>* >& superset_values) const {

			superset_values.reserve(this->size());

			if(final_depth == 0){
				if(!strict && !this->emptyset->is_null){
					superset_values.emplace_back(this->emptyset);
				}
				elements(this->root, superset_values);
			}else{
				--final_depth;
				size_t depth = 0;
				supersets_of(set, final_depth, superset_values, depth, this->root, true, strict);
			}
		}

		void supersets_of_by_cardinality(const std::bitset<N>& set, size_t final_depth, const bool& strict, std::unordered_map<size_t, std::vector<set_N_value<T, N>* > >& superset_values) const {

			superset_values.reserve(N);

			for (size_t i = 0; i < N; ++i) {
				superset_values[i].reserve(this->cardinality_distribution_non_null[i]);
			}

			if(final_depth == 0){
				if(!strict && !this->emptyset->is_null){
					superset_values[0].emplace_back(this->emptyset);
				}
				elements_by_cardinality(this->root, superset_values);
			}else{
				--final_depth;
				size_t depth = 0;
				supersets_of(set, final_depth, superset_values, depth, this->root, true, strict);
			}
		}
	};
}	// namespace efficient_DST

#endif // EFFICIENT_DST_POWERSET_BTREE_HPP
