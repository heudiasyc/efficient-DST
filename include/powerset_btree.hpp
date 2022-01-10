#ifndef EFFICIENT_DST_POWERSET_BTREE_HPP
#define EFFICIENT_DST_POWERSET_BTREE_HPP

#include "macros.hpp"

#include <algorithm>
#include <math.h>
#include <vector>
#include <string>
#include <vector>
#include <time.h>
#include <boost/functional/hash.hpp>
#include <iomanip>

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

		set_N_value() :
			is_null(true),
			set(0),
			cardinality(0),
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
		std::unordered_map<std::bitset<N>, set_N_value<T, N>* > manifest;

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
			this->emptyset = create_disjunction_node(std::bitset<N>(0), 0, nullptr, nullptr, nullptr);
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

		const size_t& size() const {
			return this->number_of_non_null_values;
		}

		const size_t& get_nb_sets_of_cardinality(const size_t card) const {
			return this->cardinality_distribution_non_null[card];
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
			this->manifest.clear();
			this->node_pool.clear();
			this->number_of_non_null_values = 0;
			for (size_t c = 0; c <= N; ++c){
				this->cardinality_distribution_non_null[c] = 0;
			}
			init_tree();
		}

		/////////////////////////////////////////

		set_N_value<T, N>* insert(const std::vector<std::string>& labels, T value){
			return insert(this->fod->to_elements(labels), value);
		}

		set_N_value<T, N>* insert(const std::vector<fod_element*>& fod_elements, T value){
			return insert(this->fod->to_set(fod_elements), value);
		}

		set_N_value<T, N>* insert(const std::bitset<N>& set, const T& value) {

			if (set == 0){
				if(this->emptyset->is_null){
					this->emptyset->is_null = false;
					++this->number_of_non_null_values;
					this->cardinality_distribution_non_null[0] = 1;
				}
				this->emptyset->value = value;
				return this->emptyset;
			}

			node<T, N>* inserted_node;
			const size_t& final_depth = get_final_element_number(set);
			size_t depth = 0;
			std::bitset<N> cursor = 1;
			node<T, N>* leaf = this->root;
			bool is_left_child;

			for(;;){
				if(depth < leaf->depth){
					// take skipped depths into account
					cursor <<= leaf->depth - depth;
					// Search for elements except the one at leaf->depth that are in leaf->set but not in set
					if ((~(set | cursor) & leaf->set) != 0){
						// if leaf->set (ignoring current depth) is not a subset of set, then we already know that our search stops here
						node<T, N>* parent_node = leaf->parent;
						if ((~leaf->set & set) == 0){
							// if leaf->set is a superset of set
							inserted_node = create_node(set, value, final_depth, parent_node, nullptr, leaf);
							leaf->parent = inserted_node;
							if (is_left_child){
								parent_node->left = inserted_node;
							}else{
								parent_node->right = inserted_node;
							}
						}else{
							// if leaf->set is not a superset of set either
							const std::bitset<N>& div_bits = cursor ^ (set ^ leaf->set);
							const size_t& first_div = div_bits._Find_first();
							// first_div cannot be after leaf->depth, because that would make it a subset of set
							// first div cannot be at leaf->depth either because the cursor set bit is at leaf->depth and we ignored it in the search for first_div
							// Thus, we have first_div < leaf->depth
							// Also, we cannot have first_div after final_depth, because that would make it a subset of leaf->set
							// That being said, we can have first_div == final_depth and first_div < final_depth
							// In any case, create a regular node at final_depth corresponding to set.
							// if first_div == final_depth, then put leaf at the left of this new node.
							// otherwise, create a disjunction node at first_div
							// and check which of set or leaf has first_div. The one that has it must be at the right of this node, the other at its left.
							if (first_div == final_depth){
								inserted_node = create_node(set, value, final_depth, parent_node, leaf, nullptr);
								leaf->parent = inserted_node;
								if (is_left_child){
									parent_node->left = inserted_node;
								}else{
									parent_node->right = inserted_node;
								}
							}else{
								// if first_div is both < final_depth and < leaf->depth
								////////
								std::bitset<N> mask = 1;
								tile_set_bit(mask, 0, first_div);
								const std::bitset<N>& disjunction_set = mask & set;
								////////
								inserted_node = create_node(set, value, final_depth, nullptr, nullptr, nullptr);
								node<T, N>* disjunction_node;
								if ((div_bits & disjunction_set) != 0){
									// if set has the element at first_div
									disjunction_node = create_disjunction_node(disjunction_set, first_div, parent_node,
											inserted_node,
											leaf
									);
								}else{
									// if leaf->set has the element at first_div
									disjunction_node = create_disjunction_node(mask & leaf->set, first_div, parent_node,
											leaf,
											inserted_node
									);
								}
								inserted_node->parent = disjunction_node;
								leaf->parent = disjunction_node;
								if (is_left_child){
									parent_node->left = disjunction_node;
								}else{
									parent_node->right = disjunction_node;
								}
							}
						}
						return inserted_node;
					}
					depth = leaf->depth;
				}
				if((set & cursor) != 0){
					if(depth != final_depth){
						if(leaf->right){
							++depth;
							cursor <<= 1;
							is_left_child = false;
							//insert(inserted_node, set, final_depth, value, depth, cursor, leaf->right, false);
						}else{
							inserted_node = create_node(set, value, final_depth, leaf, nullptr, nullptr);
							leaf->right = inserted_node;
							return inserted_node;
						}
					}else{
						if(leaf->is_null){
							leaf->is_null = false;
							++this->number_of_non_null_values;
							++this->cardinality_distribution_non_null[leaf->cardinality];
						}
						leaf->value = value;
						return leaf;
					}
				}else {
					if(leaf->left){
						++depth;
						cursor <<= 1;
						is_left_child = true;
						//insert(inserted_node, set, final_depth, value, depth, cursor, leaf->left, true);
					}else{
						inserted_node = create_node(set, value, final_depth, leaf, nullptr, nullptr);
						leaf->left = inserted_node;
						return inserted_node;
					}
				}
			}
		}

		set_N_value<T, N>* update_or_insert(const std::bitset<N>& set, const T& value){
			set_N_value<T, N>* inserted_node = find(set);
			if (inserted_node){
				if (inserted_node->is_null){
					inserted_node->is_null = false;
					++this->number_of_non_null_values;
					++this->cardinality_distribution_non_null[inserted_node->cardinality];
				}
				inserted_node->value = value;
			}else{
				inserted_node = insert(set, value);
			}
			return inserted_node;
		}

		/////////////////////////////////////////

		set_N_value<T, N>* operator[](const std::vector<std::string>& labels) const {
			return operator[](this->fod->to_elements(labels));
		}

		set_N_value<T, N>* operator[](const std::vector<fod_element*>& fod_elements) const {
			return operator[](this->fod->to_set(fod_elements));
		}

		set_N_value<T, N>* find(const std::bitset<N>& set) const {
			auto occurrence = this->manifest.find(set);
			if (occurrence == this->manifest.end()){
				return nullptr;
			} else {
				return occurrence->second;
			}
		}

		set_N_value<T, N>* operator[](const std::bitset<N>& set) const {
			auto occurrence = this->manifest.find(set);
			if (occurrence == this->manifest.end() || occurrence->second->is_null){
				return nullptr;
			} else {
				return occurrence->second;
			}
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
					if(!cursor->set[depth]){
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

		std::map<size_t, std::vector<set_N_value<T, N>* > > elements_by_set_cardinality(
				const std::binary_function<size_t, size_t, bool>& comp
			) const {
			std::map<size_t, std::vector<set_N_value<T, N>* >, comp > all_values;
			DEBUG_TREE(std::clog << "\nAll elements in tree by cardinality:\nEMPTYSET : ";);
			if(this->emptyset->is_null){
				DEBUG_TREE(std::clog << "null";);
			}else{
				DEBUG_TREE(std::clog << this->emptyset->value;);
				all_values.emplace(0, (std::vector<set_N_value<T, N>* >) {this->emptyset});
			}
			DEBUG_TREE(std::clog << "\n[0]\t. : ";);
			elements_by_cardinality(this->root, all_values);
			DEBUG_TREE(std::clog << std::endl;);
			return all_values;
		}

//		std::unordered_map<size_t, std::vector<set_N_value<T, N>* > > elements_by_set_cardinality() const {
//			std::unordered_map<size_t, std::vector<set_N_value<T, N>* > > all_values;
//			DEBUG_TREE(std::clog << "\nAll elements in tree by cardinality:\nEMPTYSET : ";);
//			if(this->emptyset->is_null){
//				DEBUG_TREE(std::clog << "null";);
//			}else{
//				DEBUG_TREE(std::clog << this->emptyset->value;);
//				all_values.emplace(0, (std::vector<set_N_value<T, N>* >) {this->emptyset});
//			}
//			DEBUG_TREE(std::clog << "\n[0]\t. : ";);
//			elements_by_cardinality(this->root, all_values);
//			DEBUG_TREE(std::clog << std::endl;);
//			return all_values;
//		}

		std::vector<set_N_value<T, N>* > elements() const {
			std::vector<set_N_value<T, N>* > all_values;
			all_values.reserve(this->size());
			DEBUG_TREE(std::clog << "\nAll elements in tree:\nEMPTYSET : ";);
			if(this->emptyset->is_null){
				DEBUG_TREE(std::clog << "null";);
			}else{
				DEBUG_TREE(std::clog << this->emptyset->value;);
				all_values.emplace_back(this->emptyset);
			}
			DEBUG_TREE(std::clog << "\n[0]\t. : ";);
			elements(this->root, all_values);
			DEBUG_TREE(std::clog << std::endl;);
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
			if(!this->emptyset->is_null){
				subset_values.emplace_back(this->emptyset);
			}
			if(set != 0){
				subset_values.reserve(this->size());
				subsets_of(set, get_final_element_number(set), subset_values, 0, (std::bitset<N>) 1, this->root);
			}
			return subset_values;
		}

		std::unordered_map<size_t, std::vector<set_N_value<T, N>* > > subsets_of_by_cardinality(const std::bitset<N>& set) const {
			std::unordered_map<size_t, std::vector<set_N_value<T, N>* > > subset_values;
			if(!this->emptyset->is_null){
				subset_values[0].emplace_back(this->emptyset);
			}
			if(set != 0){
				subset_values.reserve(N);

				for (size_t i = 0; i < N; ++i) {
					subset_values[i].reserve(this->cardinality_distribution_non_null[i]);
				}
				subsets_of(set, get_final_element_number(set), subset_values, 0, (std::bitset<N>) 1, this->root);
			}
			return subset_values;
		}

		/////////////////////////////////////////

		std::vector<set_N_value<T, N>* > supersets_of(const std::bitset<N>& set) const {
			std::vector<set_N_value<T, N>* > superset_values;
			superset_values.reserve(this->size());

			if(set == 0){
				if(!this->emptyset->is_null){
					superset_values.emplace_back(this->emptyset);
				}
				elements(this->root, superset_values);
			}else{
				supersets_of(set, get_final_element_number(set), superset_values, 0, (std::bitset<N>) 1, this->root);
			}
			return superset_values;
		}

		std::unordered_map<size_t, std::vector<set_N_value<T, N>* > > supersets_of_by_cardinality(const std::bitset<N>& set) const {
			std::unordered_map<size_t, std::vector<set_N_value<T, N>* > > superset_values;
			superset_values.reserve(N);

			for (size_t i = 0; i < N; ++i) {
				superset_values[i].reserve(this->cardinality_distribution_non_null[i]);
			}

			if(set == 0){
				if(!this->emptyset->is_null){
					superset_values[0].emplace_back(this->emptyset);
				}
				elements_by_cardinality(this->root, superset_values);
			}else{
				supersets_of(set, get_final_element_number(set), superset_values, 0, (std::bitset<N>) 1, this->root);
			}
			return superset_values;
		}

	protected:

		//Complexity O(P), where P = number of nodes assigned by user
		void fill_with_first_powerset(node<T, N>* leaf, std::function<T(const T&, const T&)> operation, const T& default_value){
			if(!leaf)
				return;

			if(!leaf->is_null)
				insert(leaf->set, operation(leaf->value, default_value));
			fill_with_first_powerset(leaf->left, operation, default_value);
			fill_with_first_powerset(leaf->right, operation, default_value);
		}

		void fill_with_second_powerset(node<T, N>* leaf, std::function<T(const T&, const T&)> operation, const T& default_value){
			if(!leaf)
				return;

			if(!leaf->is_null)
				insert(leaf->set, operation(default_value, leaf->value));
			fill_with_second_powerset(leaf->left, operation, default_value);
			fill_with_second_powerset(leaf->right, operation, default_value);
		}

		//Complexity O(P1 + P2), where P1 = number of nodes in powerset1, P2 = number of nodes in powerset2
		void fill_with_union_of_powersets(
				node<T, N>* leaf1,
				node<T, N>* leaf2,
				std::function<T(const T&, const T&)> operation,
				size_t depth,
				const T& default_value){

			if(depth < leaf1->depth && depth < leaf2->depth){
				// take skipped depths into account
				--depth;
				//////////////////////
				//size_t first_div = (set ^ cursor->set)._Find_next(depth);
				//////////////////////
				size_t next_bit_leaf1 = leaf1->set._Find_next(depth);
				size_t next_bit_leaf2 = leaf2->set._Find_next(depth);
				while(true){
					if(next_bit_leaf1 != next_bit_leaf2){
						fill_with_first_powerset(leaf1, operation, default_value);
						fill_with_second_powerset(leaf2, operation, default_value);
						return;
					}
					if(next_bit_leaf1 == N){
						break;
					}
					next_bit_leaf1 = leaf1->set._Find_next(next_bit_leaf1);
					next_bit_leaf2 = leaf2->set._Find_next(next_bit_leaf2);
				}
				depth = next_bit_leaf1;
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
				this->manifest.erase(leaf->set);
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
					this->manifest.erase(sentenced_node->set);
					this->node_pool.erase(sentenced_node);
				}
			}else{
				if(chosen_one)
					chosen_one->parent = parent;

				if(parent->left == sentenced_node)
					parent->left = chosen_one;
				else
					parent->right = chosen_one;

				this->manifest.erase(sentenced_node->set);
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

				this->manifest.erase(sentenced_node->set);
				this->node_pool.erase(sentenced_node);
			}else{
				// if parent is a regular node
				replace_node_with(sentenced_node, nullptr);
			}
		}

		/////////////////////////////////////////

		static inline size_t get_final_element_number(const std::bitset<N>& set){
//		    size_t current_bit = set._Find_first();
//		    size_t previous_bit = N;
//		    while (current_bit != N){
//		    	previous_bit = current_bit;
//		    	current_bit = set._Find_next(current_bit);
//		    }
//		    return previous_bit;

		    std::bitset<N> last_bit_cursor = 1;
		    last_bit_cursor <<= (N-1);
		    size_t last_bit_set_reverse_index = 1;
		    while ((set & last_bit_cursor) == 0){
		    	++last_bit_set_reverse_index;
		    }
		    if (last_bit_set_reverse_index > N){
		    	return N;
		    } else {
		    	return N-last_bit_set_reverse_index;
		    }
		}

		/////////////////////////////////////////

		node<T, N>* create_node(const std::bitset<N>& set, T value, size_t depth, node<T, N>* parent_node, node<T, N>* left_node, node<T, N>* right_node){
			node<T, N>* new_node = this->node_pool.emplace(set, value, depth, parent_node, left_node, right_node);
			this->manifest.emplace(set, new_node);

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
			this->manifest.emplace(set, new_node);

			if(left_node)
				new_node->left->parent = new_node;
			if(right_node)
				new_node->right->parent = new_node;
			return new_node;
		}

		/////////////////////////////////////////

		/*
		 * from must be the position of a bit evaluating to 1
		 */
		static void tile_set_bit(std::bitset<N>& mask, size_t from, size_t to){
			size_t shift_diff = 1;
			size_t diff = 1 + to - from;
			while ((shift_diff << 1) <= diff){
				mask |= mask << shift_diff;
				shift_diff <<= 1;
			}
			mask |= mask << (diff-shift_diff);
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
				std::bitset<N> cursor,
				node<T, N> *leaf) {

			// if leaf->depth > final_depth, then leaf has other elements than the ones of key
			if(leaf->depth > final_depth){
				return;
			}

			if(depth < leaf->depth){
				// take skipped depths into account
				cursor <<= leaf->depth - depth;
				// Search for elements except the one at leaf->depth that are in leaf->set but not in set
				if ((~(set | cursor) & leaf->set) != 0){
					return;
				}
				depth = leaf->depth;
			}
			if((set & cursor) != 0){
			//if(set[depth]){
				if(!leaf->is_null){
					subset_values.emplace_back(leaf);
				}
				if(depth != final_depth){
					++depth;
					cursor <<= 1;
					if(leaf->left){
						subsets_of(set, final_depth, subset_values, depth, cursor, leaf->left);
					}
					if(leaf->right){
						subsets_of(set, final_depth, subset_values, depth, cursor, leaf->right);
					}
				}
			}else {
				if(leaf->left){
					++depth;
					cursor <<= 1;
					subsets_of(set, final_depth, subset_values, depth, cursor, leaf->left);
				}
			}
		}

		static void subsets_of(
				const std::bitset<N>& set,
				const size_t& final_depth,
				std::unordered_map<size_t, std::vector<set_N_value<T, N>* > >& subset_values,
				size_t depth,
				std::bitset<N> cursor,
				node<T, N> *leaf) {

			// if leaf->depth > final_depth, then leaf has other elements than the ones of key
			if(leaf->depth > final_depth){
				return;
			}

			if(depth < leaf->depth){
				// take skipped depths into account
				cursor <<= leaf->depth - depth;
				// Search for elements except the one at leaf->depth that are in leaf->set but not in set
				if (~(set | cursor) & leaf->set != 0){
					return;
				}
				depth = leaf->depth;
			}
			if(set & cursor != 0){
			//if(set[depth]){
				if(!leaf->is_null){
					subset_values[leaf->cardinality].emplace_back(leaf);
				}
				if(depth != final_depth){
					++depth;
					cursor <<= 1;
					if(leaf->left){
						subsets_of(set, final_depth, subset_values, depth, cursor, leaf->left);
					}
					if(leaf->right){
						subsets_of(set, final_depth, subset_values, depth, cursor, leaf->right);
					}
				}
			}else {
				if(leaf->left){
					++depth;
					cursor <<= 1;
					subsets_of(set, final_depth, subset_values, depth, cursor, leaf->left);
				}
			}
		}

		/////////////////////////////////////////

		/*
		 * Complexity O(S), where S is the number of subsets assigned by user of the set defined by key
		 * All subsets of key are of depth <= final_depth and don't contain any other elements than the ones of key
		 */
		static void supersets_of(
				const std::bitset<N>& set,
				const size_t& final_depth,
				std::vector<set_N_value<T, N>* >& superset_values,
				size_t depth,
				std::bitset<N> cursor,
				node<T, N> *leaf) {

			if(depth < leaf->depth){
				////////////////////////////
				// take skipped depths into account
				std::bitset<N> mask = (const std::bitset<N>) cursor;
				tile_set_bit(mask, depth, leaf->depth);
				cursor <<= leaf->depth - depth;
				if ((mask & set & ~(leaf->set | cursor)) != 0){
					return;
				}
				////////////////////////////
				// cursor <<= leaf->depth - depth;
				///////////////////////////
				depth = leaf->depth;
			}
			if(final_depth <= depth){
				////////////////////////////
				// Search for elements except the one at leaf->depth that are in set but not in leaf->set
				//if (set & ~leaf->set != 0){
				//	return;
				//}
				////////////////////////////
				if(!leaf->is_null){
					superset_values.emplace_back(leaf);
				}
				if (final_depth != depth || ((set & cursor) == 0)){
					++depth;
					cursor <<= 1;
					if(leaf->left){
						supersets_of(set, final_depth, superset_values, depth, cursor, leaf->left);
					}
					if(leaf->right){
						supersets_of(set, final_depth, superset_values, depth, cursor, leaf->right);
					}
				}else{
					if(leaf->right){
						++depth;
						cursor <<= 1;
						supersets_of(set, final_depth, superset_values, depth, cursor, leaf->right);
					}
				}
			}else{
				if ((set & cursor) == 0){
					++depth;
					cursor <<= 1;
					if(leaf->left){
						supersets_of(set, final_depth, superset_values, depth, cursor, leaf->left);
					}
					if(leaf->right){
						supersets_of(set, final_depth, superset_values, depth, cursor, leaf->right);
					}
				}else{
					if(leaf->right){
						++depth;
						cursor <<= 1;
						supersets_of(set, final_depth, superset_values, depth, cursor, leaf->right);
					}
				}
			}
		}


		/*
		 * Complexity O(S), where S is the number of subsets assigned by user of the set defined by key
		 * All subsets of key are of depth <= final_depth and don't contain any other elements than the ones of key
		 */
		static void supersets_of(
				const std::bitset<N>& set,
				const size_t& final_depth,
				std::unordered_map<size_t, std::vector<set_N_value<T, N>* > >& superset_values,
				size_t depth,
				std::bitset<N> cursor,
				node<T, N> *leaf) {

			if(depth < leaf->depth){
				////////////////////////////
				// take skipped depths into account
				std::bitset<N> mask = (const std::bitset<N>) cursor;
				tile_set_bit(mask, depth, leaf->depth);
				cursor <<= leaf->depth - depth;
				if ((mask & set & ~(leaf->set | cursor)) != 0){
					return;
				}
				////////////////////////////
				// cursor <<= leaf->depth - depth;
				///////////////////////////
				depth = leaf->depth;
			}
			if(final_depth <= depth){
				////////////////////////////
				// Search for elements except the one at leaf->depth that are in set but not in leaf->set
				//if (set & ~leaf->set != 0){
				//	return;
				//}
				////////////////////////////
				if(!leaf->is_null){
					superset_values[leaf->cardinality].emplace_back(leaf);
				}
				if (final_depth != depth || ((set & cursor) == 0)){
					++depth;
					cursor <<= 1;
					if(leaf->left){
						supersets_of(set, final_depth, superset_values, depth, cursor, leaf->left);
					}
					if(leaf->right){
						supersets_of(set, final_depth, superset_values, depth, cursor, leaf->right);
					}
				}else{
					if(leaf->right){
						++depth;
						cursor <<= 1;
						supersets_of(set, final_depth, superset_values, depth, cursor, leaf->right);
					}
				}
			}else{
				if ((set & cursor) == 0){
					++depth;
					cursor <<= 1;
					if(leaf->left){
						supersets_of(set, final_depth, superset_values, depth, cursor, leaf->left);
					}
					if(leaf->right){
						supersets_of(set, final_depth, superset_values, depth, cursor, leaf->right);
					}
				}else{
					if(leaf->right){
						++depth;
						cursor <<= 1;
						supersets_of(set, final_depth, superset_values, depth, cursor, leaf->right);
					}
				}
			}
		}
	};
}	// namespace efficient_DST

#endif // EFFICIENT_DST_POWERSET_BTREE_HPP
