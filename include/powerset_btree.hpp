#ifndef OW_BFT_POWERSET_BTREE_HPP
#define OW_BFT_POWERSET_BTREE_HPP

#include <algorithm>
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <fod.hpp>
#include <memory_pool.hpp>
#include <fod_element.hpp>
#include <powerset_function.hpp>


namespace ow_bft{

	template <class T=double>
	class set_N_value{
	public:
		bool is_null;
		T value;
		std::vector<fod_element*> fod_elements;	// store fod_elements instead of dynamic_bitset because the btree structure is dynamic

		set_N_value(T _value, std::vector<fod_element*> _fod_elements) :
			is_null(false),
			value(_value),
			fod_elements(_fod_elements)
		{}

		set_N_value(std::vector<fod_element*> _fod_elements) :
			is_null(true),
			value(0),
			fod_elements(_fod_elements)
		{}
	};

	template <class T=double>
	class node : public set_N_value<T>{
	public:
		size_t depth;			// having a copy instead of a pointer ensures best performance as the depth
								// (that is heavily used in searches)
								// is at the same place in memory than all other information necessary to browse the tree
		node<T>* parent;
		node<T>* left;
		node<T>* right;

		node(T _value, std::vector<fod_element*> _fod_elements) :
			set_N_value<T>(_value, _fod_elements),
			depth(0),
			parent(nullptr),
			left(nullptr),
			right(nullptr)
		{}

		node(
				T _value,
				size_t _depth,
				node<T>* _parent_node,
				node<T>* _left_node,
				node<T>* _right_node,
				std::vector<fod_element*> _fod_elements) :
			set_N_value<T>(_value, _fod_elements),
			depth(_depth),
			parent(_parent_node),
			left(_left_node),
			right(_right_node)
		{}

		node(std::vector<fod_element*> _fod_elements) :
			set_N_value<T>(_fod_elements),
			depth(0),
			parent(nullptr),
			left(nullptr),
			right(nullptr)
		{}

		node(
				size_t _depth,
				node<T>* _parent_node,
				node<T>* _left_node,
				node<T>* _right_node,
				std::vector<fod_element*> _fod_elements) :
			set_N_value<T>(_fod_elements),
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
	template <class T=double>
	class powerset_btree {
	protected:
		memory_pool<node<T> > node_pool;
		node<T> *root;
		set_N_value<T> *emptyset;
		size_t number_of_non_null_values;

		/*
		non_null_nodes_by_depth structure (cardinality ordering at each depth) :
		{
			{}, 																							: fod = {}
			{c1}, 																							: fod = {c1}
			{c2}, {c1, c2}, 																				: fod = {c1, c2}
			{c3}, {c1, c3}, {c2, c3}, {c1, c2, c3}, 														: fod = {c1, c2, c3}
			{c4}, {c1, c4}, {c2, c4}, {c3, c4}, {c1, c2, c4}, {c1, c3, c4}, {c2, c3, c4}, {c1, c2, c3, c4},	: fod = {c1, c2, c3, c4}
			...
		}

		where the powerset natural order is slightly different :
		{
			{}, 																							: fod = {}
			{c1}, 																							: fod = {c1}
			{c2}, {c1, c2}, 																				: fod = {c1, c2}
			{c3}, {c1, c3}, {c2, c3}, {c1, c2, c3}, 														: fod = {c1, c2, c3}
			{c4}, {c1, c4}, {c2, c4}, {c1, c2, c4}, {c3, c4}, {c1, c3, c4}, {c2, c3, c4}, {c1, c2, c3, c4},	: fod = {c1, c2, c3, c4}
			...
		}
		*/
		//std::vector<std::vector<node<T>* > > non_null_nodes_by_depth;

	public:
		//const T null = -1;
		const FOD *fod;
		const size_t block_size;

		/*
		 * block_size is the number of predicted values to store in tree
		 * i.e. the number of nodes that should be reserved (for contiguous objects in memory)
		 * if the number of nodes exceeds block_size, this class will reserve another block of block_size nodes without reallocating the previous one
		 * in order to preserve pointers and references
		 */
		powerset_btree(const FOD& _fod, const size_t _block_size) :
			node_pool(_block_size),
			number_of_non_null_values(0),
			fod(&_fod),
			block_size(_block_size)
		{
			this->root = create_disjunction_node(0, nullptr, nullptr, nullptr, {this->fod->elements()[0]});

			std::vector<fod_element*> empty_fod_elements;
			this->emptyset = this->node_pool.emplace(empty_fod_elements);
		}

		powerset_btree(const powerset_btree<T>& p) : // @suppress("Class members should be properly initialized")
			powerset_btree(*(p.fod), p.node_pool.get_bock_size())
		{
			copy(p);
		}

		/////////////////////////////////////////

		~powerset_btree(){
			destroy_tree(this->root);
		}

		/////////////////////////////////////////

		void copy(const powerset_btree<T>& p){

			if(p.fod != this->fod){
				std::vector<set_N_value<T>* > elem = p.elements();
				for (size_t i = 0; i < elem.size(); ++i) {
					insert(this->fod->to_labels(elem[i]->fod_elements), elem[i]->value);
				}
			}else{
				std::vector<set_N_value<T>* > elem = p.elements();
				for (size_t i = 0; i < elem.size(); ++i) {
					insert(elem[i]->fod_elements, elem[i]->value);
				}
			}
		}

		size_t size() const {
			/*
			int s = this->non_null_nodes_by_depth[0].size();
			for (int d = 1; d < this->non_null_nodes_by_depth.size(); ++d) {
				s += this->non_null_nodes_by_depth[d].size();
			}
			return s;
			*/
			return this->number_of_non_null_values;
		}

		void nullify(set_N_value<T>* s_N_v){

			if(!s_N_v)
				return;

			if(s_N_v->is_null)
				// if n is already NULL, then it is simply a disjunction node => nothing to do
				return;

			if(s_N_v == this->emptyset){
				s_N_v->is_null = true;
				--this->number_of_non_null_values;
				return;
			}

			node<T>* n = (node<T>*) s_N_v;
			if(n->left && n->right){
				// if n has 2 children, transform it into a disjunction node
				n->is_null = true;
				--this->number_of_non_null_values;
				//insert_disjunction_node(n);
			}else if(n->left || n->right){
				// if n has exactly one child, then erase n and link its child to its parent
				node<T>* child;
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
		}

		/////////////////////////////////////////

		void erase_elements_containing_fod_element(const size_t final_depth){
			/*if(final_depth == 0 && this->size() == 1)
				// avoids to leave this->root as a nullptr
				nullify();
			else{*/
				//erase_nodes_of_depth(final_depth);
			std::clog << "\nROOT";
			erase_nodes_containing_element_at(this->root, &final_depth);
			//}

			// update all depth that was after final_depth
			std::vector<set_N_value<T>* > el = elements();
			size_t i;
			if(this->emptyset->is_null)
				// if emptyset value is null, then it isn't in elements(),
				// so the first element of elements() is of type node<T> and has a depth member
				i = 0;
			else
				// otherwise, it is necessary to skip the first element (emptyset) as it is out of the tree
				i = 1;

			for (; i < el.size(); ++i) {
				node<T>* n = (node<T>*) el[i];
				if(n->depth > final_depth)	// n-> depth == final_depth is impossible since all nodes of this depth have been erased
					--(n->depth);
			}
		}

		/////////////////////////////////////////

		set_N_value<T>* set_value_of_singleton_of_fod_element(const size_t final_depth, T value){
			std::vector<fod_element*> fod_elements{this->fod->elements()[final_depth]};
			return insert(fod_elements, value);
		}

		set_N_value<T>* set_value_of_sub_fod_of_size(const size_t cardinality, T value){

			if(cardinality == 0){
				if(this->emptyset->is_null){
					this->emptyset->is_null = false;
					++this->number_of_non_null_values;
				}
				this->emptyset->value = value;
				return this->emptyset;
			}else{
				std::vector<fod_element*> fod_elements;
				for (size_t j = 0; j < cardinality; ++j) {
					fod_elements.emplace_back(this->fod->elements()[j]);
				}
				return insert(fod_elements, value);
			}
		}

		/////////////////////////////////////////

		set_N_value<T>* insert(const std::vector<std::string>& labels, T value){
			return insert(this->fod->to_elements(labels), value);
		}

		set_N_value<T>* insert(const std::vector<fod_element*>& fod_elements, T value){
			return insert(this->fod->to_set(fod_elements), fod_elements, value);
		}

		set_N_value<T>* insert(const boost::dynamic_bitset<>& key, T value){
			return insert(key, this->fod->to_elements(key), value);
		}

		set_N_value<T>* insert(const boost::dynamic_bitset<>& key, const std::vector<fod_element*>& fod_elements, T value){
			size_t depth = 0;
			size_t final_depth = get_final_element_number(key);
			if(final_depth > 0){
				--final_depth;
				return insert(key, fod_elements, &final_depth, &value, &depth, nullptr, this->root);
			}else{
				if(this->emptyset->is_null){
					this->emptyset->is_null = false;
					++this->number_of_non_null_values;
				}
				this->emptyset->value = value;
				return this->emptyset;
			}
		}

		/////////////////////////////////////////

		set_N_value<T>* operator[](const std::vector<std::string>& labels) const {
			return operator[](this->fod->to_elements(labels));
		}

		set_N_value<T>* operator[](const std::vector<fod_element*>& fod_elements) const {
			return operator[](this->fod->to_set(fod_elements));
		}

		set_N_value<T>* operator[](const boost::dynamic_bitset<>& key) const {
			return find(key, get_final_element_number(key));
		}

		/////////////////////////////////////////

		set_N_value<T>* singleton_of_fod_element(const size_t final_depth) const {
			node<T> *cursor = this->root;
			size_t depth = 0;
			while (true) {
				if(cursor->depth > final_depth){
					//return this->default_value_copy();
					return nullptr;
				}
				if(depth < cursor->depth){
					// take skipped depths into account
					const boost::dynamic_bitset<>& set = this->fod->to_set(cursor->fod_elements);
					while(depth < cursor->depth){
						if(set[depth]){
							// if set has the element at depth (which isn't final_depth)
							// then there is no non null value for the singleton containing the element at final_deph in fod
							return nullptr;
							//return this->default_value_copy();
						}
						++depth;
					}
				}
				//depth = cursor->depth;
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

		set_N_value<T>* sub_fod_of_size(const size_t size) const {
			if(size == 0){
				if(this->emptyset->is_null)
					return nullptr;
				else
					return this->emptyset;
			}

			node<T> *cursor = this->root;
			size_t depth = 0;
			const size_t final_depth = size-1;
			while (true) {
				if(cursor->depth > final_depth){
					return nullptr;
					//return this->default_value_copy();
				}
				if(depth < cursor->depth){
					// take skipped depths into account
					const boost::dynamic_bitset<>& set = this->fod->to_set(cursor->fod_elements);
					while(depth < cursor->depth){
						if(!set[depth]){
							// if set hasn't the element at depth
							// then there is no non null value for a fod of size superior to depth
							return nullptr;
							//return this->default_value_copy();
						}
						++depth;
					}
				}
				//depth = cursor->depth;
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

		std::vector<set_N_value<T>* > singletons() const {
			node<T> *cursor = this->root;
			std::vector<set_N_value<T>* > singletons;
			singletons.reserve(this->fod->size());

			size_t depth = 0;
			while (true) {
				if(depth < cursor->depth){
					// take skipped depths into account
					const boost::dynamic_bitset<>& set = this->fod->to_set(cursor->fod_elements);
					while(depth < cursor->depth){
						if(set[depth]){
							// if set has the element at depth (which isn't its final_depth)
							// then there is no non null value for the singleton containing the element at final_deph in fod
							return singletons;
						}
						++depth;
					}
				}
				//depth = cursor->depth;
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

		std::vector<set_N_value<T>* > sub_fods() const {
			node<T> *cursor = this->root;
			std::vector<set_N_value<T>* > fods;
			fods.reserve(this->fod->size()+1);

			if(!this->emptyset->is_null)
				fods.emplace_back(this->emptyset);

			size_t depth = 0;
			while (true) {
				if(depth < cursor->depth){
					// take skipped depths into account
					const boost::dynamic_bitset<>& set = this->fod->to_set(cursor->fod_elements);
					while(depth < cursor->depth){
						if(!set[depth]){
							// if set hasn't the element at depth
							// then there is no non null value for a fod of size superior to depth
							return fods;
						}
						++depth;
					}
				}
				//depth = cursor->depth;
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
				const powerset_btree<T>& powerset1,
				const powerset_btree<T>& powerset2,
				T (*f)(T, T)
		){
			T val = (*f)(powerset1.emptyset->value, powerset2.emptyset->value);
			if(val > 0)
				this->emptyset->value = val;

			fill_with_union_of_powersets(
					(node<T>*) powerset1.sub_fod_of_size(1),
					(node<T>*) powerset2.sub_fod_of_size(1),
					f,
					powerset2.fod,
					0
				);
		}

		std::unordered_map<size_t, std::vector<set_N_value<T>* > > elements_by_set_cardinality() const {

			//std::vector<std::vector<set_N_value<T>* > > all_values(this->fod->size()+1);
			std::unordered_map<size_t, std::vector<set_N_value<T>* > > all_values;
			// very rough approximation of the number of slots that should be reserved for each cardinality
			/*
			float alloc_per_cardinality = ceil((float) this->node_pool.get_bock_size() / (this->fod->size()-1));
			for (size_t i = 1; i < this->fod->size(); ++i) {
				// always only one emptyset and one FOD
				all_values[i].reserve(alloc_per_cardinality);
			}
			*/
			if(!this->emptyset->is_null){
				std::clog << "\nEMPTYSET : " << this->emptyset->value;
				//all_values[0].emplace_back(this->emptyset);
				all_values.emplace(0, {this->emptyset});
			}

			std::clog << "\n" << "[" << 0 << "]\t" << "ROOT : ";
			elements(this->root, all_values, add_to_values_by_cardinality);
			std::clog << std::endl;

			return all_values;
		}

		std::vector<set_N_value<T>* > elements() const {
			std::vector<set_N_value<T>* > all_values;
			all_values.reserve(this->size());

			if(!this->emptyset->is_null){
				std::clog << "\nEMPTYSET : " << this->emptyset->value;
				all_values.emplace_back(this->emptyset);
			}
			std::clog << "\n" << "[" << 0 << "]\t" << "ROOT : ";
			elements(this->root, all_values, add_to_values);
			std::clog << std::endl;

			return all_values;
		}

		/////////////////////////////////////////

		std::vector<set_N_value<T>* > subsets_of(const boost::dynamic_bitset<>& key) const {
			return subsets_of(key, get_final_element_number(key), false);
		}

		std::vector<set_N_value<T>* > subsets_of(const std::vector<fod_element*>& fod_elements) const {
			return subsets_of(this->fod->to_set(fod_elements));
		}

		std::vector<set_N_value<T>* > subsets_of(const std::vector<std::string>& labels) const {
			return subsets_of(this->fod->to_elements(labels));
		}

		std::vector<set_N_value<T>* > strict_subsets_of(const boost::dynamic_bitset<>& key) const {
			return subsets_of(key, get_final_element_number(key), true);
		}

		std::vector<set_N_value<T>* > strict_subsets_of(const std::vector<fod_element*>& fod_elements) const {
			return strict_subsets_of(this->fod->to_set(fod_elements));
		}

		std::vector<set_N_value<T>* > strict_subsets_of(const std::vector<std::string>& labels) const {
			return strict_subsets_of(this->fod->to_elements(labels));
		}

		std::unordered_map<size_t, std::vector<set_N_value<T>* > > strict_subsets_of_by_cardinality(const std::vector<fod_element*>& fod_elements) const {
			const boost::dynamic_bitset<>& key = this->fod->to_set(fod_elements);
			return subsets_of_by_cardinality(key, get_final_element_number(key), true);
		}

		/////////////////////////////////////////

		std::vector<set_N_value<T>* > supersets_of(const boost::dynamic_bitset<>& key) const {
			return supersets_of(key, get_final_element_number(key), false);
		}

		std::vector<set_N_value<T>* > supersets_of(const std::vector<fod_element*>& fod_elements) const {
			return supersets_of(this->fod->to_set(fod_elements));
		}

		std::vector<set_N_value<T>* > supersets_of(const std::vector<std::string>& labels) const {
			return supersets_of(this->fod->to_elements(labels));
		}

		std::unordered_map<size_t, std::vector<set_N_value<T>* > > supersets_of_by_cardinality(const std::vector<fod_element*>& fod_elements) const {
			const boost::dynamic_bitset<>& key = this->fod->to_set(fod_elements);
			return supersets_of_by_cardinality(key, get_final_element_number(key), false);
		}

		std::vector<set_N_value<T>* > strict_supersets_of(const boost::dynamic_bitset<>& key) const {
			return supersets_of(key, get_final_element_number(key), true);
		}

		std::vector<set_N_value<T>* > strict_supersets_of(const std::vector<fod_element*>& fod_elements) const {
			return strict_supersets_of(this->fod->to_set(fod_elements));
		}

		std::vector<set_N_value<T>* > strict_supersets_of(const std::vector<std::string>& labels) const {
			return strict_supersets_of(this->fod->to_elements(labels));
		}

		std::unordered_map<size_t, std::vector<set_N_value<T>* > > strict_supersets_of_by_cardinality(const std::vector<fod_element*>& fod_elements) const {
			const boost::dynamic_bitset<>& key = this->fod->to_set(fod_elements);
			return supersets_of_by_cardinality(key, get_final_element_number(key), true);
		}

		/////////////////////////////////////////

	private:

		void destroy_tree(node<T> *leaf){
			if(leaf){
				destroy_tree(leaf->left);
				destroy_tree(leaf->right);

				if(!leaf->is_null){
					leaf->is_null = true;
					--this->number_of_non_null_values;
					//erase_node_from_non_null_nodes(leaf);
				}
				this->node_pool.erase(leaf);
			}
		}

		/*
		void erase_node_from_non_null_nodes(node<T>* null_node){
			std::vector<set_N_value<T>* >& nodes_of_same_depth = this->non_null_nodes_by_depth[null_node->depth];
			auto it = std::find(nodes_of_same_depth.begin(), nodes_of_same_depth.end(), null_node);
			nodes_of_same_depth.erase(it);
		}

		void erase_nodes_of_depth_from_non_null_nodes(int depth){
			this->non_null_nodes_by_depth[depth].clear();
			this->non_null_nodes_by_depth.erase(this->non_null_nodes_by_depth.begin() + depth + 1);
		}
		*/

		void replace_node_with(node<T>* sentenced_node, node<T>* chosen_one, bool keep_root=true){
			node<T>* parent = sentenced_node->parent;

			if(!sentenced_node->is_null){
				sentenced_node->is_null = true;
				--this->number_of_non_null_values;
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

		void erase_node_without_child(node<T> *sentenced_node, const size_t root_element_position=0){
			node<T>* parent = sentenced_node->parent;

			if(!parent){
				// if sentenced_node is the tree root, just nullify it without erasing it
				if(!sentenced_node->is_null){
					sentenced_node->is_null = true;
					--this->number_of_non_null_values;
				}

				sentenced_node->fod_elements = {this->fod->elements()[root_element_position]};
				return;
			}

			if(parent->is_null && parent->parent){
				// if parent is a disjunction node other than root, erase also parent
				node<T>* parent_other_child;
				if(parent->left == sentenced_node)
					parent_other_child = parent->right;
				else
					parent_other_child = parent->left;

				replace_node_with(parent, parent_other_child);

				if(!sentenced_node->is_null){
					sentenced_node->is_null = true;
					--this->number_of_non_null_values;
				}

				this->node_pool.erase(sentenced_node);
			}else{
				// if parent is a regular node
				replace_node_with(sentenced_node, nullptr);
			}
		}

		void erase_node_because_of_its_depth(node<T> *sentenced_node, const size_t root_element_position=0) {
			// first node of a branch that doesn't contain the element at sentenced_node->depth in fod
			node<T> *left_child = sentenced_node->left;

			std::clog << "\ndestruction of right child " << std::endl;
			destroy_tree(sentenced_node->right);

			std::clog << "self destruction ";
			if(left_child){
				std::clog << "with left child ";
				// if left_child exists, then replace sentenced_node by left_child
				replace_node_with(sentenced_node, left_child, false);
			}else{
				std::clog << "with nullptr ";
				// if left_child doesn't exists
				erase_node_without_child(sentenced_node, root_element_position);
			}
		}

		void erase_nodes_containing_element_at(node<T>* leaf, const size_t *final_depth){
			if(!leaf){
				std::clog << "null";
				return;
			}

			std::string tab = "";
			for (size_t i = 0; i < leaf->depth; ++i) {
				tab += "\t->";
			}

			if(!leaf->is_null){
				std::clog << "\n" << tab << " " << leaf->value;
			}

			if(leaf->depth > *final_depth){
				// take skipped depths into account
				const boost::dynamic_bitset<>& next_set = this->fod->to_set(leaf->fod_elements);
				if(next_set[*final_depth]){
					// if next_set contains the element at depth in fod
					destroy_tree(leaf->left);
					destroy_tree(leaf->right);
					erase_node_without_child(leaf, *final_depth);
				}
				return;
			}
			if(leaf->depth != *final_depth){
				std::clog << "\n" << tab << " left : ";
				erase_nodes_containing_element_at(leaf->left, final_depth);
				std::clog << "\n" << tab << " right : ";
				erase_nodes_containing_element_at(leaf->right, final_depth);
			}else{
				erase_node_because_of_its_depth(leaf, *final_depth);
			}
		}

		/////////////////////////////////////////

		static size_t get_final_element_number(const boost::dynamic_bitset<>& key){
			size_t final_element_number = 0;
			size_t i = key.size()-1;
			for(; i > 0; --i){
				if(key[i]){
					final_element_number = i + 1;
					break;
				}
			}
			// out of loop because -1 converted to size_t (which is unsigned) gives a positive number which causes an infinite loop
			if(i == 0 && key[i]){
				final_element_number = 1;
			}
			return final_element_number;
		}

		static size_t get_final_element_number(const std::vector<fod_element*>& fod_elements){
			size_t final_element_number = 0;
			for(size_t i = 0; i < fod_elements.size(); ++i){
				if(fod_elements[i]->position_in_fod >= final_element_number){
					final_element_number = fod_elements[i]->position_in_fod + 1;
				}
			}
			return final_element_number;
		}

		/////////////////////////////////////////

		node<T>* create_node(T value, size_t depth, node<T>* parent_node, node<T>* left_node, node<T>* right_node, const std::vector<fod_element*>& fod_elements){
			node<T>* new_node = this->node_pool.emplace(value, depth, parent_node, left_node, right_node, fod_elements);
			//insert_node(new_node);
			++this->number_of_non_null_values;
			//this->non_null_nodes_by_depth[depth+1].push_back(new_node);
			return new_node;
		}
		/*
		void insert_node(node<T>* new_node){

			//int missing_floors = new_node->depth + 1 - this->nodes_by_depth.size();
			//while(missing_floors > 0){
			//	this->nodes_by_depth.emplace_back();
			//	missing_floors -= 1;
			//}

			// ordered insertion to keep the powerset elements ordered by set cardinalities at each depth
			int position = 0;
			std::vector<node<T>* > nodes = this->nodes_by_depth[new_node->depth];
			while(position < nodes.size()
					&& nodes[position].fod_elements.size() < new_node->fod_elements.size()){
				++position;
			}
			if(position >= nodes.size()){
				--position;
			}
			// push instead of emplace because we heavily use pointers and addresses changes in case of reallocation
			this->nodes_by_depth[new_node->depth].push(position, new_node);
		}
		*/
		node<T>* create_disjunction_node(size_t depth, node<T>* parent_node, node<T>* left_node, node<T>* right_node, const std::vector<fod_element*>& fod_elements){
			node<T>* new_node = this->node_pool.emplace(depth, parent_node, left_node, right_node, fod_elements);
			//insert_disjunction_node(new_node);
			return new_node;
		}
		/*
		void insert_disjunction_node(node<T>* new_node){

			//int missing_floors = new_node->depth + 1 - this->disjunction_nodes_by_depth.size();
			//while(missing_floors > 0){
			//	this->disjunction_nodes_by_depth.emplace_back();
			//	missing_floors -= 1;
			//}

			// push instead of emplace because we heavily use pointers and addresses changes in case of reallocation
			this->disjunction_nodes_by_depth[new_node->depth].push_back(new_node);
		}

		void move_disjunction_to_regular_node(node<T>* disjunction_node){
			insert_node(disjunction_node);
			erase_node_from(disjunction_node, this->disjunction_nodes_by_depth);
		}

		void move_regular_to_disjunction_node(node<T>* regular_node){
			insert_disjunction_node(regular_node);
			erase_node_from(regular_node, this->nodes_by_depth);
		}

		void move_regular_to_disjunction_node(node<T>* regular_node, int index){
			insert_disjunction_node(regular_node);
			this->nodes_by_depth[regular_node->depth].erase(index);
		}
		*/
		/////////////////////////////////////////

		/*
		 * Complexity O(N), where N = fod size
		 */
		set_N_value<T>* insert(const boost::dynamic_bitset<>& key, const std::vector<fod_element*>& fod_elements, size_t *final_depth,
				T *value, size_t *depth, node<T>* parent_node, node<T>*& leaf){

			if(*depth < leaf->depth){
				// take skipped depths into account
				const boost::dynamic_bitset<>& next_set = this->fod->to_set(leaf->fod_elements);
				while(*depth < leaf->depth){
					if(key[*depth] != next_set[*depth]){
						// if there is a disjunction between key and next_set for the element at *depth in fod
						if(key[*depth]){
							// if key has the first different bit
							if(*depth != *final_depth){
								node<T> *right_node = create_node(*value, *final_depth, nullptr, nullptr, nullptr, fod_elements);
								boost::dynamic_bitset<> disjunction_set = key;
								for (size_t k = *depth + 1; k < disjunction_set.size(); ++k) {
									disjunction_set[k] = false;
								}
								leaf = create_disjunction_node(*depth, parent_node, leaf, right_node, this->fod->to_elements(disjunction_set));
								leaf->right->parent = leaf;
								leaf->left->parent = leaf;
								return right_node;
							}else{
								// if (next_set U first different bit) & key == key
								leaf = create_node(*value, *final_depth, parent_node, leaf, nullptr, fod_elements);
								leaf->left->parent = leaf;
								return leaf;
							}
						}else{
							// if next_set & key != key and next_set has the first different bit
							// Given that *depth < leaf->depth, leaf can't be a parent of the node described by key
							// Hence this only possible case :
							node<T> *left_node = create_node(*value, *final_depth, nullptr, nullptr, nullptr, fod_elements);
							boost::dynamic_bitset<> disjunction_set = next_set;
							for (size_t k = *depth + 1; k < disjunction_set.size(); ++k) {
								disjunction_set[k] = false;
							}
							leaf = create_disjunction_node(*depth, parent_node, left_node, leaf, this->fod->to_elements(disjunction_set));
							leaf->right->parent = leaf;
							leaf->left->parent = leaf;
							return left_node;
						}
					}
					if(*depth == *final_depth){
						// if next_set & key == key
						leaf = create_node(*value, *final_depth, parent_node, nullptr, leaf, fod_elements);
						leaf->right->parent = leaf;
						return leaf;
					}
					++(*depth);
				}
			}
			//*depth = leaf->depth;
			if(key[*depth]){
				if(*depth != *final_depth){
					if(leaf->right){
						++(*depth);
						return insert(key, fod_elements, final_depth, value, depth, leaf, leaf->right);
					}else{
						// link this node to the position of the fod element corresponding to *final_depth
						leaf->right = create_node(*value, *final_depth, leaf, nullptr, nullptr, fod_elements);
						return leaf->right;
					}
				}else{
					if(leaf->is_null){
						// disjunction node
						//move_disjunction_to_regular_node(leaf);
						leaf->is_null = false;
						++this->number_of_non_null_values;
						//this->non_null_nodes_by_depth[leaf->depth + 1].push_back(leaf);
					}
					leaf->value = *value;
					return leaf;
				}
			}else {
				if(leaf->left){
					++(*depth);
					return insert(key, fod_elements, final_depth, value, depth, leaf, leaf->left);
				}else{
					// link this node to the position of the fod element corresponding to *final_depth
					leaf->left = create_node(*value, *final_depth, leaf, nullptr, nullptr, fod_elements);
					return leaf->left;
				}
			}
		}

		/////////////////////////////////////////

		set_N_value<T>* find(const boost::dynamic_bitset<>& key, size_t *depth, const size_t *final_depth, node<T> *leaf) const {
			// if leaf->depth > *final_depth, then leaf has other elements than the ones of key
			if(leaf->depth > *final_depth){
				return nullptr;
			}
			if(*depth < leaf->depth){
				// take skipped depths into account
				const boost::dynamic_bitset<>& next_set = this->fod->to_set(leaf->fod_elements);
				while(*depth < leaf->depth){
					if(key[*depth] != next_set[*depth]){
						// if there is a disjunction between key and next_set for the element at *depth in fod
						return nullptr;
					}
					++(*depth);
				}
			}
			//*depth = leaf->depth;
			if(key[*depth]){
				if(*depth != *final_depth){
					if(!leaf->right){
						return nullptr;
					}
					++(*depth);
					return find(key, depth, final_depth, leaf->right);
				}else{
					if(leaf->is_null)
						return nullptr;
					else
						return leaf;
				}
			}else {
				if(!leaf->left){
					return nullptr;
				}
				++(*depth);
				return find(key, depth, final_depth, leaf->left);
			}
		}

		set_N_value<T>* find(const boost::dynamic_bitset<>& key, size_t final_depth) const {
			/*
			if(key.size() <= *depth){
				std::string key_str = "";
				for (int i = 0; i < key.size(); ++i) {
					std::string bit_str = "0";
					if(key[i]){
						bit_str = "1";
					}
					key_str = bit_str + key_str;
				}
				throw "\nCardinality parameter given for key " << key_str << " superior to its actual cardinality\n";
			}
			*/
			size_t depth = 0;
			if(final_depth != 0){
				--final_depth;
				return find(key, &depth, &final_depth, this->root);
			}else{
				if(this->emptyset->is_null)
					return nullptr;
				else
					return this->emptyset;
			}
		}

		/////////////////////////////////////////

		static void add_to_values(std::vector<set_N_value<T>* >& values, node<T>* leaf){
			values.emplace_back(leaf);
		}

		static void add_to_values_by_cardinality(std::unordered_map<size_t, std::vector<set_N_value<T>* > >& values, node<T>* leaf){
			if(values.find(leaf->fod_elements.size()) != values.end()){
				values[leaf->fod_elements.size()].emplace_back(leaf);
			}else{
				values.emplace(leaf->fod_elements.size(), {leaf});
			}
		}

		/*
		 * Complexity O(P), where P = number of nodes assigned by user
		 */
		template<typename return_type>
		void elements(node<T>* leaf, return_type& all_values, std::function<void(return_type&, node<T>*)> add_to_values_func) const {
			/*if(!leaf){
				std::clog << "null";
				return;
			}*/

			if(!leaf->is_null){
				std::clog << leaf->value;
				add_to_values_func(all_values, leaf);
			}else{
				std::clog << "null";
			}


			if(leaf->left){
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
				elements(leaf->left, all_values, add_to_values_func);
			}

			if(leaf->right){
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
				elements(leaf->right, all_values, add_to_values_func);
			}
		}

		/*
		 * Complexity O(P), where P = number of nodes assigned by user
		void elements_by_set_cardinality(node<T>* leaf, std::unordered_map<size_t, std::vector<set_N_value<T>* > >& all_values) const {
			if(!leaf)
				return;

			if(!leaf->is_null){
				all_values[leaf->fod_elements.size()].emplace_back(leaf);
			}
			elements_by_set_cardinality(leaf->left, all_values);
			elements_by_set_cardinality(leaf->right, all_values);
		}*/

		/*
		 * Complexity O(S), where S is the number of subsets assigned by user of the set defined by key
		 * All subsets of key are of depth <= final_depth and don't contain any other elements than the ones of key
		 */
		template<typename return_type>
		void subsets_of(const boost::dynamic_bitset<>& key, const size_t& final_depth,
				return_type& subset_values, size_t depth, node<T> *leaf, const bool& is_key, const bool& strict, std::function<void(return_type&, node<T>*)> add_to_values_func) const {

			// if leaf->depth > *final_depth, then leaf has other elements than the ones of key
			if(leaf->depth > final_depth){
				return;
			}
			if(depth < leaf->depth){
				// take skipped depths into account
				const boost::dynamic_bitset<>& next_set = this->fod->to_set(leaf->fod_elements);
				while(depth < leaf->depth){
					if(!key[depth] && next_set[depth]){
						// if next_set has an element that key doesn't
						return;
					}
					++depth;
				}
			}
			//depth = leaf->depth;
			if(key[depth]){
				if(depth != final_depth){
					// get value only if leaf doesn't correspond to key (only strict subsets)
					if(!leaf->is_null){
						add_to_values_func(subset_values, leaf);
					}
					++depth;
					if(leaf->left){
						subsets_of(key, final_depth, subset_values, depth, leaf->left, false, strict, add_to_values_func);
					}
					if(leaf->right){
						subsets_of(key, final_depth, subset_values, depth, leaf->right, is_key, strict, add_to_values_func);
					}
				}else if(!(strict && is_key)){//this->fod->is_subset_of(this->fod->to_set(leaf->fod_elements), key)){
					// get value if you don't only want strict subsets
					// OR if this set is a *strict* subset of key
					if(!leaf->is_null){
						add_to_values_func(subset_values, leaf);
					}
				}
			}else {
				if(leaf->left){
					++depth;
					subsets_of(key, final_depth, subset_values, depth, leaf->left, is_key, strict, add_to_values_func);
				}
			}
		}
		/*
		std::vector<T&> non_null_values_of_subsets_of(boost::dynamic_bitset<> key, int final_depth){
			std::vector<T&> subset_values{this->emptyset->value};
			int depth = 0;
			if(final_depth >= 0)
				non_null_values_of_subsets_of(key, &final_depth, subset_values, depth, this->root);
			return subset_values;
		}
		*/
		std::vector<set_N_value<T>* > subsets_of(const boost::dynamic_bitset<>& key, size_t final_depth, const bool& strict) const {
			std::vector<set_N_value<T>* > subset_values;
			subset_values.reserve(this->size());

			if(final_depth > 0){
				if(!this->emptyset->is_null){
					subset_values.emplace_back(this->emptyset);
				}

				--final_depth;
				size_t depth = 0;
				subsets_of<std::vector<set_N_value<T>* > >(
						key, final_depth, subset_values, depth, this->root, true, strict, add_to_values);
			}else{
				if(!this->emptyset->is_null && !strict){
					subset_values.emplace_back(this->emptyset);
				}
			}
			return subset_values;
		}

		std::unordered_map<size_t, std::vector<set_N_value<T>* > > subsets_of_by_cardinality(
				const boost::dynamic_bitset<>& key, size_t final_depth, const bool& strict) const {

			std::unordered_map<size_t, std::vector<set_N_value<T>* > > subset_values;

			//std::vector<std::vector<set_N_value<T>* > > subset_values(this->fod->size()+1);
			// very rough approximation of the number of slots that should be reserved for each cardinality
			//float alloc_per_cardinality = ceil((float) this->size() / (this->fod->size()-1));
			/*
			for (int i = 1; i < key.count(); ++i) {
				// always only one emptyset and one FOD,
				// and no set of cardinality greater than key.count() could be a subset of it
				subset_values[i].reserve(alloc_per_cardinality);
			}*/

			if(!this->emptyset->is_null){
				//subset_values[0].emplace_back(this->emptyset);
				subset_values.emplace(0, {this->emptyset});
			}
			if(final_depth > 0){
				--final_depth;
				size_t depth = 0;
				subsets_of<std::unordered_map<size_t, std::vector<set_N_value<T>* > > >(
						key, final_depth, subset_values, depth, this->root, true, strict, add_to_values_by_cardinality);
			}
			return subset_values;
		}

		/////////////////////////////////////////

		/*
		 * Complexity O(S), where S is the number of supersets assigned by user of the set defined by key
		 * All supersets of key are of depth >= final_depth and contain every element of key
		 */
		template<typename return_type>
		void supersets_of(const boost::dynamic_bitset<>& key, const size_t *final_depth,
				return_type& superset_values, size_t depth, node<T> *leaf, const bool strict, std::function<void(return_type&, node<T>*)> add_to_values_func) const {

			if(depth < leaf->depth){
				// take skipped depths into account
				const boost::dynamic_bitset<>& next_set = this->fod->to_set(leaf->fod_elements);
				while(depth < leaf->depth){
					if(key[depth] && !next_set[depth]){
						// if key has an element that next_set doesn't
						return;
					}
					++depth;
				}
			}
			//depth = leaf->depth;
			if(depth < *final_depth){
				if(key[depth++]){
					// if key has an element at depth, so its supersets
					if(leaf->right){
						supersets_of(key, final_depth, superset_values, depth, leaf->right, strict, add_to_values_func);
					}
				}else{
					// if key has no element at depth, its supersets have one or not
					if(leaf->left){
						supersets_of(key, final_depth, superset_values, depth, leaf->left, strict, add_to_values_func);
					}
					if(leaf->right){
						supersets_of(key, final_depth, superset_values, depth, leaf->right, strict, add_to_values_func);
					}
				}
			}else if(depth == *final_depth){

				if(!strict || this->fod->is_superset_of(this->fod->to_set(leaf->fod_elements), key)){
					// get value if you don't only want strict supersets
					// OR if this set is a *strict* superset of key
					if(!leaf->is_null){
						add_to_values_func(superset_values, leaf);
					}
				}
				++depth;
				if(leaf->right){
					supersets_of(key, final_depth, superset_values, depth, leaf->right, strict, add_to_values_func);
				}
			}else{
				// if depth > *final_depth, then we are visiting strict supersets of key
				if(!leaf->is_null){
					add_to_values_func(superset_values, leaf);
				}

				++depth;
				// Now, no matter what additional elements they have, they are all supersets of key
				if(leaf->left){
					supersets_of(key, final_depth, superset_values, depth, leaf->left, strict, add_to_values_func);
				}
				if(leaf->right){
					supersets_of(key, final_depth, superset_values, depth, leaf->right, strict, add_to_values_func);
				}
			}
		}

		std::vector<set_N_value<T>* > supersets_of(const boost::dynamic_bitset<>& key, size_t final_depth, const bool strict) const {

			if(final_depth == 0){
				std::vector<set_N_value<T>* > superset_values = this->elements();
				if(strict){
					superset_values.erase(superset_values.begin());
				}
				return superset_values;
			}else{
				--final_depth;
				size_t depth = 0;
				std::vector<set_N_value<T>* > superset_values;
				superset_values.reserve(this->size());

				supersets_of<std::vector<set_N_value<T>* > >(
						key, &final_depth, superset_values, depth, this->root, strict, add_to_values);

				return superset_values;
			}
		}

		std::unordered_map<size_t, std::vector<set_N_value<T>* > > supersets_of_by_cardinality(const boost::dynamic_bitset<>& key, size_t final_depth, const bool strict) const {

			if(final_depth == 0){
				std::unordered_map<size_t, std::vector<set_N_value<T>* > > superset_values = this->elements_by_set_cardinality();
				if(strict){
					superset_values.erase(0);
				}
				return superset_values;
			}else{
				--final_depth;
				size_t depth = 0;
				std::vector<std::vector<set_N_value<T>* > > superset_values(this->fod->size()+1);
				// very rough approximation of the number of slots that should be reserved for each cardinality
				float alloc_per_cardinality = ceil((float) this->node_pool.get_bock_size() / (this->fod->size()-1));
				for (size_t i = key.count()+1; i < this->fod->size(); ++i) {
					// always only one emptyset and one FOD,
					// and no set of cardinality less than key.count() could be a superset of it
					superset_values[i].reserve(alloc_per_cardinality);
				}

				supersets_of<std::vector<std::vector<set_N_value<T>* > > >(
						key, &final_depth, superset_values, depth, this->root, strict, add_to_values_by_cardinality);

				return superset_values;
			}
		}
	};
}	// namespace ow_bft

#endif // OW_BFT_POWERSET_BTREE_HPP
