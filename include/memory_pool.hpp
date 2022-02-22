#ifndef MEMORY_POOL_HPP
#define MEMORY_POOL_HPP

#include <unordered_map>
#include <vector>
#include <iostream>

#include <macros.hpp>

namespace efficient_DST{

	/*class slot_indices {
	public:
		size_t block;
		size_t index;

		slot_indices(size_t _block, size_t _index) : block(_block), index(_index)
		{}
	};*/

	/*
	 * Stores contiguous elements with vectors, but split them in a vector to avoid reallocations (which invalidate pointers)
	 */
	template <class T>
	class memory_pool{
	protected:
		//const size_t tolerance = 10;
		size_t block_size;
		std::vector<std::vector<T> > slots;
		std::vector<std::vector<T*> > reusable_slots;
		std::unordered_map<T*, size_t> block_by_slot;

	public:

		memory_pool(const size_t _block_size) : block_size(_block_size){
			add_new_block();
		}

		~memory_pool(){
			DEBUG(std::clog << "memory_pool destroyed." << std::endl;);
			this->clear();
		}

		/////////////////////////////////////////

		const size_t& get_block_size() const {
			return this->block_size;
		}

		size_t number_of_occupants() const {
			size_t size = 0;
			for (size_t b = 0; b < this->slots.size(); ++b){
				size += this->slots[b].size() - this->reusable_slots[b].size();
			}
			return size;
		}

		template<typename... Ts>
		T* emplace(Ts... args){
			// if there is at least one slot returned to memory_pool

			// get the first contiguous block with a free reusable slot
			size_t b = 0;
			while(b < this->reusable_slots.size() && this->reusable_slots[b].size() == 0){
				++b;
			}
			if(b == this->slots.size()){
				// if there is none, create a new slot

				--b;
				if(this->slots[b].size() == this->block_size){
					// if the last block is full, create a new one
					add_new_block();
					++b;
				}
				size_t end_index = this->slots[b].size();
				// initialize slot with a new instance with parameters
				this->slots[b].emplace_back(args...);
				// create entry in block_by_slot with its block index in this->slots
				// (for reverse querying)
				this->block_by_slot.emplace(&this->slots[b][end_index], b);

				// return a pointer to this free space in slots
				return &this->slots[b][end_index];

			}else{
				// otherwise, simply reuse the last one of them (LIFO)
				T* object_p = this->reusable_slots[b].back();
				// remove it from reusable slots
				this->reusable_slots[b].pop_back();
				*object_p = T(args...);
				// return a pointer to this free space in slots
				return object_p;
			}
		}

		void erase(T* slot){

			// if this address corresponds to the slot at the end of this->slots
			if(slot == &this->slots.back().back()){
				// erase that slot
				this->block_by_slot.erase(slot);
				this->slots.back().pop_back();

				// Also, if there is more than one block
				// AND the last block is entirely free
				// (AND the before last block has at least this->tolerance free_indices)
				if(this->slots.size() > 1 && this->slots.back().size() == this->reusable_slots.back().size()
						//&& this->slots[this->slots.size()-2].size() <= this->block_size - this->tolerance
						){
					// then remove this last block
					for(size_t i = 0; i < this->slots.back().size(); ++i){
						this->block_by_slot.erase(&this->slots.back()[i]);
					}
					this->slots.pop_back();
					this->reusable_slots.pop_back();
				}
			}else{
				// otherwise, reuse that slot
				try {
					this->reusable_slots[this->block_by_slot.at(slot)].emplace_back(slot);
				}catch(const std::out_of_range& oor){
					// if this is not one of the address used in this memory pool
					std::cerr << "[memory_pool::erase] Out of Range error: " << oor.what() << '\n';
				}
			}
			/*
			try {
				// get indices by querying index_by_slot with a slot address
				slot_indices sa = this->index_by_slot.at(slot);
				// if this address corresponds to the slot at the end of this->slots
				if(sa.block == this->slots.size()-1 && sa.index == this->slots[sa.block].size()-1){
					// erase that slot
					this->index_by_slot.erase(slot);
					this->slots.back().pop_back();

					// Also, if there is more than one block
					// AND the last block is entirely free
					// (AND the before last block has at least this->tolerance free_indices)
					if(this->slots.size() > 1 && this->slots.back().size() == 0
							//&& this->slots[this->slots.size()-2].size() <= this->block_size - this->tolerance
							){
						// then remove this last block
						this->slots.pop_back();
						this->reusable_slots.pop_back();
					}
				}else{
					// otherwise, reuse that slot
					this->reusable_slots[sa.block].emplace_back(slot);
				}
			}catch(const std::out_of_range& oor){
				// if this is not one of the address used in this memory pool
				std::cerr << "[memory_pool::erase] Out of Range error: " << oor.what() << '\n';
			}*/
		}

		void clear(){
			this->slots.clear();
			this->reusable_slots.clear();
			this->block_by_slot.clear();
			this->add_new_block();
		}

	protected:

		void add_new_block(){
			this->slots.emplace_back();
			this->slots.back().reserve(this->block_size);
			this->reusable_slots.emplace_back();
			this->reusable_slots.back().reserve(this->block_size);
		}
	};

}	// namespace efficient_DST

#endif // MEMORY_POOL_HPP
