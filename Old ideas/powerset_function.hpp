/*
 * Copyright (C) 2019-2023  Maxime Chaveroche (maxime.chaveroche@gmail.com)
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the CeCILL License, either version 2.1 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * CeCILL License for more details.
 * 
 * You should have received a copy of the CeCILL License
 * along with this program. If not, see <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html>.
 */
 
#ifndef EFFICIENT_DST_POWERSET_FUNCTION_HPP
#define EFFICIENT_DST_POWERSET_FUNCTION_HPP

namespace efficient_DST{

	class powerset_function {
	protected:

		std::vector<powerset_function*> dependent_powerset_functions;

	public:

		void push_back_powerset_function(powerset_function*& p_function){
			this->dependent_powerset_functions.emplace_back(p_function);
		}

		void erase_powerset_function(const powerset_function& p_function){
			auto it = std::find(this->dependent_powerset_functions.begin(), this->dependent_powerset_functions.end(), &p_function);
			this->dependent_powerset_functions.erase(it);
		}

		virtual void erase_elements_containing_fod_element(const size_t position) {
			for (size_t i = 0; i < this->dependent_powerset_functions.size(); ++i) {
				this->dependent_powerset_functions[i]->erase_elements_containing_fod_element(position);
			}
		}

		virtual ~powerset_function()
		{}

	};
}	// namespace efficient_DST

#endif // EFFICIENT_DST_POWERSET_FUNCTION_HPP
