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

#ifndef EFFICIENT_DST_MACROS_HPP
#define EFFICIENT_DST_MACROS_HPP

//#define DEBUG_BUILD
//#define DEBUG_BUILD_TREE

#ifdef DEBUG_BUILD
#  define DEBUG(x) do x while (0)
#else
#  define DEBUG(x) do {} while (0)
#endif

#ifdef DEBUG_BUILD_TREE
#  define DEBUG_TREE(x) do x while (0)
#else
#  define DEBUG_TREE(x) do {} while (0)
#endif

#endif
