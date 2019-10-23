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
