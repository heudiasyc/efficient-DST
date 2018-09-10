# efficient-DST

This is a C++ implementation of the Dempster-Shafer Theory (DST) that allows for successful computations even with a substantial amount of elements N in the frame of discernment (FOD) by reducing time complexities from <img src="https://latex.codecogs.com/gif.latex?O(2^{2N})" /> to <img src="https://latex.codecogs.com/gif.latex?O(F^2)" />, where <img src="https://latex.codecogs.com/gif.latex?F " /> is the number of focal elements in mass function, and spatial complexities from <img src="https://latex.codecogs.com/gif.latex?O(2^{N})" /> to <img src="https://latex.codecogs.com/gif.latex?O(F)" />.

The core of this implementation uses the class dynamic\_bitset (a dynamic contiguous sequence of bits) of the library BOOST to represent sets, along with a special dynamic binary tree of my design allowing for the search of a set with a complexity linear in N in the worst case and the storage of all representations (mass, belief, commonality, conjunctive decomposition, etc) with a spatial complexity equal or close to <img src="https://latex.codecogs.com/gif.latex?O(F)" /> in average.

Concerning the methods used to fuse and compute representations other than mass functions, there were few adaptations to be made, except for the case of the conjunctive/disjunctive decomposition.

To the best of my knowledge, there was no article or implementation featuring a method for the computation of these decompositions in time <img src="https://latex.codecogs.com/gif.latex?O(F^2)" /> and space <img src="https://latex.codecogs.com/gif.latex?O(F)" />. So I had to create one. 