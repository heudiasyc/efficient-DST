# efficient-DST

C++ framework allowing one to perform most transformations required in Dempster-Shafer Theory (DST) (if there is some extra functionality that you would like to be implemented, please let me know).
This framework is based on my works around complexity reduction for Möbius transforms, including my Efficient Möbius Transformations (EMT).
These works led to algorithms featuring complexities always better than the Fast Möbius Transform (FMT) one. More specifically, they are all between <img src="https://latex.codecogs.com/gif.latex?O(F)" />, where <img src="https://latex.codecogs.com/gif.latex?F " /> is the number of focal sets in the considered mass function,
and <img src="https://latex.codecogs.com/gif.latex?O(N.2^{N})" />, exploiting both the information in focal sets and the structure of the Boolean lattice.
Spatial complexities are between <img src="https://latex.codecogs.com/gif.latex?O(F)" /> and <img src="https://latex.codecogs.com/gif.latex?O(2^{N})" />. See [] for more details.
In fact, you are no longer limited by the number of elements in the frame of discernment (FOD) <img src="https://latex.codecogs.com/gif.latex?N" />. 
But, beware of the number of focal sets resulting from multiple fusions. Using techniques approximating the mass function based on its focal sets is always wise.


The core of this implementation uses the class dynamic\_bitset (a dynamic contiguous sequence of bits) of the library BOOST to represent sets, along with a special dynamic binary tree of my design.

# Usage

- Copy the include repository into your C++ project,
- Link it to your compilation process,
- See src/main.cpp to get an example of usage.