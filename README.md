# efficient-DST

C++ framework allowing one to perform most transformations required in Dempster-Shafer Theory (DST), such as computing the commonality function q from a mass function m (and back) or from a conjunctive decomposition (and back), or the implicability/belief function b (and plausibility function p) from a mass function (and back) or from a disjunctive decomposition (and back).
This framework is based on my works around complexity reduction for Möbius transforms, including my Efficient Möbius Transformations (EMT).

This project implements each set as a std::bitset (a template building a contiguous array of N bits). The core of this library is `mobius_inversion.hpp`, which contains the templates implementing the *Efficient Möbius Transform (EMT)* algorithms published in this [paper](https://arxiv.org/pdf/2107.07359.pdf) (extension of the [paper](https://link.springer.com/chapter/10.1007/978-3-030-35514-2_29) published in the proceedings of Scalable Uncertainty Management: 13th International Conference, SUM 2019).
These algorithms feature complexities always better than the Fast Möbius Transform (FMT) one (but with overheads in real time). More specifically, they have time complexities between O(F) and O(N.2^N), where F is the number of focal sets in the considered powerset function and N is the number of elements in the frame of discernment (FOD). In fact, you are no longer limited by N for the computation of Möbius transforms such as the weight functions in the conjunctive and disjunctive decompositions. These algorithms exploit both the information in focal sets and the structure of distributive lattices.
Their spatial complexities are between O(F) and O(2^N).

But, beware of the number of focal sets resulting from multiple fusions. Using techniques approximating the mass function based on its focal sets is always wise.

In addition, you will find in the `figures` directory some plots comparing the **Fast Möbius Transform (FMT)** vs **naive (direct) computation with focal sets** vs **semilattice EMT** vs **lattice EMT**, under different circumstances: number of random focal sets (support), special focal sets structures (almost-bayesian, almost-consonant), powerset size (2^N, with N ranging from 10 to 1600). Left vertical axis represents the elapsed real time to compute the Möbius transform of some random powerset function definition. For the EMT algorithms, this includes the computation of extra support elements. The right vertical axis represents the total number of these *scheme* support elements (focal sets + extra support elements). Note: the number of scheme support elements for *direct* and *semilattice* EMT always overlap because they rely on the same extra elements (focal points), while the *lattice* EMT requires the smallest lattice containing these focal points (necessarily bigger). The *FMT* scheme always relies on the whole powerset. Thus, it always has the worst space complexity, but not always the worst time complexity (it can even have the best one). It depends on the characteristics of the given powerset function. 
For example, if you have a powerset function defined on a powerset of size 2^1600 with only 10 random focal sets, then the FMT is not even tractable, while it only take a fraction of second for the direct and semilattice schemes (direct is best). On the contrary, if you have 500 random focal sets in a powerset of size 2^17, the FMT will be the fastest. In between, if you have 100 random focal sets in a powerset of size 2^N, with N > 20, the semilattice EMT scheme will be the fastest and the most space-efficient.

The SVG file `class_diagram.svg` contains a light partial class diagram, giving an overview of the structure of this framework.

# Usage

- Copy the include repository into your C++ project,
- Link it to your compilation process,
- See src/demo.hpp to get an example of usage.

# Future updates

Currently, I am re-implementing the main fusion functions of DST (working branch `fusion_rules`), but I've been quite busy since the end of my PhD. So, feel free to contribute!

If there is some extra functionality that you would like to be implemented, please let me know.
