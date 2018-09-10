// Marek's implementation
#include <boost/bft/mass.hpp>
#include <boost/bft/commonality.hpp>
#include <boost/bft/to_commonality.hpp>
#include <boost/bft/to_mass.hpp>
#include <boost/bft/discounting.hpp>
#include <boost/bft/rule_dubois_prade.hpp>
#include <boost/bft/rule_conjunctive.hpp>

// My implementation
#include <n_discounting.hpp>
#include <mass.hpp>
#include <commonality.hpp>
#include <belief.hpp>
#include <conjunctive_decomposition.hpp>
#include <converter_to_mass_focal_elements.hpp>
#include <implicability.hpp>
#include <fod.hpp>
#include <powerset_btree.hpp>
#include <rule_conjunctive.hpp>
#include <iostream>
#include <fstream>

#include <string>
//#include <matplotlibcpp.h>
#include <bitset>
#include <vector>

#include <sstream>


BOOST_BFT_DEFINE_CLASS(C1);
BOOST_BFT_DEFINE_CLASS(C2);
BOOST_BFT_DEFINE_CLASS(C3);
BOOST_BFT_DEFINE_CLASS(C4);
BOOST_BFT_DEFINE_CLASS(C5);
BOOST_BFT_DEFINE_CLASS(C6);
BOOST_BFT_DEFINE_CLASS(C7);
BOOST_BFT_DEFINE_CLASS(C8);
BOOST_BFT_DEFINE_CLASS(Lambda);

typedef boost::bft::fod<C1, C2> fod2;		// max number of elements with this framework = 10
typedef boost::bft::fod<C1, C2, C3> fod3;
typedef boost::bft::fod<C1, C2, C3, C4> fod4;
/*
powerset(fod) structure = natural order, i.e.
{
	{}, 																							: fod = {}
	{c1}, 																							: fod = {c1}
	{c2}, {c1, c2}, 																				: fod = {c1, c2}
	{c3}, {c1, c3}, {c2, c3}, {c1, c2, c3}, 														: fod = {c1, c2, c3}
	{c4}, {c1, c4}, {c2, c4}, {c1, c2, c4}, {c3, c4}, {c1, c3, c4}, {c2, c3, c4}, {c1, c2, c3, c4},	: fod = {c1, c2, c3, c4}
	...
}
*/

int main(){
/*
	boost::bft::rule_conjunctive rule;

    typedef fod4 fod_t;

    const std::string output_path = "/home/max/eclipse-workspace/DST_experiments/data/output/";

    const boost::bft::mass<fod_t>::container_type ma1 = {0, 0.3, 0, 0.7};       // mass assignment
    // => it would have been smarter to only index the focal elements given their natural order index...
    const boost::bft::mass<fod_t> m1(ma1);                                      // mass definition

    // to do so without explicit indexing of focal elements :
    //std::array<int,4> A = {10,20,30,40};
    //std::array<int,4> B = A; //copy array A into array B
    //double ar[4] = { [1] = 0.3, [3] = 0.7};		// fill an array with 0s, except at 1 and 3 (not implemented with standard gcc compiler)
    const boost::bft::mass<fod_t>::container_type ma2 = {0, 0.3, 0, 0.7};
    const boost::bft::mass<fod_t> m2(ma2);

    std::ofstream myfile;

    boost::bft::mass<fod_t> m12 = m1.apply(rule, m2);     // mass fusion with respect to rule
*/
    using namespace ow_bft;

    FOD fod({"e", "f", "g", "h", "i", "j", "l", "o"});
    //FOD fod({"a", "b", "c"});

    mass<> m(fod);

   // m.set_value({"b"}, 0.3);
    m.set_value({"e", "f", "g", "i", "j", "o"}, 0.05);
    m.set_value({"i", "l", "o"}, 0.48);
    m.set_value({"e", "f", "g", "l", "j"}, 0.27);
    m.set_fod_value(0.2);
/*
    m.set_value({"b", "c", "d", "e"}, 0.3);
    //m.set_value({"a", "b"}, 0.3);
    //m.set_value({"b", "c"}, 0.4);
    m.set_value({"e", "f", "g", "h", "i", "j", "k", "l"}, 0.4);
    m.set_fod_value(0.3);*/

    std::clog << "\n============================================\n";

    commonality<> q(m);

    std::clog << "\nCommonality values " << std::endl;

    print<>(std::clog, q.get_special_elements());

	std::clog << "\n============================================\n";

	belief<> b2(q);

	std::clog << "\nBelief values " << std::endl;

    print<>(std::clog, b2.get_special_elements());

	m.erase_elements_containing_fod_element("i");

	std::clog << "\nMass values " << std::endl;

    print<>(std::clog, m.get_focal_elements());

	std::clog << "\n============================================\n";

	conjunctive_decomposition<> w(q);

	std::clog << "\nWeight values " << std::endl;

    print<>(std::clog, w.get_special_elements());

	std::clog << "\n============================================\n";

	std::clog << "\nBack to commonality values " << std::endl;

	powerset_btree<> q_tree(w.get_FOD(), w.block_size);
	conjunctive_decomposition<>::to_commonality_special_elements(w.get_special_elements(), q_tree);

    print<>(std::clog, q_tree);

	std::clog << "\n============================================\n";

	std::clog << "\nBack to mass values " << std::endl;

	const powerset_btree<>& m_tree = converter_to_mass_focal_elements<commonality<>>(w.precision).convert(q_tree);

    print<>(std::clog, m_tree);
/*
	std::clog << "\n============================================\n";

	std::clog << "\nMass fusion values " << std::endl;

    mass<> m2(m);

    rule_conjunctive rule;
    mass<> m12 = m.apply(rule, m2);     // mass fusion with respect to rule

    print<>(std::clog, m12.get_focal_elements());*/

	/*
	std::clog << std::endl;
    fod.erase("c");
    std::vector<fod_element*> f_e = fod.elements();
    for (size_t i = 0; i < f_e.size(); ++i) {
    	std::clog << f_e[i]->position_in_fod << " ";
    }
    std::clog << std::endl;
    std::clog << "CLEAR" << std::endl;
    values = m.focal_elements->elements();
    std::clog << "powerset.size = " << m.focal_elements->size() << std::endl;
    std::clog << "focal_elements.size = " << values.size() << std::endl;
	for (size_t i = 0; i < values.size(); ++i) {
		std::string k = "";
		if(!values[i]){
			std::clog << "null" << std::endl;
			continue;
		}
		for (size_t j = 0; j < values[i]->fod_elements.size(); ++j) {
			k += " " + patch::to_string<size_t>(values[i]->fod_elements[j]->position_in_fod);
		}
		std::clog << patch::to_string<double>(values[i]->value) << " " << k << std::endl;
	}

    std::clog << m[{}] << "\n";
    std::clog << m[{"a", "b"}] << "\n";
    std::clog << m[{"a", "b", "c"}] << "\n";
    std::clog << m[{"a", "b"}] << "\n";
    std::clog << m[{"a"}] << "\n";
    fod.push_back("b");
    std::clog << m[{"a", "b"}] << "\n";*/


/*
    int n = 10;
    double alpha = 0.1;
    mass<fod_t> m = m1;
    const commonality<fod_t> q1 = to_commonality(m1);
    commonality<fod_t> q;
    n_discounting dn = n_discounting(alpha);
    myfile.open(output_path + "m_discount_n" + patch::to_string(n));

    for (int i = 0; i < 10; ++i) {
    	m = dn.operator ()(i, m1);
    	for (std::size_t j = 0; j < fod_t::powerset_size; ++j) {          // std::size_t = unsigned integer type returned by the sizeof operator
    		myfile << m.values().at(j) << " ";
		}
		myfile << "\n";
	}
    myfile.close();

    */
/*
    const int N = 64;	// 64 is the minimum guaranteed size of unsigned long long int, the largest number type allowed.
    std::bitset<N> omega;
    std::bitset<N> omega_barre;
    //omega.set();
    omega_barre = omega;
    omega_barre.set();
    omega_barre.set(5, false);
    std::clog << omega_barre.to_ullong() << "\n";
    std::vector<std::bitset<N>> oui;
    oui.emplace_back();		// constructs an instance of std::bitset<N> and places it at the end of oui
    						// => simply adds N bits at the end of this vector (1 bit = 1 element of fod)
    oui.emplace_back(omega);
    oui.emplace_back();
    oui[2].set();
    std::clog << oui[0].to_string() << " " << oui[1].to_string() << " " << oui[2].to_string();
    exit(0);
*/
    /*
    std::ofstream myfile;

    mass<fod_t> m12 = m1.apply(rule, m2);     // mass fusion with respect to rule

    int n = 10;
    double alpha = 0.1;
    mass<fod_t> m = m1;
    const commonality<fod_t> q1 = to_commonality(m1);
    commonality<fod_t> q;
    n_discounting dn = n_discounting(alpha);
    myfile.open(output_path + "m_discount_n" + patch::to_string(n));

    for (int i = 0; i < 10; ++i) {
    	m = dn.operator ()(i, m1);
    	for (std::size_t j = 0; j < fod_t::powerset_size; ++j) {          // std::size_t = unsigned integer type returned by the sizeof operator
    		myfile << m.values().at(j) << " ";
		}
		myfile << "\n";
	}
    myfile.close();
	myfile.open (output_path + "q_discount_n" + patch::to_string(n));
	for (int i = 0; i < 10; ++i) {
    	q = dn.operator ()(i, q1);
    	m = to_mass(q);
    	for (std::size_t j = 0; j < fod_t::powerset_size; ++j) {          // std::size_t = unsigned integer type returned by the sizeof operator
    		myfile << q.values().at(j) << " ";
		}
		myfile << "\n";
	}
	myfile.close();

	m = m1;
	myfile.open (output_path + "m_discount" + patch::to_string(n));
	for (std::size_t j = 0; j < fod_t::powerset_size; ++j) {          // std::size_t = unsigned integer type returned by the sizeof operator
		myfile << m.values().at(j) << " ";
	}
	myfile << "\n";
	discounting d = discounting(alpha);
	for (int i = 0; i < 9; ++i) {
    	d.operator ()(m);
    	for (std::size_t j = 0; j < fod_t::powerset_size; ++j) {          // std::size_t = unsigned integer type returned by the sizeof operator
    		myfile << m.values().at(j) << " ";
		}
		myfile << "\n";
	}
	myfile.close();
	*/
/*
 	bool u;
    for (std::size_t i = 0; i < fod_t::powerset_size; ++i) {          // std::size_t = unsigned integer type returned by the sizeof operator
    	u = is_subset_of(i, 12);
    	std::clog << u << " ";
    }
*/
}
