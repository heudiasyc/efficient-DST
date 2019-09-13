#include <iostream>
#include <string>
#include <vector>

#include <mass.hpp>
#include <commonality.hpp>
#include <belief.hpp>
#include <plausibility.hpp>
#include <pignistic_probability.hpp>
#include <conjunctive_weight.hpp>
#include <disjunctive_weight.hpp>
#include <implicability.hpp>
#include <fod.hpp>
#include <rule_conjunctive.hpp>
#include <rule_conjunctive_cautious.hpp>
#include <rule_disjunctive_bold.hpp>
#include <rule_disjunctive.hpp>
#include <rule_dempster.hpp>
#include <rule_dubois_prade.hpp>


int main(){
    using namespace efficient_DST;

    FOD fod({"e", "f", "g", "h", "i", "j", "l", "o"});

    conjunctive_weight<> w0(fod);

    w0.set_emptyset_value(0.82);
    w0.set_value({"f"}, 803.12195);
    w0.set_value({"f", "g"}, 0.14286);
    w0.set_value({"f", "j"}, 0.02381);
    w0.set_value({"e", "f", "i", "o"}, 0.125);
    w0.set_value({"f", "l"}, 0.03571);

    commonality<> q0(w0);

    std::cout << "\nCommonality values from w0" << std::endl;

    print<>(std::cout, q0.get_definition());

	std::cout << "\n============================================\n";

	mass<> m0(q0);

	std::cout << "\nMass values from w0 " << std::endl;

    print<>(std::cout, m0.get_definition());

	std::cout << "\n============================================\n";

	conjunctive_weight<> w0_back(q0);

	std::cout << "\nConjunctive weights back from q0 " << std::endl;

    print<>(std::cout, w0_back.get_definition());

	std::cout << "\n============================================\n";

    mass<> m(fod);

    m.set_emptyset_value(0.18);
    m.set_value({"f", "g"}, 0.06);
    m.set_value({"f", "j"}, 0.41);
    m.set_value({"e", "f", "i", "o"}, 0.07);
    m.set_value({"f", "l"}, 0.27);
    m.set_fod_value(0.01);

/*
    m.set_emptyset_value(0.18);
    m.set_value({"e", "f", "g", "i", "j", "l", "o"}, 0.06);
    m.set_value({"e", "f", "g", "i"}, 0.38);
    m.set_value({"e", "f", "j", "l", "o"}, 0.1);
    m.set_value({"f", "g", "i", "j", "l"}, 0.27);
    m.set_fod_value(0.01);

	m.set_emptyset_value(0.01);
    m.set_value({"h"}, 0.06);
    m.set_value({"h", "j", "l", "o"}, 0.38);
    m.set_value({"g", "h", "i"}, 0.1);
    m.set_value({"e", "h", "o"}, 0.27);
    m.set_fod_value(0.18);
	*/
    FOD fod_vec({"a", "b", "c", "d"});

    //q_vec corresponds to the mass function {0.01, 0.06, 0, 0, 0.38, 0, 0.05, 0, 0, 0, 0, 0.05, 0.27, 0, 0, 0.18}.
    std::vector<double> q_vec = {1, 0.29, 0.28, 0.23, 0.88, 0.18, 0.23, 0.18, 0.5, 0.23, 0.23, 0.23, 0.45, 0.18, 0.18, 0.18};

    commonality<> q_from_vec(q_vec, fod_vec);

    std::cout << "\nCommonality values " << std::endl;

    print<>(std::cout, q_from_vec.get_definition());

	std::cout << "\n============================================\n";

	mass<> m_back_vec(q_from_vec);

	std::cout << "\nMass values from commonality " << std::endl;

    print<>(std::cout, m_back_vec.get_definition());

	std::cout << "\n============================================\n";

    plausibility<> pl(m);

	std::cout << "\nPlausibility contour from mass " << std::endl;

	std::vector<double> contour_pl = pl.get_contour();
	for (size_t i = 0; i < contour_pl.size(); ++i){
		std::cout << contour_pl[i] << std::endl;
	}

	std::cout << "\n============================================\n";

    pignistic_probability<> bet_p(m);

	std::cout << "\nPignistic probability contour from mass " << std::endl;

	std::vector<double> contour_bet = bet_p.get_contour();
	for (size_t i = 0; i < contour_bet.size(); ++i){
		std::cout << contour_bet[i] << std::endl;
	}

    std::cout << "\n============================================\n";

    commonality<> q(m);

    std::cout << "\nCommonality values " << std::endl;

    print<>(std::cout, q.get_definition());

	std::cout << "\n============================================\n";

	std::cout << q[{"f", "i", "l"}] << std::endl;
	std::cout << q[{"f", "j"}] << std::endl;
	std::cout << q[{"f", "i", "o"}] << std::endl;
	std::cout << q[{"h", "o"}] << std::endl;

	mass<> m_back(q);

	std::cout << "\nMass values from commonality " << std::endl;

    print<>(std::cout, m_back.get_definition());

	std::cout << "\n============================================\n";

	//std::clog << q[{"j", "l"}];

	conjunctive_weight<> w(q);

	std::cout << "\nConjunctive weight values " << std::endl;

    print<>(std::cout, w.get_definition());

	std::cout << "\n============================================\n";

    commonality<> q2(w);

    std::cout << "\nCommonality values from weights " << std::endl;

    print<>(std::cout, q2.get_definition());

	std::cout << "\n============================================\n";

	mass<> m2(q2);

	std::cout << "\nMass values from conjunctive weights " << std::endl;

    print<>(std::cout, m2.get_definition());

	std::cout << "\n============================================\n";


	implicability<> b(m);

	std::cout << "\nImplicability values " << std::endl;

	print<>(std::cout, b.get_definition());

	std::cout << "\n============================================\n";

	mass<> m_back2(b);

	std::cout << "\nMass values from implicability " << std::endl;

    print<>(std::cout, m_back2.get_definition());

	std::cout << "\n============================================\n";

	disjunctive_weight<> v(b);

	std::cout << "\nDisjunctive weight values " << std::endl;

    print<>(std::cout, v.get_definition());

	std::cout << "\n============================================\n";

    implicability<> b2(v);

    std::cout << "\nImplicability values from weights " << std::endl;

    print<>(std::cout, b2.get_definition());

	std::cout << "\n============================================\n";

	mass<> m3(b2);

	std::cout << "\nMass values from disjunctive weights " << std::endl;

    print<>(std::cout, m3.get_definition());

	std::cout << "\n============================================\n";

	conjunctive_weight<> w4(w);

    conjunctive_weight<> w14 = w.apply<rule_Dempster<> >(w4);

	std::cout << "\nMass from conjunctive weight fusion values " << std::endl;

    print<>(std::cout, mass<>(commonality<>(w14)).get_definition());

	std::cout << "\n============================================\n";

	commonality<> q4(q);

	commonality<> q14 = q.apply<rule_conjunctive_cautious<> >(q4);

	std::cout << "\nMass from commonality fusion values " << std::endl;

    print<>(std::cout, mass<>(q14).get_definition());

	std::cout << "\n============================================\n";

	mass<> m4(m);

	mass<> m14 = m.apply<rule_Dubois_Prade<> >(m4);

	std::cout << "\nMass fusion values " << std::endl;

    print<>(std::cout, m14.get_definition());

	std::cout << "\n============================================\n";
}