/*
 * diagram.hpp
 *
 *  Created on: Jun 3, 2026
 *      Author: connor
 */

#ifndef DIAGRAM_HPP_
#define DIAGRAM_HPP_

#include <list>
#include <string>
#include <tuple>
#include <variant>
namespace diagram {

// Based on the work of Kallay and Surjan.
class Diagram {
public:
	using triplet = std::tuple<int, int, int>;
	constexpr Diagram(int interaction, triplet const &t1,
					  triplet const &t2 = {0, 0, 0},
					  triplet const &t3 = {0, 0, 0},
					  triplet const &t4 = {0, 0, 0}) :
		interaction{interaction},
		triple_1{t1},
		triple_2{t2},
		triple_3{t3},
		triple_4{t4} {}

	constexpr int get_interaction() const { return interaction; }

	constexpr triplet get_t1() const { return triple_1; }
	constexpr triplet get_t2() const { return triple_2; }
	constexpr triplet get_t3() const { return triple_3; }
	constexpr triplet get_t4() const { return triple_4; }

	static std::string triple_string(triplet const &triple, int start_i,
									 int start_a, int start_m, int start_e);

	std::string to_string() const;

private:
	int interaction;
	triplet triple_1, triple_2, triple_3, triple_4;
};

std::list<Diagram> compute_cc_residual(int residual_exc, int max_exc,
									   bool t1_transform, bool canonical);

} // namespace diagram

#endif /* DIAGRAM_HPP_ */
