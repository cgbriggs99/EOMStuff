/*
 * diagram.cpp
 *
 *  Created on: Jun 3, 2026
 *      Author: connor
 */

#include "diagram.hpp"
#include <sstream>

namespace diagram {

std::string Diagram::triple_string(triplet const &triple, int start_i,
								   int start_a, int start_m, int start_e) {

	std::ostringstream out;
	if (std::get<0>(triple) > 0) {
		out << "t(";

		int internal_lines = std::get<1>(triple);

		int internal_particles = std::get<2>(triple);

		int internal_holes = internal_lines - internal_particles;

		int external_particles = std::get<0>(triple) - internal_particles;

		int external_holes = std::get<0>(triple) - internal_holes;

		bool begin = true;

		for (int i = 0; i < external_holes; i++) {
			if (begin) {
				out << "i" << i + start_i;
				begin = false;
			} else {
				out << ", i" << i + start_i;
			}
		}

		for (int i = 0; i < internal_holes; i++) {
			if (begin) {
				out << "m" << i + start_m;
				begin = false;
			} else {
				out << ", m" << i + start_m;
			}
		}

		for (int a = 0; a < external_particles; a++) {
			out << ", a" << a + start_a;
		}

		for (int a = 0; a < internal_particles; a++) {
			out << ", e" << a + start_e;
		}

		out << ")";
	}

	return out.str();
}

std::string Diagram::to_string() const {

	std::ostringstream out;

	if (interaction < 0 || interaction > 12) {
		throw std::runtime_error(
			"Interaction id is outside of the expected range!");
	}

	int curr_i_ind = 1;
	int curr_a_ind = 1;
	int curr_m_ind = 1;
	int curr_e_ind = 1;

	out << triple_string(triple_1, curr_i_ind, curr_a_ind, curr_m_ind,
						 curr_e_ind);

	int internal_lines = std::get<1>(triple_1);

	int internal_particles = std::get<2>(triple_1);

	int internal_holes = internal_lines - internal_particles;

	int external_particles = std::get<0>(triple_1) - internal_particles;

	int external_holes = std::get<0>(triple_1) - internal_holes;

	curr_i_ind += external_holes;
	curr_a_ind += external_particles;
	curr_m_ind += internal_holes;
	curr_i_ind += internal_particles;

	out << triple_string(triple_2, curr_i_ind, curr_a_ind, curr_m_ind,
						 curr_e_ind);

	internal_lines = std::get<1>(triple_2);

	internal_particles = std::get<2>(triple_2);

	internal_holes = internal_lines - internal_particles;

	external_particles = std::get<0>(triple_2) - internal_particles;

	external_holes = std::get<0>(triple_2) - internal_holes;

	curr_i_ind += external_holes;
	curr_a_ind += external_particles;
	curr_m_ind += internal_holes;
	curr_i_ind += internal_particles;

	out << triple_string(triple_3, curr_i_ind, curr_a_ind, curr_m_ind,
						 curr_e_ind);

	internal_lines = std::get<1>(triple_3);

	internal_particles = std::get<2>(triple_3);

	internal_holes = internal_lines - internal_particles;

	external_particles = std::get<0>(triple_3) - internal_particles;

	external_holes = std::get<0>(triple_3) - internal_holes;

	curr_i_ind += external_holes;
	curr_a_ind += external_particles;
	curr_m_ind += internal_holes;
	curr_i_ind += internal_particles;

	out << triple_string(triple_4, curr_i_ind, curr_a_ind, curr_m_ind,
						 curr_e_ind);
	internal_lines = std::get<1>(triple_4);

	internal_particles = std::get<2>(triple_4);

	internal_holes = internal_lines - internal_particles;

	external_particles = std::get<0>(triple_4) - internal_particles;

	external_holes = std::get<0>(triple_4) - internal_holes;

	curr_i_ind += external_holes;
	curr_a_ind += external_particles;
	curr_m_ind = 1;
	curr_e_ind = 1;

	switch (interaction) {
	case 0:
		out << "f(a" << curr_a_ind << ", e1)";
		break;
	case 1:
		out << "f(m1, i" << curr_i_ind << ")";
		break;
	case 2:
		out << "f(m1, e1)";
		break;
	case 3:
		out << "<a" << curr_a_ind << " a" << curr_a_ind + 1 << "||e1 e2>";
		break;
	case 4:
		out << "<m1 m2||i" << curr_i_ind << " i" << curr_i_ind + 1 << ">";
		break;
	case 5:
		out << "<m1 a" << curr_a_ind << "||i" << curr_i_ind << " e1>";
		break;
	case 6:
		out << "<m1 a" << curr_a_ind << "||e1 e2>";
		break;
	case 7:
		out << "<m1 m2||e1 i" << curr_i_ind << ">";
		break;
	case 8:
		out << "<a" << curr_a_ind << " a" << curr_a_ind + 1 << "||i"
			<< curr_i_ind << " m1>";
		break;
	case 9:
		out << "<a" << curr_a_ind << " m1||i" << curr_i_ind << " i"
			<< curr_i_ind + 1 << ">";
		break;
	case 10:
		out << "<m1 m2||e1 e2>";
		break;
	case 11:
		out << "f(a" << curr_a_ind << ", i" << curr_i_ind << ")";
		break;
	case 12:
		out << "<a" << curr_a_ind << " a" << curr_a_ind + 1 << "||i"
			<< curr_i_ind << " i" << curr_i_ind + 1 << ">";
		break;
	default:
		throw std::runtime_error(
			"Interaction id is outside of the expected range!");
	}

	return out.str();
}

// static int exc_level[] = {0, 0, -1, 0, 0, 0, -1, -1, 1, 1, -2, 1, 2};
// static int connections[] = {1, 1, 2, 2, 2, 2, 3, 3, 1, 1, 4, 0, 0};
// static int particle_lines[] = {1, 0, 1, 2, 0, 1, 2, 1, 1, 0, 2, 0, 0};
// static int hole_lines[] = {0, 1, 1, 0, 2, 1, 1, 2, 0, 1, 2, 0, 0};

static std::tuple<int, int, int> zero_tuple{0, 0, 0};

static void assign_Fab_terms(int residual, int max_excitation,
							 std::list<Diagram> *out, bool t1_transform) {
	if (residual >= ((t1_transform) ? 2 : 1) && residual <= max_excitation) {
		out->emplace_back(0, std::make_tuple(residual, 1, 1), zero_tuple,
						  zero_tuple, zero_tuple);
	}
}

static void assign_Fij_terms(int residual, int max_excitation,
							 std::list<Diagram> *out, bool t1_transform) {
	if (residual >= ((t1_transform) ? 2 : 1) && residual <= max_excitation) {
		out->emplace_back(1, std::make_tuple(residual, 1, 0), zero_tuple,
						  zero_tuple, zero_tuple);
	}
}

static void assign_Fia_t_terms(int residual, int max_excitation,
							   std::list<Diagram> *out, bool t1_transform) {
	if (residual + 1 <= max_excitation && residual >= 0 &&
		!(residual == 0 && t1_transform)) {
		out->emplace_back(2, std::make_tuple(residual + 1, 2, 1), zero_tuple,
						  zero_tuple, zero_tuple);
	}
}

static void assign_Fia_t2_terms(int residual, int max_excitation,
								std::list<Diagram> *out, bool t1_transform) {
	for (int x = std::max(residual - max_excitation, (t1_transform) ? 2 : 1);
		 x <= max_excitation; x++) {
		int y = residual - x + 1;

		if (y > max_excitation) {
			continue;
		}

		if (y < ((t1_transform) ? 2 : 1)) {
			continue;
		}

		if (x >= y) {
			out->emplace_back(2, std::make_tuple(x, 1, 1),
							  std::make_tuple(y, 1, 0), zero_tuple, zero_tuple);
		}

		if (x > y) {
			out->emplace_back(2, std::make_tuple(x, 1, 0),
							  std::make_tuple(y, 1, 1), zero_tuple, zero_tuple);
		}
	}
}

static void assign_Fia_terms(int residual, int max_excitation,
							 std::list<Diagram> *out, bool t1_transform) {

	assign_Fia_t_terms(residual, max_excitation, out, t1_transform);

	assign_Fia_t2_terms(residual, max_excitation, out, t1_transform);
}

static void assign_Wabcd_t_terms(int residual, int max_excitation,
								 std::list<Diagram> *out, bool t1_transform) {
	if (max_excitation >= 2 && residual >= 2 && residual <= max_excitation) {
		out->emplace_back(3, std::make_tuple(residual, 2, 2), zero_tuple,
						  zero_tuple, zero_tuple);
	}
}

static void assign_Wabcd_t2_terms(int residual, int max_excitation,
								  std::list<Diagram> *out, bool t1_transform) {
	for (int x = std::max(residual - max_excitation, (t1_transform) ? 2 : 1);
		 x <= max_excitation; x++) {
		int y = residual - x;

		if (y > max_excitation) {
			continue;
		}

		if (y < ((t1_transform) ? 2 : 1)) {
			continue;
		}

		if (x >= y) {
			out->emplace_back(3, std::make_tuple(x, 1, 1),
							  std::make_tuple(y, 1, 0), zero_tuple, zero_tuple);
		}
	}
}

static void assign_Wabcd_terms(int residual, int max_excitation,
							   std::list<Diagram> *out, bool t1_transform) {

	assign_Wabcd_t_terms(residual, max_excitation, out, t1_transform);
	assign_Wabcd_t2_terms(residual, max_excitation, out, t1_transform);
}

static void assign_Wijkl_t_terms(int residual, int max_excitation,
								 std::list<Diagram> *out, bool t1_transform) {
	if (max_excitation >= 2 && residual >= 2 && residual <= max_excitation) {
		out->emplace_back(4, std::make_tuple(residual, 2, 0), zero_tuple,
						  zero_tuple, zero_tuple);
	}
}

static void assign_Wijkl_t2_terms(int residual, int max_excitation,
								  std::list<Diagram> *out, bool t1_transform) {
	for (int x = std::max(residual - max_excitation, (t1_transform) ? 2 : 1);
		 x <= max_excitation; x++) {
		int y = residual - x;

		if (y > max_excitation) {
			continue;
		}

		if (y < ((t1_transform) ? 2 : 1)) {
			continue;
		}

		if (x >= y) {
			out->emplace_back(4, std::make_tuple(x, 1, 0),
							  std::make_tuple(y, 1, 0), zero_tuple, zero_tuple);
		}
	}
}

static void assign_Wijkl_terms(int residual, int max_excitation,
							   std::list<Diagram> *out, bool t1_transform) {
	assign_Wijkl_t_terms(residual, max_excitation, out, t1_transform);

	assign_Wijkl_t2_terms(residual, max_excitation, out, t1_transform);
}

static void assign_Waijb_t_terms(int residual, int max_excitation,
								 std::list<Diagram> *out, bool t1_transform) {
	if (residual > 0 && !(t1_transform && residual == 1) &&
		residual <= max_excitation) {
		out->emplace_back(5, std::make_tuple(residual, 2, 1), zero_tuple,
						  zero_tuple, zero_tuple);
	}
}

static void assign_Waijb_t2_terms(int residual, int max_excitation,
								  std::list<Diagram> *out, bool t1_transform) {
	for (int x = std::max(residual - max_excitation, (t1_transform) ? 2 : 1);
		 x <= max_excitation; x++) {
		int y = residual - x;

		if (y > max_excitation) {
			continue;
		}

		if (y < ((t1_transform) ? 2 : 1)) {
			continue;
		}

		if (x >= y) {
			out->emplace_back(5, std::make_tuple(x, 1, 1),
							  std::make_tuple(y, 1, 0), zero_tuple, zero_tuple);
		}

		if (x > y) {
			out->emplace_back(5, std::make_tuple(x, 1, 0),
							  std::make_tuple(y, 1, 1), zero_tuple, zero_tuple);
		}
	}
}

static void assign_Waijb_terms(int residual, int max_excitation,
							   std::list<Diagram> *out, bool t1_transform) {
	assign_Waijb_t_terms(residual, max_excitation, out, t1_transform);
	assign_Waijb_t2_terms(residual, max_excitation, out, t1_transform);
}

static void assign_Wiabc_t_terms(int residual, int max_excitation,
								 std::list<Diagram> *out, bool t1_transform) {
	if (max_excitation >= 2 && residual >= 1 &&
		residual + 1 <= max_excitation) {
		out->emplace_back(6, std::make_tuple(residual + 1, 3, 2), zero_tuple,
						  zero_tuple, zero_tuple);
	}
}

static void assign_Wiabc_t2_terms(int residual, int max_excitation,
								  std::list<Diagram> *out, bool t1_transform) {

	for (int x = std::max(residual - max_excitation, (t1_transform) ? 2 : 1);
		 x <= max_excitation; x++) {
		int y = residual - x + 1;

		if (y > max_excitation) {
			continue;
		}

		if (y < ((t1_transform) ? 2 : 1)) {
			continue;
		}

		if (x >= y) {
			if (x >= 2) {
				out->emplace_back(6, std::make_tuple(x, 2, 2),
								  std::make_tuple(y, 1, 0), zero_tuple,
								  zero_tuple);
			}
			out->emplace_back(6, std::make_tuple(x, 2, 1),
							  std::make_tuple(y, 1, 1), zero_tuple, zero_tuple);
		}

		if (x > y) {
			out->emplace_back(6, std::make_tuple(x, 1, 1),
							  std::make_tuple(y, 2, 1), zero_tuple, zero_tuple);
			if (y >= 2) {
				out->emplace_back(6, std::make_tuple(x, 1, 0),
								  std::make_tuple(y, 2, 2), zero_tuple,
								  zero_tuple);
			}
		}
	}
}

static void assign_Wiabc_t3_terms(int residual, int max_excitation,
								  std::list<Diagram> *out, bool t1_transform) {
	for (int z = t1_transform ? 2 : 1;
		 z <= std::min(max_excitation, residual - ((t1_transform) ? 4 : 2));
		 z++) {

		for (int y = t1_transform ? 2 : 1;
			 y <= std::min(z, residual - z - ((t1_transform) ? 2 : 1)); y++) {
			int x = residual - z - y + 1;

			if (x > max_excitation) {
				continue;
			}

			if (x < ((t1_transform) ? 2 : 1)) {
				continue;
			}

			if (x >= y && y >= z) {
				out->emplace_back(6, std::make_tuple(x, 1, 1),
								  std::make_tuple(y, 1, 1),
								  std::make_tuple(z, 1, 0), zero_tuple);
			}

			if (x >= y && y > z) {
				out->emplace_back(6, std::make_tuple(x, 1, 1),
								  std::make_tuple(y, 1, 0),
								  std::make_tuple(z, 1, 1), zero_tuple);
			}

			if (x > y && y >= z) {
				out->emplace_back(6, std::make_tuple(x, 1, 0),
								  std::make_tuple(y, 1, 1),
								  std::make_tuple(z, 1, 1), zero_tuple);
			}
		}
	}
}

static void assign_Wiabc_terms(int residual, int max_excitation,
							   std::list<Diagram> *out, bool t1_transform) {
	assign_Wiabc_t_terms(residual, max_excitation, out, t1_transform);

	assign_Wiabc_t2_terms(residual, max_excitation, out, t1_transform);

	assign_Wiabc_t3_terms(residual, max_excitation, out, t1_transform);
}

static void assign_Wijak_t_terms(int residual, int max_excitation,
								 std::list<Diagram> *out, bool t1_transform) {
	if (max_excitation >= 2 && residual >= 1 &&
		residual + 1 <= max_excitation) {
		out->emplace_back(7, std::make_tuple(residual + 1, 3, 1), zero_tuple,
						  zero_tuple, zero_tuple);
	}
}

static void assign_Wijak_t2_terms(int residual, int max_excitation,
								  std::list<Diagram> *out, bool t1_transform) {
	for (int x = std::max(residual - max_excitation, (t1_transform) ? 2 : 1);
		 x <= max_excitation; x++) {
		int y = residual - x + 1;

		if (y > max_excitation) {
			continue;
		}

		if (y < ((t1_transform) ? 2 : 1)) {
			continue;
		}

		if (x >= y) {
			out->emplace_back(7, std::make_tuple(x, 2, 1),
							  std::make_tuple(y, 1, 0), zero_tuple, zero_tuple);
			out->emplace_back(7, std::make_tuple(x, 2, 0),
							  std::make_tuple(y, 1, 1), zero_tuple, zero_tuple);
		}

		if (x > y) {
			out->emplace_back(7, std::make_tuple(x, 1, 1),
							  std::make_tuple(y, 2, 0), zero_tuple, zero_tuple);
			out->emplace_back(7, std::make_tuple(x, 1, 0),
							  std::make_tuple(y, 2, 1), zero_tuple, zero_tuple);
		}
	}
}

static void assign_Wijak_t3_terms(int residual, int max_excitation,
								  std::list<Diagram> *out, bool t1_transform) {
	for (int z = t1_transform ? 2 : 1;
		 z <= std::min(max_excitation, residual - ((t1_transform) ? 4 : 2));
		 z++) {

		for (int y = t1_transform ? 2 : 1;
			 y <= std::min(z, residual - z - ((t1_transform) ? 2 : 1)); y++) {
			int x = residual - z - y + 1;

			if (x > max_excitation) {
				continue;
			}

			if (x < ((t1_transform) ? 2 : 1)) {
				continue;
			}

			if (x >= y && y >= z) {
				out->emplace_back(7, std::make_tuple(x, 1, 1),
								  std::make_tuple(y, 1, 0),
								  std::make_tuple(z, 1, 0), zero_tuple);
			}

			if (x >= y && y > z) {
				out->emplace_back(7, std::make_tuple(x, 1, 0),
								  std::make_tuple(y, 1, 0),
								  std::make_tuple(z, 1, 1), zero_tuple);
			}

			if (x > y && y >= z) {
				out->emplace_back(7, std::make_tuple(x, 1, 0),
								  std::make_tuple(y, 1, 1),
								  std::make_tuple(z, 1, 0), zero_tuple);
			}
		}
	}
}

static void assign_Wijak_terms(int residual, int max_excitation,
							   std::list<Diagram> *out, bool t1_transform) {
	assign_Wijak_t_terms(residual, max_excitation, out, t1_transform);

	assign_Wijak_t2_terms(residual, max_excitation, out, t1_transform);

	assign_Wijak_t3_terms(residual, max_excitation, out, t1_transform);
}

static void assign_Wabic_terms(int residual, int max_excitation,
							   std::list<Diagram> *out, bool t1_transform) {
	if (residual - 1 <= max_excitation && !(residual == 2 && t1_transform) &&
		residual >= 2) {
		out->emplace_back(8, std::make_tuple(residual - 1, 1, 1), zero_tuple,
						  zero_tuple, zero_tuple);
	}
}

static void assign_Waijk_terms(int residual, int max_excitation,
							   std::list<Diagram> *out, bool t1_transform) {
	if (residual - 1 <= max_excitation && !(residual == 2 && t1_transform) &&
		residual >= 2) {
		out->emplace_back(9, std::make_tuple(residual - 1, 1, 0), zero_tuple,
						  zero_tuple, zero_tuple);
	}
}

static void assign_Wijab_t_terms(int residual, int max_excitation,
								 std::list<Diagram> *out, bool t1_transform) {
	if (residual + 2 <= max_excitation && max_excitation >= 2) {
		out->emplace_back(10, std::make_tuple(residual + 2, 4, 2), zero_tuple,
						  zero_tuple, zero_tuple);
	}
}
static void assign_Wijab_t2_terms(int residual, int max_excitation,
								  std::list<Diagram> *out, bool t1_transform) {

	for (int x = std::max(residual - max_excitation, (t1_transform) ? 2 : 1);
		 x <= max_excitation; x++) {
		int y = residual - x + 2;

		if (y > max_excitation) {
			continue;
		}

		if (y < ((t1_transform) ? 2 : 1)) {
			continue;
		}

		if (x >= y) {
			if (x >= 2) {
				out->emplace_back(10, std::make_tuple(x, 3, 2),
								  std::make_tuple(y, 1, 0), zero_tuple,
								  zero_tuple);
				out->emplace_back(10, std::make_tuple(x, 3, 1),
								  std::make_tuple(y, 1, 1), zero_tuple,
								  zero_tuple);

				out->emplace_back(10, std::make_tuple(x, 2, 2),
								  std::make_tuple(y, 2, 0), zero_tuple,
								  zero_tuple);
			}
			out->emplace_back(10, std::make_tuple(x, 2, 1),
							  std::make_tuple(y, 2, 1), zero_tuple, zero_tuple);
		}

		if (x > y && y >= 2) {
			out->emplace_back(10, std::make_tuple(x, 1, 0),
							  std::make_tuple(y, 3, 2), zero_tuple, zero_tuple);
			out->emplace_back(10, std::make_tuple(x, 1, 1),
							  std::make_tuple(y, 3, 1), zero_tuple, zero_tuple);
			out->emplace_back(10, std::make_tuple(x, 2, 0),
							  std::make_tuple(y, 2, 2), zero_tuple, zero_tuple);
		}
	}
}
static void assign_Wijab_t3_terms(int residual, int max_excitation,
								  std::list<Diagram> *out, bool t1_transform) {

	for (int z = t1_transform ? 2 : 1;
		 z <= std::min(max_excitation, residual - ((t1_transform) ? 4 : 2));
		 z++) {

		for (int y = t1_transform ? 2 : 1;
			 y <= std::min(z, residual - z - ((t1_transform) ? 2 : 1)); y++) {
			int x = residual - z - y + 2;

			if (x > max_excitation) {
				continue;
			}

			if (x < ((t1_transform) ? 2 : 1)) {
				continue;
			}

			if (x >= y && y >= z) {
				if (x >= 2) {
					out->emplace_back(10, std::make_tuple(x, 2, 2),
									  std::make_tuple(y, 1, 0),
									  std::make_tuple(z, 1, 0), zero_tuple);
				}
				out->emplace_back(10, std::make_tuple(x, 2, 1),
								  std::make_tuple(y, 1, 1),
								  std::make_tuple(z, 1, 0), zero_tuple);
				out->emplace_back(10, std::make_tuple(x, 2, 0),
								  std::make_tuple(y, 1, 1),
								  std::make_tuple(z, 1, 1), zero_tuple);
			}

			if (x >= y && y > z) {
				out->emplace_back(10, std::make_tuple(x, 2, 1),
								  std::make_tuple(y, 1, 0),
								  std::make_tuple(z, 1, 1), zero_tuple);
				if (z >= 2) {
					out->emplace_back(10, std::make_tuple(x, 1, 0),
									  std::make_tuple(y, 1, 0),
									  std::make_tuple(z, 2, 2), zero_tuple);
				}
				out->emplace_back(10, std::make_tuple(x, 1, 1),
								  std::make_tuple(y, 1, 0),
								  std::make_tuple(z, 2, 1), zero_tuple);
				out->emplace_back(10, std::make_tuple(x, 1, 1),
								  std::make_tuple(y, 1, 1),
								  std::make_tuple(z, 2, 0), zero_tuple);
			}

			if (x > y && y >= z) {
				if (y >= 2) {
					out->emplace_back(10, std::make_tuple(x, 1, 0),
									  std::make_tuple(y, 2, 2),
									  std::make_tuple(z, 1, 0), zero_tuple);
				}
				out->emplace_back(10, std::make_tuple(x, 1, 1),
								  std::make_tuple(y, 2, 1),
								  std::make_tuple(z, 1, 0), zero_tuple);
				out->emplace_back(10, std::make_tuple(x, 1, 0),
								  std::make_tuple(y, 2, 1),
								  std::make_tuple(z, 1, 1), zero_tuple);
				out->emplace_back(10, std::make_tuple(x, 1, 1),
								  std::make_tuple(y, 2, 0),
								  std::make_tuple(z, 1, 1), zero_tuple);
			}

			if (x > y && y > z) {
				out->emplace_back(10, std::make_tuple(x, 1, 0),
								  std::make_tuple(y, 1, 1),
								  std::make_tuple(z, 2, 1), zero_tuple);
			}
		}
	}
}
static void assign_Wijab_t4_terms(int residual, int max_excitation,
								  std::list<Diagram> *out, bool t1_transform) {

	for (int w = t1_transform ? 2 : 1;
		 w <= std::min(max_excitation, residual - ((t1_transform) ? 6 : 3));
		 w++) {
		for (int z = t1_transform ? 2 : 1;
			 z <=
			 std::min(max_excitation, residual - w - ((t1_transform) ? 4 : 2));
			 z++) {
			for (int y = t1_transform ? 2 : 1;
				 y <= std::min(max_excitation,
							   residual - w - z - ((t1_transform) ? 2 : 1));
				 y++) {
				int x = residual - w - z - y + 2;

				if (x > max_excitation) {
					continue;
				}

				if (x < ((t1_transform) ? 2 : 1)) {
					continue;
				}

				if (x >= y && y >= z && z >= w) {
					out->emplace_back(
						10, std::make_tuple(x, 1, 1), std::make_tuple(y, 1, 1),
						std::make_tuple(z, 1, 0), std::make_tuple(w, 1, 0));
				}

				if (x >= y && y > z && z >= w) {
					out->emplace_back(
						10, std::make_tuple(x, 1, 1), std::make_tuple(y, 1, 0),
						std::make_tuple(z, 1, 1), std::make_tuple(w, 1, 0));
				}
				if (x >= y && y >= z && z > w) {
					out->emplace_back(
						10, std::make_tuple(x, 1, 1), std::make_tuple(y, 1, 0),
						std::make_tuple(z, 1, 0), std::make_tuple(w, 1, 1));
				}
				if (x > y && y >= z && z >= w) {
					out->emplace_back(
						10, std::make_tuple(x, 1, 0), std::make_tuple(y, 1, 1),
						std::make_tuple(z, 1, 1), std::make_tuple(w, 1, 0));
				}
				if (x > y && y >= z && z > w) {
					out->emplace_back(
						10, std::make_tuple(x, 1, 0), std::make_tuple(y, 1, 1),
						std::make_tuple(z, 1, 0), std::make_tuple(w, 1, 1));
				}
				if (x >= y && y > z && z >= w) {
					out->emplace_back(
						10, std::make_tuple(x, 1, 0), std::make_tuple(y, 1, 0),
						std::make_tuple(z, 1, 1), std::make_tuple(w, 1, 1));
				}
			}
		}
	}
}

static void assign_Wijab_terms(int residual, int max_excitation,
							   std::list<Diagram> *out, bool t1_transform) {
	assign_Wijab_t_terms(residual, max_excitation, out, t1_transform);

	assign_Wijab_t2_terms(residual, max_excitation, out, t1_transform);

	assign_Wijab_t3_terms(residual, max_excitation, out, t1_transform);

	assign_Wijab_t4_terms(residual, max_excitation, out, t1_transform);
}

static void assign_Fai_terms(int residual, std::list<Diagram> *out) {
	if (residual == 1) {
		out->emplace_back(11, zero_tuple, zero_tuple, zero_tuple, zero_tuple);
	}
}

static void assign_Wabij_terms(int residual, std::list<Diagram> *out) {
	if (residual == 2) {
		out->emplace_back(12, zero_tuple, zero_tuple, zero_tuple, zero_tuple);
	}
}

std::list<Diagram> compute_cc_residual(int residual, int max_excitation,
									   bool t1_transform, bool canonical) {

	std::list<Diagram> out;

	if (residual < 0 || max_excitation < 1 || residual > max_excitation) {
		throw std::runtime_error("Invalid choice for the residual and "
								 "the excitation level!");
	}

	assign_Fab_terms(residual, max_excitation, &out, t1_transform);

	assign_Fij_terms(residual, max_excitation, &out, t1_transform);

	if (t1_transform || !canonical) {
		assign_Fia_terms(residual, max_excitation, &out, t1_transform);
	}

	assign_Wabcd_terms(residual, max_excitation, &out, t1_transform);

	assign_Wijkl_terms(residual, max_excitation, &out, t1_transform);

	assign_Waijb_terms(residual, max_excitation, &out, t1_transform);

	assign_Wiabc_terms(residual, max_excitation, &out, t1_transform);

	assign_Wijak_terms(residual, max_excitation, &out, t1_transform);

	assign_Wabic_terms(residual, max_excitation, &out, t1_transform);

	assign_Waijk_terms(residual, max_excitation, &out, t1_transform);

	assign_Wijab_terms(residual, max_excitation, &out, t1_transform);

	if (t1_transform || !canonical) {
		assign_Fai_terms(residual, &out);
	}

	assign_Wabij_terms(residual, &out);

	return out;
}

std::list<Diagram> compute_ccn_residual(int residual, int max_excitation,
										bool t1_transform, bool canonical) {

	if (residual < max_excitation) {
		return compute_cc_residual(residual, max_excitation, t1_transform,
								   canonical);
	} else {
		std::list<Diagram> out;

		// When the residual and max excitation match, we need to use the
		// perturbative expansion.
		assign_Fab_terms(residual, max_excitation, &out, true);
		assign_Fij_terms(residual, max_excitation, &out, true);

		if (!canonical) {
			assign_Fia_terms(residual, max_excitation, &out, true);
		}

		assign_Wabcd_terms(residual, max_excitation - 1, &out, t1_transform);

		assign_Wijkl_terms(residual, max_excitation - 1, &out, t1_transform);

		assign_Waijb_terms(residual, max_excitation - 1, &out, t1_transform);

		assign_Wiabc_terms(residual, max_excitation - 2, &out, t1_transform);

		assign_Wijak_terms(residual, max_excitation - 2, &out, t1_transform);

		assign_Wabic_terms(residual, max_excitation - 1, &out, t1_transform);

		assign_Waijk_terms(residual, max_excitation - 1, &out, t1_transform);

		assign_Wijab_terms(residual, max_excitation - 1, &out, t1_transform);

		if (!canonical) {
			assign_Fai_terms(residual, &out);
		}

		assign_Wabij_terms(residual, &out);

		return out;
	}
}

std::list<Diagram> compute_ccp_residual(int residual, int max_excitation,
										bool t1_transform, bool canonical) {
	if (residual == 0) {
		auto out =
			compute_cc_residual(0, max_excitation, t1_transform, canonical);
	}
	if (residual < max_excitation) {
		return compute_cc_residual(residual, max_excitation, t1_transform,
								   canonical);
	} else {
		std::list<Diagram> out;
		
		// CC(n) is correct to n + 1 order. We need to correct out to 2n + 3 order.
		// The CC(n)(n + 1) will have a contribution starting at n order.

		return out;
	}
}

} // namespace diagram
