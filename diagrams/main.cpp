/*
 * main.cpp
 *
 *  Created on: Jun 3, 2026
 *      Author: connor
 */

#include "diagram.hpp"
#include <stdexcept>

#include <cstdio>

using namespace std;
using namespace diagram;

int main(void) {
	int theory;

	std::printf("Enter the type of theory (0 for CC(n), 1 for CC(n)(n-1), 2 "
				"for CCn): ");

	std::scanf("%d", &theory);

	std::printf("Enter a level of theory: ");

	int level = 0;

	std::scanf("%d", &level);

	int t1_trans;

	std::printf("T1 transform? ");

	std::scanf("%d", &t1_trans);

	for (int i = 0; i <= level; i++) {
		std::list<Diagram> terms;
		switch (theory) {
		case 0:
			terms = std::move(compute_cc_residual(i, level, t1_trans, true));
			break;
		case 1:
			break;
		case 2:
			terms = std::move(diagram::compute_ccn_residual(i, level, t1_trans, true));
			break;
		default:
			throw std::runtime_error("Invalid theory value!");
		}

		if (i == 0) {
			std::printf("\nECC = ");
		} else {
			std::printf("\nR%d = ", i);
		}

		if (terms.size() == 0) {
			std::printf("0\n");
		} else {
			bool plus = false;
			for (auto const &term : terms) {
				if (plus) {
					std::printf(" + ");
				} else {
					plus = true;
				}
				std::printf("%s", term.to_string().c_str());
			}
		}
	}
	std::printf("\n");

	return 0;
}
