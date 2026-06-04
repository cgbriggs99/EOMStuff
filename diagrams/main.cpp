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
	std::printf("Enter a level of theory: ");

	int level = 0;

	std::scanf("%d", &level);

	int t1_trans;

	std::printf("T1 transform? ");

	std::scanf("%d", &t1_trans);

	for (int i = 0; i <= level; i++) {
		auto terms = compute_cc_residual(i, level, t1_trans, true);

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

	return 0;
}
