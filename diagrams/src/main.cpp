/*
 * main.cpp
 *
 *  Created on: Jun 3, 2026
 *      Author: connor
 */

#include <diagrams/diagram.hpp>
#include <diagrams/theory-input.hpp>

#include <stdexcept>

#include <cstdio>
#include <string>
#include <cctype>

using namespace std;
using namespace diagram;

int main(void) {
    std::printf("Enter a level of theory:\n");

    char str[256] = { 0 };

    bool front_whitespace = true;

    int i;

    for (i = 0; i < 255; i++) {

        auto ch = getc(stdin);

        if (ch == '\n' || ch == 0) {
            break;
        }

        if (front_whitespace && isspace(ch)) {
            continue;
        } else {
            front_whitespace = false;
        }
        str[i] = ch;
        if (ch == 0) {
            break;
        }
    }

    // Strip whitespace from the back.

    for (int j = i - 1; j >= 0; j--) {
        if (isspace(str[j])) {
            str[j] = 0;
        } else {
            break;
        }
    }

    // Parse string.
    diagram::Theory theory { std::string(str) };

    int t1_trans;

    if (!theory.do_skip_singles() && theory.get_type() != Theory::LINEAR) {

        std::printf("T1 transform? ");

        std::scanf("%d", &t1_trans);
    } else {
        t1_trans = true;
    }

    printf("Excitation level: %d\n", theory.get_level());

    for (int i = 0; i <= theory.get_level(); i++) {
        if (i == 1 && theory.do_skip_singles()) {
            continue;
        }
        std::list<Diagram> terms;
        switch (theory.get_type()) {
        case Theory::FULL:
            terms = std::move(compute_cc_residual(i, theory.get_level(), t1_trans || theory.do_skip_singles(), !t1_trans));
            break;
        case Theory::PARENTHESIS:
            terms = std::move(compute_ccp_residual(i, theory.get_level(), t1_trans || theory.do_skip_singles(), !t1_trans));
            break;
        case Theory::PERTURB:
            terms = std::move(compute_ccn_residual(i, theory.get_level(), t1_trans || theory.do_skip_singles(), !t1_trans));
            break;
        case Theory::LINEAR:
            terms = std::move(compute_lcc_residual(i, theory.get_level(), theory.do_skip_singles(), true));
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
