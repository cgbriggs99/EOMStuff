/*
 * test_parser.cpp
 *
 *  Created on: Jun 22, 2026
 *      Author: connor
 */

#include <general-dlpno.hpp>

#include <cstdio>
#include <cctype>

static void print_type(psi::general_dlpno::Theory::TheoryType type) {
    using namespace psi::general_dlpno;

    switch (type) {
    case Theory::FULL:
        std::printf("Theory type: FULL\n");
        break;
    case Theory::PARENTHESIS:
        std::printf("Theory type: PARENTHESIS\n");
        break;
    case Theory::PERTURB:
        std::printf("Theory type: PERTURB\n");
        break;
    default:
        std::printf("Unknown type!");
    }
}

int main(void) {
    std::printf("Enter a theory string:\n");

    char str[256] = { 0 };

    bool front_whitespace = true;

    int i;

    for (i = 0; i < 255; i++) {

        auto ch = getc(stdin);

        if(ch == '\n' || ch == 0) {
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
    psi::general_dlpno::Theory theory{std::string(str)};

    std::printf("Max excitation: %d\n", theory.get_level());

    print_type(theory.get_type());

    return 0;

}

