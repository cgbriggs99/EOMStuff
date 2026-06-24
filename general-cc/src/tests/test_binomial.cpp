/*
 * test_binomial.cpp
 *
 *  Created on: Jun 24, 2026
 *      Author: connor
 */

#include <catch2/catch_test_macros.hpp>
#include <general_cc/utils.hpp>

using namespace psi::general_cc;

TEST_CASE("Binomial") {
    REQUIRE(binomial_coefficient(39, 33) == 3262623UL);
    REQUIRE(binomial_coefficient(64, 14) == 47855699958816UL);
    REQUIRE(binomial_coefficient(67, 15) == 345780890878896UL);
}
