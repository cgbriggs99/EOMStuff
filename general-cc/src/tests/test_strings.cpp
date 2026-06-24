/*
 * test_strings.cpp
 *
 *  Created on: Jun 24, 2026
 *      Author: connor
 */

#include <general_cc/orbit-string.hpp>
#include <catch2/catch_test_macros.hpp>
#include <testing/strings.hpp>

using namespace psi::general_cc;

TEST_CASE("Orbital Strings") {
    SECTION("Create from vector") {
        std::vector<unsigned short> orbitals { 0, 5, 6, 10 };

        OrbitString string { orbitals, 20, OrbitString::ALPHA };

        INFO(string.to_string());

        REQUIRE(string.get_num_indices() == 4);
        REQUIRE(string.get_indices() != nullptr);
        REQUIRE(string.get_max_orbital() == 20);
        REQUIRE(string.get_spin() == OrbitString::ALPHA);
        REQUIRE(string.to_index() == 240);

        for (int i = 0; i < 4; i++) {
            REQUIRE(string.get_index(i) == orbitals[i]);
        }

        REQUIRE_THROWS(string.get_index(4));

        OrbitString string2 = string;

        string2++;

        REQUIRE(string.get_num_indices() == 4);
        REQUIRE(string.get_indices() != nullptr);
        REQUIRE(string.get_max_orbital() == 20);
        REQUIRE(string.get_spin() == OrbitString::ALPHA);
        REQUIRE(string.to_index() == 240);
        REQUIRE(string2.get_num_indices() == 4);
        REQUIRE(string2.get_indices() != nullptr);
        REQUIRE(string2.get_indices() != string.get_indices());
        REQUIRE(string2.get_max_orbital() == 20);
        REQUIRE(string2.get_spin() == OrbitString::ALPHA);
        REQUIRE(string2.to_index() == 241);

        REQUIRE(string.get_index(0) + 1 == string2.get_index(0));
        REQUIRE(string.get_index(1) == string2.get_index(1));
        REQUIRE(string.get_index(2) == string2.get_index(2));
        REQUIRE(string.get_index(3) == string2.get_index(3));

    }

    SECTION("Increment test") {
        OrbitString string1 { 10, 20, OrbitString::ALPHA };

        const size_t combs = string1.number_of_strings();

        for (size_t i = 0; i < combs; i++) {
            OrbitString string2 { i, 10, 20, OrbitString::ALPHA };

            REQUIRE(string1 == string2);
            if (i != combs - 1) {
                REQUIRE(string1.has_next());
                ++string1;
            } else {
                REQUIRE(!string1.has_next());
                REQUIRE_THROWS(++string1);
            }
        }
    }
}
