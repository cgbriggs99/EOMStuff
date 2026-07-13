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

TEST_CASE("Product Tables") {
    auto[mults, prods] = generate_product_table(5, 2, 2);

    SECTION("Check accuracy") {
        einsums::Tensor < ptrdiff_t, 2 > expected_prods { "Expected products", 10, 10 };
        einsums::Tensor<signed char, 2> expected_mults { "Expected mults", 10, 10 };

        expected_mults.zero();
        expected_prods.zero();

        std::vector < ptrdiff_t > prods_data { 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 0, 0, 1, 0, 0, 4, 0,
                0, 0, 0, 0, 0, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 2, 0, 4, 0, 0, 0, 0, 0, 0, 0, 3, 4, 0, 0, 0, 0, 1, 0, 2, 3, 0, 0, 0, 0, 0, 1,
                0, 2, 0, 4, 0, 0, 0, 0, 1, 0, 0, 3, 4, 0, 0, 0, 0, 0, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0 };
        std::vector<signed char> mults_data { 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, -1, 0, 0, -1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0,
                0, 1, 0, 0, 0, 0, -1, -1, 0, 0, -1, 0, 0, 0, 0, 1, 0, -1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0,
                -1, 0, -1, 0, 1, 0, 0, 0, 0, 1, 0, 0, -1, -1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0 };

        expected_prods.vector_data() = prods_data;
        expected_mults.vector_data() = mults_data;

        println (expected_prods);
        println(prods);

        println(expected_mults);
        println(mults);

        REQUIRE(mults.dim(0) == expected_mults.dim(0));
        REQUIRE(mults.dim(1) == expected_mults.dim(1));
        REQUIRE(prods.dim(0) == expected_prods.dim(0));
        REQUIRE(prods.dim(1) == expected_prods.dim(1));

        for (int i = 0; i < mults.dim(0); i++) {
            for (int j = 0; j < mults.dim(1); j++) {
                REQUIRE((prods(i, j) == expected_prods(i, j) || expected_mults(i, j) == 0));
                REQUIRE(mults(i, j) == expected_mults(i, j));
            }
        }
    }

    SECTION("Check read/write") {
        write_product_table(mults, prods, "test");

        auto[mults2, prods2] = read_product_table("test");

        REQUIRE(mults2.dim(0) == mults.dim(0));
        REQUIRE(mults2.dim(1) == mults.dim(1));
        REQUIRE(prods2.dim(0) == prods.dim(0));
        REQUIRE(prods2.dim(1) == prods.dim(1));

        for (int i = 0; i < prods.dim(0); i++) {
            for (int j = 0; j < prods.dim(0); j++) {
                REQUIRE(mults2(i, j) == mults(i, j));
                REQUIRE(prods2(i, j) == prods(i, j));
            }
        }
    }
}
