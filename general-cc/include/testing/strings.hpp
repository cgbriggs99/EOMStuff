/*
 * strings.hpp
 *
 *  Created on: Jun 24, 2026
 *      Author: connor
 */

#ifndef GENERAL_CC_INCLUDE_TESTING_STRINGS_HPP_
#define GENERAL_CC_INCLUDE_TESTING_STRINGS_HPP_

#include <catch2/catch_test_macros.hpp>
#include <general_cc/orbit-string.hpp>
#include <testing/strings.hpp>

namespace Catch {
template<>
struct StringMaker<psi::general_cc::OrbitString> {
    static std::string convert(psi::general_cc::OrbitString const &value) {
        return value.to_string();
    }
};

}

#endif /* GENERAL_CC_INCLUDE_TESTING_STRINGS_HPP_ */
