/*
 * utils.hpp
 *
 *  Created on: Jun 24, 2026
 *      Author: connor
 */

#ifndef GENERAL_CC_INCLUDE_GENERAL_CC_UTILS_HPP_
#define GENERAL_CC_INCLUDE_GENERAL_CC_UTILS_HPP_

#include <unordered_map>
#include <tuple>

namespace psi::general_cc {

namespace detail {

/**
 *  @brief Compute a unique hash for two values where the first value is greater than the second.
 *
 *  We don't have to worry about either value being 0 or 1, or the second value being greater than
 *  or equal to the first, since these base cases are handled by the driver.
 *  The equation I chose for this is to use @$\frac{n(n + 1)}{2} + k@$.
 *
 *  @param key The inputs to the binomial coefficient.
 *  @return A unique value for the pair of inputs.
 */
struct pair_hash {
    size_t operator()(std::pair<unsigned short, unsigned short> const &key) const;
};

using binomial_storage_type = std::unordered_map<std::pair<unsigned short, unsigned short>, size_t, pair_hash>;

/**
 *  @brief Driver function for the binomial coefficient.
 *
 *  This does the actual work for the binomial coefficient function, including handling the
 *  dynamic programming implementation of the algorithm. It is provided for users to potentially
 *  allow for the reuse of a dynamic storage list to avoid reallocations. To call,
 */
size_t binomial_coefficient_driver(unsigned short n, unsigned short k, binomial_storage_type *dynamic_list);
}

inline size_t binomial_coefficient(unsigned short n, unsigned short k) {
    detail::binomial_storage_type dynamic_list;

    return detail::binomial_coefficient_driver(n, k, &dynamic_list);
}

constexpr size_t constexpr_binomial_coefficient(unsigned short n, unsigned short k) {
    if (n == 0) {
        return 0;
    }

    if (k > n) {
        return 0;
    }

    if (k == 0 || k == n) {
        return 1;
    }

    int new_k = k;

    if (k > n - k) {
        new_k = n - k;
    }

    if (new_k == 1) {
        return n;
    }

    if (new_k == 2) {
        return ((size_t) n * ((size_t) n - 1)) / 2;
    }

    return constexpr_binomial_coefficient(n - 1, k) + constexpr_binomial_coefficient(n - 1, k - 1);
}

}

#endif /* GENERAL_CC_INCLUDE_GENERAL_CC_UTILS_HPP_ */
