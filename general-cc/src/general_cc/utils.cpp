/*
 * utils.cpp
 *
 *  Created on: Jun 24, 2026
 *      Author: connor
 */

#include <general_cc/utils.hpp>
#include <unordered_map>
#include <tuple>

namespace psi::general_cc {
namespace detail {

size_t pair_hash::operator()(std::pair<unsigned short, unsigned short> const &key) const {
    size_t n, k;

    n = std::get < 0 > (key);
    k = std::get < 1 > (key);

    return (n * (n + 1)) / 2 + k;
}

size_t binomial_coefficient_driver(unsigned short n, unsigned short in_k, binomial_storage_type *dynamic_list) {

    unsigned short k = in_k;

    if (n < k) {
        return 0;
    }

    if (n == 0) {
        return 0;
    }

    if (k > n - k) {
        k = n - k;
    }

    // Base case: (n C 0) = 1 and (n C n) = 1
    if (k == 0 || k == n) {
        return 1;
    }

    // Optimized case: (n C 1) = n
    if (k == 1) {
        return n;
    }

    // Dynamic step: see if we've already handled this case.
    auto has_pair = dynamic_list->find(std::pair<unsigned short, unsigned short> { n, k });

    if (has_pair != dynamic_list->cend()) {
        return std::get<1>(*has_pair);
    }

    size_t retval = binomial_coefficient_driver(n - 1, k, dynamic_list) + binomial_coefficient_driver(n - 1, k - 1, dynamic_list);

    (*dynamic_list)[std::pair<unsigned short, unsigned short> { n, k }] = retval;

    return retval;
}
}
}
