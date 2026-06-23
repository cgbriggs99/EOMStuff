/*
 * theory-input.hpp
 *
 *  Created on: Jun 23, 2026
 *      Author: connor
 */

#ifndef DIAGRAMS_THEORY_INPUT_HPP_
#define DIAGRAMS_THEORY_INPUT_HPP_

#include <string>
#include <utility>

namespace diagram {
class Theory {
public:

    enum TheoryType {
        FULL = 0, // This is full treatment of the given excitations.
        PERTURB = 1, // Treat all excitations but the highest at the last perturbative order that only includes the considered excitations. Treat the highest excitations at the first non-vanishing perturbative order.
        // This is the CCn formalism. Singles are always treated exactly. CC2 treats doubles at 1st order. CC3 treats doubles at 3rd order, triples at 2nd order. CC4 treats doubles and triples at 4th order,
        // quadruples at 3rd order.
        PARENTHESIS = 2, // This is full treatment of all excitations but the highest one. The highest is treated perturbatively through the Wigner 2n + 1 rule. This is the parenthesis formalism.
        // When the order is 2, this gives MP2. That is, CC(D) is essentially MP2.
        BRACKET = 3, // The difference between Parenthesis and Bracket is in the energy expression. The amplitudes are the same.
        LINEAR = 4 // Only consider terms with one T amplitude.
    };

    Theory(std::string const &theory_str);

    // 0 is hf. 1 is ccs if full, otherwise, it doesn't make sense. 2 is ccsd/mp2/cc2 depending on the type. 3 and up is what would be expected.
    // 3 would be ccsdt/ccsd(t)/cc3, 4 would be ccsdtq/ccsdt(q)/cc4, etc.
    constexpr Theory(unsigned int level, TheoryType type, bool skip_singles = false) : level_ { level }, type_ { type }, skip_singles_ {
            skip_singles } {

        // Part of the CCn theory is using the T1-transformed formalism.
        // We can't skip singles.
        if (skip_singles_ && type_ == PERTURB) {
            skip_singles_ = false;
        }

    }

    constexpr unsigned int get_level() const {
        return level_;
    }
    constexpr TheoryType get_type() const {
        return type_;
    }

    constexpr bool do_skip_singles() const {
        return skip_singles_;
    }

    std::string method_string() const;

    bool operator==(Theory const &other) const;

    // This is mostly subjective. The orderings used for this are:
    // If levelA < levelB, then true.
    // If levelA == levelB, then consider the approximation.
    // Skip singles < no skips
    // Linear < Bracket < Parenthesis (by design) < CCn (due to being an iterative method)
    bool operator<=(Theory const &other) const;

    inline bool operator>(Theory const &other) const {
        return !(*this <= other);
    }

    inline bool operator!=(Theory const &other) const {
        return !(*this == other);
    }

    inline bool operator>=(Theory const &other) const {
        return other <= *this;
    }

    inline bool operator<(Theory const &other) const {
        return other > *this;
    }

private:
    unsigned int level_;

    TheoryType type_;

    bool skip_singles_;
};
}

#endif /* DIAGRAMS_THEORY_INPUT_HPP_ */
