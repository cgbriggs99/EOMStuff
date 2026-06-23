/*
 * theory-input.cpp
 *
 *  Created on: Jun 22, 2026
 *      Author: connor
 */

#include <diagrams/theory-input.hpp>
#include <diagrams/utils.hpp>

#include <string>
#include <stdexcept>
#include <tuple>
#include <cstdio>
#include <sstream>

namespace diagram {

static std::pair<unsigned int, Theory::TheoryType> process_parentheses(std::string const &method, int start_pos) {
    int open_parens = 0;
    bool is_open = false;
    char number[32] = { 0 };
    int number_pos = 0;

    int first_number = 0, second_number = 0;

    for (int i = 0; i + start_pos < method.size(); i++) {
        if (method[i + start_pos] == '(') {
            if (is_open) {
                throw std::invalid_argument(
                        "Mismatched parentheses! Value of the method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number, and non-perturbative theories can be prefixed with L.");
            }
            if (open_parens >= 2) {
                throw std::invalid_argument(
                        "Too many parentheses found! Value of the method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number, and non-perturbative theories can be prefixed with L.");
            }
            is_open = true;
            open_parens++;
        } else if (method[i + start_pos] == ')') {
            if (!is_open) {
                throw std::invalid_argument(
                        "Mismatched parentheses! Value of the method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number, and non-perturbative theories can be prefixed with L.");
            }

            is_open = false;

            if (number_pos == 0) {
                throw std::invalid_argument(
                        "No number given! Value of the method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number, and non-perturbative theories can be prefixed with L.");
            }

            if (open_parens == 1) {
                std::sscanf(number, "%d", &first_number);
            } else {
                std::sscanf(number, "%d", &second_number);
            }
            for (int j = 0; j <= number_pos; j++) {
                number[j] = 0;
            }
            number_pos = 0;
        } else if (isdigit(method[i + start_pos])) {
            if (number_pos == 32) {
                throw std::invalid_argument(
                        "Number argument too long! Value of the method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number, and non-perturbative theories can be prefixed with L.");
            }
            number[number_pos] = method[i + start_pos];
            number_pos++;
        } else {
            throw std::invalid_argument(
                    "Obtained a non-digit character when expecting a digit! Value of the method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                            "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number, and non-perturbative theories can be prefixed with L.");
        }
    }

    if (first_number <= 0) {
        throw std::invalid_argument(
                "No valid numbers found! Value of the method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                        "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number, and non-perturbative theories can be prefixed with L.");
    }

    if (second_number <= 0) {
        return {first_number, Theory::FULL};
    } else if (second_number == first_number + 1) {
        return {second_number, Theory::PARENTHESIS};
    } else {
        throw std::invalid_argument(
                "Incompatible values! Value of the method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                        "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number, and non-perturbative theories can be prefixed with L.");
    }
}

static char letters[] = "SDTQPH";

Theory::Theory(std::string const &theory_str) {
    char number[32] = { 0 };

    // Search through the string. Start by seeing if the string is one of HF or MP2.
    std::string upper = theory_str;

    to_upper(&upper);

    if (upper == "HF" || upper == "SCF") {
        level_ = 0;
        type_ = FULL;
    } else if (upper == "MP2") {
        level_ = 2;
        type_ = PARENTHESIS;
        skip_singles_ = true;
    } else {
        int num_c = 0;

        int curr_pos = 0;
        type_ = FULL;
        if (upper[0] == 'L') {
            type_ = LINEAR;
        } else if (upper[0] != 'C') {
            throw std::invalid_argument(
                    "Value of the method string is not acceptable! It should be of one of the following forms: CCD, CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                            "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number, and non-perturbative theories can be prefixed with L.");
        }
        for (curr_pos = 0; curr_pos < upper.size() && num_c < 2; curr_pos++) {
            if (upper[curr_pos] == 'C') {
                num_c++;
            }
        }

        if (num_c < 2 || curr_pos == upper.size()) {
            throw std::invalid_argument(
                    "Value of the method string is not acceptable! It should be of one of the following forms: CCD, CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                            "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number, and non-perturbative theories can be prefixed with L.");
        }

        if (upper[curr_pos] == '(') {
            auto[level, type] = process_parentheses(upper, curr_pos);
            level_ = level;
            type_ = type;
        } else if (isdigit(upper[curr_pos])) {
            for (int i = 0; i < 31 && curr_pos + i < upper.size(); i++) {
                number[i] = upper[i + curr_pos];
            }

            std::sscanf(number, "%d", &level_);
            if (type_ == LINEAR) {
                throw std::invalid_argument(
                        "Cannot do linear perturbation! Value of the method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number, and non-perturbative theories can be prefixed with L.");
            }
            type_ = PERTURB;
        } else {
            int i, skip_single = 0;
            int found_open = 0;
            int found_closed = 0;
            for (i = 0; i < 6 - skip_single + found_open + found_closed && curr_pos + i < upper.size(); i++) {
                if (upper[i + curr_pos] == '(') {
                    found_open++;
                    if (i + curr_pos < upper.size() - 2) {
                        if (upper[i + 1 + curr_pos] != letters[i]) {
                            throw std::invalid_argument(
                                    "Unknown letter! Value of the method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                            "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number, and non-perturbative theories can be prefixed with L.");
                        }
                        if (upper[i + 2 + curr_pos] != ')') {
                            throw std::invalid_argument(
                                    "Malformed string! Value of the method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                            "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number, and non-perturbative theories can be prefixed with L.");
                        }
                    }
                } else if (upper[i + curr_pos] == ')') {
                    found_closed++;
                    if (i + curr_pos != upper.size() - 1) {
                        throw std::invalid_argument(
                                "Malformed string! Value of the method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                        "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number, and non-perturbative theories can be prefixed with L.");
                    }
                } else if (upper[i + curr_pos] != letters[i - found_open - found_closed + skip_single]) {
                    if (i == 0 && upper[i + curr_pos] == 'D') {
                        skip_single = 1;
                    } else {
                        throw std::invalid_argument(
                                "Unknown letter! Value of the method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                        "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number, and non-perturbative theories can be prefixed with L.");
                    }
                }
                level_ = i - found_open - found_closed + 1 + skip_single;
            }

            if (skip_single == 1) {
                skip_singles_ = true;
            }

            if (i == 0) {
                throw std::invalid_argument(
                        "No level specified after CC! Value of the method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number, and non-perturbative theories can be prefixed with L.");
            }

            if (found_open == 1 && found_closed == 1) {
                if (type_ == LINEAR) {
                    throw std::invalid_argument(
                            "Cannot do linear perturbation! Value of the method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                    "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number, and non-perturbative theories can be prefixed with L.");
                }
                type_ = PARENTHESIS;
            } else if (found_open == 0 && found_closed == 0) {
                if (type_ != LINEAR) {
                    type_ = FULL;
                }
            } else {
                throw std::invalid_argument(
                        "Malformed string! Value of the method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number, and non-perturbative theories can be prefixed with L.");
            }
        }
    }
}

std::string Theory::method_string() const {
    std::ostringstream out;

    // First, check for special cases.
    if (level_ == 0) {
        return "SCF";
    }

    if (level_ == 2 && (type_ == PARENTHESIS || type_ == PERTURB || type_ == BRACKET) && skip_singles_) {
        return "MP2";
    }

    // Next, check to see if we have a perturbative method.
    if (type_ == PERTURB) {
        out << "CC" << level_;

        return out.str();
    }

    // Next, check to see if we have a small method.
    if (level_ <= 6) {
        if (type_ == LINEAR) {
            out << 'L';
        }
        out << "CC";

        int highest_full = level_;

        if (type_ == PARENTHESIS || type_ == BRACKET) {
            highest_full--;
        }

        // Now, we put letters in the stream.
        for (int i = 0; i < highest_full; i++) {
            if (i == 0 && skip_singles_) {
                continue;
            }
            out << letters[i];
        }

        if (type_ == PARENTHESIS) {
            out << '(' << letters[level_ - 1] << ')';
        } else if (type_ == BRACKET) {
            out << '[' << letters[level_ - 1] << ']';
        }

        return out.str();
    }

    // Finally, we output the default case.
    if (type_ == LINEAR) {
        out << 'L';
    }

    out << "CC";

    switch (type_) {
    case FULL:
    case LINEAR:
        out << '(' << level_ << ')';
        break;
    case PARENTHESIS:
        out << '(' << (level_ - 1) << ")(" << level_ << ')';
        break;
    case BRACKET:
        out << '(' << (level_ - 1) << ")[" << level_ << ']';
        break;
    default:
        throw std::logic_error("Unexpected enum value!");
    }
    if (skip_singles_) {
        out << "(no-S)";
    }
    return out.str();
}

bool Theory::operator==(Theory const &other) const {
    if (level_ == 0 && other.get_level() == 0) {
        return true;
    }

    if (level_ == other.get_level()) {
        // For these, different theories can give MP2.
        // We'll consider these equal.
        if (level_ == 2) {
            bool is_mp2 = skip_singles_ && (type_ == PARENTHESIS || type_ == BRACKET);
            bool other_is_mp2 = other.do_skip_singles() && (other.get_type() == PARENTHESIS || other.get_type() == BRACKET);

            if (is_mp2 && other_is_mp2) {
                return true;
            }
        }

        // At levels less than 2, the different types of perturbations become the same.
        if (level_ <= 2) {
            if ((type_ == PARENTHESIS || type_ == BRACKET) && (other.get_type() == PARENTHESIS || other.get_type() == BRACKET)) {
                return skip_singles_ == other.do_skip_singles();
            }
        }

        return type_ == other.get_type() && skip_singles_ == other.do_skip_singles();
    }

    return false;
}

bool Theory::operator<=(Theory const &other) const {
    if (level_ > other.get_level()) {
        return false;
    }

    if (level_ == other.get_level()) {
        // For these, different theories can give MP2.
        // We'll consider these equal.
        if (level_ == 2) {
            bool is_mp2 = skip_singles_ && (type_ == PARENTHESIS || type_ == BRACKET);
            bool other_is_mp2 = other.do_skip_singles() && (other.get_type() == PARENTHESIS || other.get_type() == BRACKET);

            if (is_mp2 && other_is_mp2) {
                return true;
            }
        }

        // At levels less than 2, the different types of perturbations become the same.
        if (level_ <= 2) {
            if ((type_ == PARENTHESIS || type_ == BRACKET) && (other.get_type() == PARENTHESIS || other.get_type() == BRACKET)) {
                return skip_singles_ || !other.do_skip_singles();
            }
        }
        if (skip_singles_ && !other.do_skip_singles()) {
            return true;
        }
        if (!skip_singles_ && other.do_skip_singles()) {
            return false;
        }

        return type_ >= other.get_type(); // I know it's greater or equal. It's because the enum order is backwards.
    }

    return true;
}

}
