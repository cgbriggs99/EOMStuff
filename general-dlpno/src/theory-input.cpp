/*
 * theory-input.cpp
 *
 *  Created on: Jun 22, 2026
 *      Author: connor
 */

#include <general-dlpno.hpp>
#include <utils.hpp>

#include <string>
#include <stdexcept>
#include <tuple>
#include <cstdio>

namespace psi::general_dlpno {

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
                        "Mismatched parentheses! Value of the DLPNO method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number and any of these can be prefixed with DLPNO-.");
            }
            if (open_parens >= 2) {
                throw std::invalid_argument(
                        "Too many parentheses found! Value of the DLPNO method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number and any of these can be prefixed with DLPNO-.");
            }
            is_open = true;
            open_parens++;
        } else if (method[i + start_pos] == ')') {
            if (!is_open) {
                throw std::invalid_argument(
                        "Mismatched parentheses! Value of the DLPNO method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number and any of these can be prefixed with DLPNO-.");
            }

            is_open = false;

            if (number_pos == 0) {
                throw std::invalid_argument(
                        "No number given! Value of the DLPNO method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number and any of these can be prefixed with DLPNO-.");
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
                        "Number argument too long! Value of the DLPNO method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number and any of these can be prefixed with DLPNO-.");
            }
            number[number_pos] = method[i + start_pos];
            number_pos++;
        } else {
            throw std::invalid_argument(
                    "Obtained a non-digit character when expecting a digit! Value of the DLPNO method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                            "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number and any of these can be prefixed with DLPNO-.");
        }
    }

    if (first_number <= 0) {
        throw std::invalid_argument(
                "No valid numbers found! Value of the DLPNO method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                        "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number and any of these can be prefixed with DLPNO-.");
    }

    if (second_number <= 0) {
        return {first_number, Theory::FULL};
    } else if (second_number == first_number + 1) {
        return {second_number, Theory::PARENTHESIS};
    } else {
        throw std::invalid_argument(
                "Incompatible values! Value of the DLPNO method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                        "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number and any of these can be prefixed with DLPNO-.");
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
    } else if (upper == "MP2" || upper == "DLPNO-MP2" || upper == "LPNO-MP2") {
        level_ = 2;
        type_ = PARENTHESIS;
    } else {
        volatile int num_c = 0;

        volatile int curr_pos = 0;
        for (curr_pos = 0; curr_pos < upper.size() && num_c < 2; curr_pos++) {
            if (upper[curr_pos] == 'C') {
                num_c++;
            }
        }

        if (num_c < 2 || curr_pos == upper.size()) {
            throw std::invalid_argument(
                    "Value of the DLPNO method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                            "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number and any of these can be prefixed with DLPNO-.");
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
            type_ = PERTURB;
        } else {
            volatile int i;
            volatile int found_paren = 0;
            for (i = 0; i < 6 && curr_pos + i < upper.size(); i++) {
                if (upper[i + curr_pos] == '(') {
                    found_paren++;
                    if (i + curr_pos < upper.size() - 2) {
                        if (upper[i + 1 + curr_pos] != letters[i]) {
                            throw std::invalid_argument(
                                    "Unknown letter! Value of the DLPNO method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                            "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number and any of these can be prefixed with DLPNO-.");
                        }
                        if (upper[i + 2 + curr_pos] != ')') {
                            throw std::invalid_argument(
                                    "Malformed string! Value of the DLPNO method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                            "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number and any of these can be prefixed with DLPNO-.");
                        }
                    }
                } else if (upper[i + curr_pos] == ')') {
                    if (i + curr_pos != upper.size() - 1) {
                        throw std::invalid_argument(
                                "Malformed string! Value of the DLPNO method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                        "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number and any of these can be prefixed with DLPNO-.");
                    }
                } else if (upper[i + curr_pos] != letters[i - found_paren]) {
                    throw std::invalid_argument(
                            "Unknown letter! Value of the DLPNO method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                    "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number and any of these can be prefixed with DLPNO-.");
                }
                level_ = i - found_paren;
            }

            if (i == 0) {
                throw std::invalid_argument(
                        "No level specified after CC! Value of the DLPNO method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number and any of these can be prefixed with DLPNO-.");
            }

            if (found_paren == 1) {
                type_ = PARENTHESIS;
            } else if (found_paren == 0) {
                type_ = FULL;
            } else {
                throw std::invalid_argument(
                        "Malformed string! Value of the DLPNO method string is not acceptable! It should be of one of the following forms: CCSD, CCSD(T), CCSDT, CCSDT(Q), CCSDTQ, CCSDTQ(P), "
                                "CCSDTQP, CCSDTQP(H), CCSDTQPH, CC(n), CC(n - 1)(n), or CCn, where n is a number and any of these can be prefixed with DLPNO-.");
            }
        }
    }
}

}
