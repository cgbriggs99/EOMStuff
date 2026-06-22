/*
 * utils.cpp
 *
 *  Created on: Jun 22, 2026
 *      Author: connor
 */

#include <utils.hpp>

#include <string>
#include <cctype>

namespace psi::general_dlpno {
void to_upper(std::string *str) {
    for (auto &ch : *str) {
        ch = toupper(ch);
    }
}
}

