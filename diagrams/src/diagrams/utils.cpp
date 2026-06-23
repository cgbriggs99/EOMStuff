/*
 * utils.cpp
 *
 *  Created on: Jun 22, 2026
 *      Author: connor
 */

#include <diagrams/utils.hpp>

#include <string>
#include <cctype>

namespace diagram {
void to_upper(std::string *str) {
    for (auto &ch : *str) {
        ch = toupper(ch);
    }
}
}

