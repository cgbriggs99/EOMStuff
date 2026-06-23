/*
 * interaction_generators.cpp
 *
 *  Created on: Jun 23, 2026
 *      Author: connor
 */

#include <diagrams/diagram.hpp>

namespace diagram {
static std::tuple<int, int, int> zero_tuple { 0, 0, 0 };

namespace detail {
void assign_Fab_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {
    // F is zero order. We'll assume this even if we are t1-transforming.
    int term_order = residual - 1;

    if (residual == 1) {
        term_order = t1_order;
    }

    if (term_order > total_order) {
        return;
    }

    if (residual >= ((skip_t1)? 2: 1) && residual <= max_excitation) {
        out->emplace_back(Fab, std::make_tuple(residual, 1, 1), zero_tuple, zero_tuple, zero_tuple);
    }
}

void assign_Fij_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {
    int term_order = residual - 1;

    if (residual == 1) {
        term_order = t1_order;
    }

    if (term_order > total_order) {
        return;
    }
    if (residual >= ((skip_t1)? 2: 1) && residual <= max_excitation) {
        out->emplace_back(Fji, std::make_tuple(residual, 1, 0), zero_tuple, zero_tuple, zero_tuple);
    }
}

void assign_Fia_t_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {
    int term_order = residual;

    if (term_order > total_order) {
        return;
    }
    if (residual + 1 <= max_excitation && residual >= 0 && !(residual == 0 && skip_t1)) {
        out->emplace_back(Fia, std::make_tuple(residual + 1, 2, 1), zero_tuple, zero_tuple, zero_tuple);
    }
}

void assign_Fia_t2_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {
    for (int x = std::max(residual - max_excitation, (skip_t1)? 2: 1); x <= max_excitation; x++) {
        int y = residual - x + 1;

        if (y > max_excitation) {
            continue;
        }

        if (y < ((skip_t1)? 2: 1)) {
            continue;
        }

        int term_order = 0;

        if (x == 1) {
            term_order += t1_order;
        } else {
            term_order += x - 1;
        }

        if (y == 1) {
            term_order += t1_order;
        } else {
            term_order += y - 1;
        }

        if (term_order > total_order) {
            continue;
        }

        if (x >= y) {
            out->emplace_back(Fia, std::make_tuple(x, 1, 1), std::make_tuple(y, 1, 0), zero_tuple, zero_tuple);
        }

        if (x > y) {
            out->emplace_back(Fia, std::make_tuple(x, 1, 0), std::make_tuple(y, 1, 1), zero_tuple, zero_tuple);
        }
    }
}

void assign_Fia_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {

    assign_Fia_t_terms(residual, max_excitation, out, skip_t1, total_order, t1_order);

    assign_Fia_t2_terms(residual, max_excitation, out, skip_t1, total_order, t1_order);
}

void assign_Wabcd_t_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {
    // W is a first-order term.
    int term_order = residual;

    if (term_order > total_order) {
        return;
    }
    if (max_excitation >= 2 && residual >= 2 && residual <= max_excitation) {
        out->emplace_back(Wadbc, std::make_tuple(residual, 2, 2), zero_tuple, zero_tuple, zero_tuple);
    }
}

void assign_Wabcd_t2_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {
    for (int x = std::max(residual - max_excitation, (skip_t1)? 2: 1); x <= max_excitation; x++) {
        int y = residual - x;

        if (y > max_excitation) {
            continue;
        }

        if (y < ((skip_t1)? 2: 1)) {
            continue;
        }

        int term_order = 1;

        if (x == 1) {
            term_order += t1_order;
        } else {
            term_order += x - 1;
        }

        if (y == 1) {
            term_order += t1_order;
        } else {
            term_order += y - 1;
        }

        if (term_order > total_order) {
            continue;
        }

        if (x >= y) {
            out->emplace_back(Wadbc, std::make_tuple(x, 1, 1), std::make_tuple(y, 1, 0), zero_tuple, zero_tuple);
        }
    }
}

void assign_Wabcd_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {

    assign_Wabcd_t_terms(residual, max_excitation, out, skip_t1, total_order, t1_order);
    assign_Wabcd_t2_terms(residual, max_excitation, out, skip_t1, total_order, t1_order);
}

void assign_Wijkl_t_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {

    int term_order = residual;

    if (term_order > total_order) {
        return;
    }
    if (max_excitation >= 2 && residual >= 2 && residual <= max_excitation) {
        out->emplace_back(Wjkil, std::make_tuple(residual, 2, 0), zero_tuple, zero_tuple, zero_tuple);
    }
}

void assign_Wijkl_t2_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {
    for (int x = std::max(residual - max_excitation, (skip_t1)? 2: 1); x <= max_excitation; x++) {
        int y = residual - x;

        if (y > max_excitation) {
            continue;
        }

        if (y < ((skip_t1)? 2: 1)) {
            continue;
        }

        int term_order = 1;

        if (x == 1) {
            term_order += t1_order;
        } else {
            term_order += x - 1;
        }

        if (y == 1) {
            term_order += t1_order;
        } else {
            term_order += y - 1;
        }

        if (term_order > total_order) {
            continue;
        }

        if (x >= y) {
            out->emplace_back(Wjkil, std::make_tuple(x, 1, 0), std::make_tuple(y, 1, 0), zero_tuple, zero_tuple);
        }
    }
}

void assign_Wijkl_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {
    assign_Wijkl_t_terms(residual, max_excitation, out, skip_t1, total_order, t1_order);

    assign_Wijkl_t2_terms(residual, max_excitation, out, skip_t1, total_order, t1_order);
}

void assign_Waijb_t_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {
    int term_order = residual;

    if (residual == 1) {
        term_order = t1_order + 1;
    }

    if (term_order > total_order) {
        return;
    }
    if (residual > 0 && !(skip_t1 && residual == 1) && residual <= max_excitation) {
        out->emplace_back(Wjaib, std::make_tuple(residual, 2, 1), zero_tuple, zero_tuple, zero_tuple);
    }
}

void assign_Waijb_t2_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {
    for (int x = std::max(residual - max_excitation, (skip_t1)? 2: 1); x <= max_excitation; x++) {
        int y = residual - x;

        if (y > max_excitation) {
            continue;
        }

        if (y < ((skip_t1)? 2: 1)) {
            continue;
        }

        int term_order = 1;

        if (x == 1) {
            term_order += t1_order;
        } else {
            term_order += x - 1;
        }

        if (y == 1) {
            term_order += t1_order;
        } else {
            term_order += y - 1;
        }

        if (term_order > total_order) {
            continue;
        }

        if (x >= y) {
            out->emplace_back(Wjaib, std::make_tuple(x, 1, 1), std::make_tuple(y, 1, 0), zero_tuple, zero_tuple);
        }

        if (x > y) {
            out->emplace_back(Wjaib, std::make_tuple(x, 1, 0), std::make_tuple(y, 1, 1), zero_tuple, zero_tuple);
        }
    }
}

void assign_Waijb_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {
    assign_Waijb_t_terms(residual, max_excitation, out, skip_t1, total_order, t1_order);
    assign_Waijb_t2_terms(residual, max_excitation, out, skip_t1, total_order, t1_order);
}

void assign_Wiabc_t_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {

    int term_order = residual + 1;

    if (term_order > total_order) {
        return;
    }
    if (max_excitation >= 2 && residual >= 1 && residual + 1 <= max_excitation) {
        out->emplace_back(Wicab, std::make_tuple(residual + 1, 3, 2), zero_tuple, zero_tuple, zero_tuple);
    }
}

void assign_Wiabc_t2_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {

    for (int x = std::max(residual - max_excitation, (skip_t1)? 2: 1); x <= max_excitation; x++) {
        int y = residual - x + 1;

        if (y > max_excitation) {
            continue;
        }

        if (y < ((skip_t1)? 2: 1)) {
            continue;
        }

        int term_order = 1;

        if (x == 1) {
            term_order += t1_order;
        } else {
            term_order += x - 1;
        }

        if (y == 1) {
            term_order += t1_order;
        } else {
            term_order += y - 1;
        }

        if (term_order > total_order) {
            continue;
        }

        if (x >= y) {
            if (x >= 2) {
                out->emplace_back(Wicab, std::make_tuple(x, 2, 2), std::make_tuple(y, 1, 0), zero_tuple, zero_tuple);
            }
            out->emplace_back(Wicab, std::make_tuple(x, 2, 1), std::make_tuple(y, 1, 1), zero_tuple, zero_tuple);
        }

        if (x > y) {
            out->emplace_back(Wicab, std::make_tuple(x, 1, 1), std::make_tuple(y, 2, 1), zero_tuple, zero_tuple);
            if (y >= 2) {
                out->emplace_back(Wicab, std::make_tuple(x, 1, 0), std::make_tuple(y, 2, 2), zero_tuple, zero_tuple);
            }
        }
    }
}

void assign_Wiabc_t3_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {
    for (int z = skip_t1? 2: 1; z <= std::min(max_excitation, residual - ((skip_t1)? 4: 2)); z++) {

        for (int y = skip_t1? 2: 1; y <= std::min(z, residual - z - ((skip_t1)? 2: 1)); y++) {
            int x = residual - z - y + 1;

            if (x > max_excitation) {
                continue;
            }

            if (x < ((skip_t1)? 2: 1)) {
                continue;
            }

            int term_order = 1;

            if (x == 1) {
                term_order += t1_order;
            } else {
                term_order += x - 1;
            }

            if (y == 1) {
                term_order += t1_order;
            } else {
                term_order += y - 1;
            }

            if (z == 1) {
                term_order += t1_order;
            } else {
                term_order += z - 1;
            }

            if (term_order > total_order) {
                continue;
            }

            if (x >= y && y >= z) {
                out->emplace_back(Wicab, std::make_tuple(x, 1, 1), std::make_tuple(y, 1, 1), std::make_tuple(z, 1, 0), zero_tuple);
            }

            if (x >= y && y > z) {
                out->emplace_back(Wicab, std::make_tuple(x, 1, 1), std::make_tuple(y, 1, 0), std::make_tuple(z, 1, 1), zero_tuple);
            }

            if (x > y && y >= z) {
                out->emplace_back(Wicab, std::make_tuple(x, 1, 0), std::make_tuple(y, 1, 1), std::make_tuple(z, 1, 1), zero_tuple);
            }
        }
    }
}

void assign_Wiabc_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {
    assign_Wiabc_t_terms(residual, max_excitation, out, skip_t1, total_order, t1_order);

    assign_Wiabc_t2_terms(residual, max_excitation, out, skip_t1, total_order, t1_order);

    assign_Wiabc_t3_terms(residual, max_excitation, out, skip_t1, total_order, t1_order);
}

void assign_Wijak_t_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {
    int term_order = residual + 1;

    if (term_order > total_order) {
        return;
    }
    if (max_excitation >= 2 && residual >= 1 && residual + 1 <= max_excitation) {
        out->emplace_back(Wijak, std::make_tuple(residual + 1, 3, 1), zero_tuple, zero_tuple, zero_tuple);
    }
}

void assign_Wijak_t2_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {
    for (int x = std::max(residual - max_excitation, (skip_t1)? 2: 1); x <= max_excitation; x++) {
        int y = residual - x + 1;

        if (y > max_excitation) {
            continue;
        }

        if (y < ((skip_t1)? 2: 1)) {
            continue;
        }

        int term_order = 1;

        if (x == 1) {
            term_order += t1_order;
        } else {
            term_order += x - 1;
        }

        if (y == 1) {
            term_order += t1_order;
        } else {
            term_order += y - 1;
        }

        if (term_order > total_order) {
            continue;
        }

        if (x >= y) {
            out->emplace_back(Wijak, std::make_tuple(x, 2, 1), std::make_tuple(y, 1, 0), zero_tuple, zero_tuple);
            if (x >= 2) {
                out->emplace_back(Wijak, std::make_tuple(x, 2, 0), std::make_tuple(y, 1, 1), zero_tuple, zero_tuple);
            }
        }

        if (x > y) {
            if (y >= 2) {
                out->emplace_back(Wijak, std::make_tuple(x, 1, 1), std::make_tuple(y, 2, 0), zero_tuple, zero_tuple);
            }
            out->emplace_back(Wijak, std::make_tuple(x, 1, 0), std::make_tuple(y, 2, 1), zero_tuple, zero_tuple);
        }
    }
}

void assign_Wijak_t3_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {
    for (int z = skip_t1? 2: 1; z <= std::min(max_excitation, residual - ((skip_t1)? 4: 2)); z++) {

        for (int y = skip_t1? 2: 1; y <= std::min(z, residual - z - ((skip_t1)? 2: 1)); y++) {
            int x = residual - z - y + 1;

            if (x > max_excitation) {
                continue;
            }

            if (x < ((skip_t1)? 2: 1)) {
                continue;
            }

            int term_order = 1;

            if (x == 1) {
                term_order += t1_order;
            } else {
                term_order += x - 1;
            }

            if (y == 1) {
                term_order += t1_order;
            } else {
                term_order += y - 1;
            }

            if (z == 1) {
                term_order += t1_order;
            } else {
                term_order += z - 1;
            }

            if (term_order > total_order) {
                continue;
            }

            if (x >= y && y >= z) {
                out->emplace_back(Wijak, std::make_tuple(x, 1, 1), std::make_tuple(y, 1, 0), std::make_tuple(z, 1, 0), zero_tuple);
            }

            if (x >= y && y > z) {
                out->emplace_back(Wijak, std::make_tuple(x, 1, 0), std::make_tuple(y, 1, 0), std::make_tuple(z, 1, 1), zero_tuple);
            }

            if (x > y && y >= z) {
                out->emplace_back(Wijak, std::make_tuple(x, 1, 0), std::make_tuple(y, 1, 1), std::make_tuple(z, 1, 0), zero_tuple);
            }
        }
    }
}

void assign_Wijak_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {
    assign_Wijak_t_terms(residual, max_excitation, out, skip_t1, total_order, t1_order);

    assign_Wijak_t2_terms(residual, max_excitation, out, skip_t1, total_order, t1_order);

    assign_Wijak_t3_terms(residual, max_excitation, out, skip_t1, total_order, t1_order);
}

void assign_Wabic_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {

    int term_order = residual - 1;

    if (residual == 2) {
        term_order = t1_order + 1;
    }

    if (term_order > total_order) {
        return;
    }
    if (residual - 1 <= max_excitation && !(residual == 2 && skip_t1) && residual >= 2) {
        out->emplace_back(Wacib, std::make_tuple(residual - 1, 1, 1), zero_tuple, zero_tuple, zero_tuple);
    }
}

void assign_Waijk_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {
    int term_order = residual - 1;

    if (residual == 2) {
        term_order = t1_order + 1;
    }

    if (term_order > total_order) {
        return;
    }
    if (residual - 1 <= max_excitation && !(residual == 2 && skip_t1) && residual >= 2) {
        out->emplace_back(Wajik, std::make_tuple(residual - 1, 1, 0), zero_tuple, zero_tuple, zero_tuple);
    }
}

void assign_Wijab_t_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {
    int term_order = residual + 2;

    if (term_order > total_order) {
        return;
    }
    if (residual + 2 <= max_excitation && max_excitation >= 2) {
        out->emplace_back(Wijab, std::make_tuple(residual + 2, 4, 2), zero_tuple, zero_tuple, zero_tuple);
    }
}

void assign_Wijab_t2_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {

    for (int x = std::max(residual - max_excitation, (skip_t1)? 2: 1); x <= max_excitation; x++) {
        int y = residual - x + 2;

        if (y > max_excitation) {
            continue;
        }

        if (y < ((skip_t1)? 2: 1)) {
            continue;
        }

        int term_order = 1;

        if (x == 1) {
            term_order += t1_order;
        } else {
            term_order += x - 1;
        }

        if (y == 1) {
            term_order += t1_order;
        } else {
            term_order += y - 1;
        }

        if (term_order > total_order) {
            continue;
        }

        if (x >= y) {
            if (x >= 2) {
                out->emplace_back(Wijab, std::make_tuple(x, 3, 2), std::make_tuple(y, 1, 0), zero_tuple, zero_tuple);
                out->emplace_back(Wijab, std::make_tuple(x, 3, 1), std::make_tuple(y, 1, 1), zero_tuple, zero_tuple);

                if (y >= 2) {
                    out->emplace_back(Wijab, std::make_tuple(x, 2, 2), std::make_tuple(y, 2, 0), zero_tuple, zero_tuple);
                }
            }
            out->emplace_back(Wijab, std::make_tuple(x, 2, 1), std::make_tuple(y, 2, 1), zero_tuple, zero_tuple);
        }

        if (x > y && y >= 2) {
            out->emplace_back(Wijab, std::make_tuple(x, 1, 0), std::make_tuple(y, 3, 2), zero_tuple, zero_tuple);
            out->emplace_back(Wijab, std::make_tuple(x, 1, 1), std::make_tuple(y, 3, 1), zero_tuple, zero_tuple);
            out->emplace_back(Wijab, std::make_tuple(x, 2, 0), std::make_tuple(y, 2, 2), zero_tuple, zero_tuple);
        }
    }
}

void assign_Wijab_t3_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {

    for (int z = skip_t1? 2: 1; z <= std::min(max_excitation, residual - ((skip_t1)? 4: 2)); z++) {

        for (int y = skip_t1? 2: 1; y <= std::min(z, residual - z - ((skip_t1)? 2: 1)); y++) {
            int x = residual - z - y + 2;

            if (x > max_excitation) {
                continue;
            }

            if (x < ((skip_t1)? 2: 1)) {
                continue;
            }

            int term_order = 1;

            if (x == 1) {
                term_order += t1_order;
            } else {
                term_order += x - 1;
            }

            if (y == 1) {
                term_order += t1_order;
            } else {
                term_order += y - 1;
            }

            if (z == 1) {
                term_order += t1_order;
            } else {
                term_order += z - 1;
            }

            if (term_order > total_order) {
                continue;
            }

            if (x >= y && y >= z) {
                if (x >= 2) {
                    out->emplace_back(Wijab, std::make_tuple(x, 2, 2), std::make_tuple(y, 1, 0), std::make_tuple(z, 1, 0), zero_tuple);
                    out->emplace_back(Wijab, std::make_tuple(x, 2, 0), std::make_tuple(y, 1, 1), std::make_tuple(z, 1, 1), zero_tuple);
                }
                out->emplace_back(Wijab, std::make_tuple(x, 2, 1), std::make_tuple(y, 1, 1), std::make_tuple(z, 1, 0), zero_tuple);
            }

            if (x >= y && y > z) {
                out->emplace_back(Wijab, std::make_tuple(x, 2, 1), std::make_tuple(y, 1, 0), std::make_tuple(z, 1, 1), zero_tuple);

                if (z >= 2) {
                    out->emplace_back(Wijab, std::make_tuple(x, 1, 0), std::make_tuple(y, 1, 0), std::make_tuple(z, 2, 2), zero_tuple);
                    out->emplace_back(Wijab, std::make_tuple(x, 1, 1), std::make_tuple(y, 1, 1), std::make_tuple(z, 2, 0), zero_tuple);
                }
                out->emplace_back(Wijab, std::make_tuple(x, 1, 1), std::make_tuple(y, 1, 0), std::make_tuple(z, 2, 1), zero_tuple);
            }

            if (x > y && y >= z) {
                if (y >= 2) {
                    out->emplace_back(Wijab, std::make_tuple(x, 1, 0), std::make_tuple(y, 2, 2), std::make_tuple(z, 1, 0), zero_tuple);
                    out->emplace_back(Wijab, std::make_tuple(x, 1, 1), std::make_tuple(y, 2, 0), std::make_tuple(z, 1, 1), zero_tuple);
                }
                out->emplace_back(Wijab, std::make_tuple(x, 1, 1), std::make_tuple(y, 2, 1), std::make_tuple(z, 1, 0), zero_tuple);
                out->emplace_back(Wijab, std::make_tuple(x, 1, 0), std::make_tuple(y, 2, 1), std::make_tuple(z, 1, 1), zero_tuple);
            }

            if (x > y && y > z) {
                out->emplace_back(Wijab, std::make_tuple(x, 1, 0), std::make_tuple(y, 1, 1), std::make_tuple(z, 2, 1), zero_tuple);
            }
        }
    }
}

void assign_Wijab_t4_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {

    for (int w = skip_t1? 2: 1; w <= std::min(max_excitation, residual - ((skip_t1)? 6: 3)); w++) {
        for (int z = skip_t1? 2: 1; z <= std::min(max_excitation, residual - w - ((skip_t1)? 4: 2)); z++) {
            for (int y = skip_t1? 2: 1; y <= std::min(max_excitation, residual - w - z - ((skip_t1)? 2: 1)); y++) {
                int x = residual - w - z - y + 2;

                if (x > max_excitation) {
                    continue;
                }

                if (x < ((skip_t1)? 2: 1)) {
                    continue;
                }

                int term_order = 1;

                if (x == 1) {
                    term_order += t1_order;
                } else {
                    term_order += x - 1;
                }

                if (y == 1) {
                    term_order += t1_order;
                } else {
                    term_order += y - 1;
                }

                if (z == 1) {
                    term_order += t1_order;
                } else {
                    term_order += z - 1;
                }

                if (w == 1) {
                    term_order += t1_order;
                } else {
                    term_order += w - 1;
                }

                if (term_order > total_order) {
                    continue;
                }

                if (x >= y && y >= z && z >= w) {
                    out->emplace_back(Wijab, std::make_tuple(x, 1, 1), std::make_tuple(y, 1, 1), std::make_tuple(z, 1, 0),
                            std::make_tuple(w, 1, 0));
                }

                if (x >= y && y > z && z >= w) {
                    out->emplace_back(Wijab, std::make_tuple(x, 1, 1), std::make_tuple(y, 1, 0), std::make_tuple(z, 1, 1),
                            std::make_tuple(w, 1, 0));
                }
                if (x >= y && y >= z && z > w) {
                    out->emplace_back(Wijab, std::make_tuple(x, 1, 1), std::make_tuple(y, 1, 0), std::make_tuple(z, 1, 0),
                            std::make_tuple(w, 1, 1));
                }
                if (x > y && y >= z && z >= w) {
                    out->emplace_back(Wijab, std::make_tuple(x, 1, 0), std::make_tuple(y, 1, 1), std::make_tuple(z, 1, 1),
                            std::make_tuple(w, 1, 0));
                }
                if (x > y && y >= z && z > w) {
                    out->emplace_back(Wijab, std::make_tuple(x, 1, 0), std::make_tuple(y, 1, 1), std::make_tuple(z, 1, 0),
                            std::make_tuple(w, 1, 1));
                }
                if (x >= y && y > z && z >= w) {
                    out->emplace_back(Wijab, std::make_tuple(x, 1, 0), std::make_tuple(y, 1, 0), std::make_tuple(z, 1, 1),
                            std::make_tuple(w, 1, 1));
                }
            }
        }
    }
}

void assign_Wijab_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order, int t1_order) {
    assign_Wijab_t_terms(residual, max_excitation, out, skip_t1, total_order, t1_order);

    assign_Wijab_t2_terms(residual, max_excitation, out, skip_t1, total_order, t1_order);

    assign_Wijab_t3_terms(residual, max_excitation, out, skip_t1, total_order, t1_order);

    assign_Wijab_t4_terms(residual, max_excitation, out, skip_t1, total_order, t1_order);
}

void assign_Fai_terms(int residual, std::list<Diagram> *out) {
    if (residual == 1) {
        out->emplace_back(Fai, zero_tuple, zero_tuple, zero_tuple, zero_tuple);
    }
}

void assign_Wabij_terms(int residual, std::list<Diagram> *out) {
    if (residual == 2) {
        out->emplace_back(Wabij, zero_tuple, zero_tuple, zero_tuple, zero_tuple);
    }
}

}
}

