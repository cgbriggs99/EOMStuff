/*
 * cc_builder.cpp
 *
 *  Created on: Jun 23, 2026
 *      Author: connor
 */

#include <diagrams/diagram.hpp>

#include <list>
#include <stdexcept>

namespace diagram {
std::list<Diagram> compute_cc_residual(int residual, int max_excitation, bool skip_t1, bool canonical) {

    std::list<Diagram> out;

    if (residual == 0 && max_excitation == 0) {
        return out;
    }

    if (residual < 0 || max_excitation < 1 || residual > max_excitation) {
        throw std::runtime_error("Invalid choice for the residual and "
                "the excitation level!");
    }

    detail::assign_Fab_terms(residual, max_excitation, &out, skip_t1);

    detail::assign_Fij_terms(residual, max_excitation, &out, skip_t1);

    if (!canonical) {
        detail::assign_Fia_terms(residual, max_excitation, &out, skip_t1);
    }

    detail::assign_Wabcd_terms(residual, max_excitation, &out, skip_t1);

    detail::assign_Wijkl_terms(residual, max_excitation, &out, skip_t1);

    detail::assign_Waijb_terms(residual, max_excitation, &out, skip_t1);

    detail::assign_Wiabc_terms(residual, max_excitation, &out, skip_t1);

    detail::assign_Wijak_terms(residual, max_excitation, &out, skip_t1);

    detail::assign_Wabic_terms(residual, max_excitation, &out, skip_t1);

    detail::assign_Waijk_terms(residual, max_excitation, &out, skip_t1);

    detail::assign_Wijab_terms(residual, max_excitation, &out, skip_t1);

    if (!canonical) {
        detail::assign_Fai_terms(residual, &out);
    }

    detail::assign_Wabij_terms(residual, &out);

    return out;
}

std::list<Diagram> compute_ccn_residual(int residual, int max_excitation, bool skip_t1, bool canonical) {

    if (residual < max_excitation) {
        return compute_cc_residual(residual, max_excitation, skip_t1, canonical);
    } else {
        std::list<Diagram> out;

        // When the residual and max excitation match, we need to use the
        // perturbative expansion.
        if (residual == 1) {
            detail::assign_Fab_terms(residual, max_excitation, &out, false, max_excitation - 1, 0);
            detail::assign_Fij_terms(residual, max_excitation, &out, false, max_excitation - 1, 0);
        } else {
            detail::assign_Fab_terms(residual, max_excitation, &out, true, max_excitation - 1, 0);
            detail::assign_Fij_terms(residual, max_excitation, &out, true, max_excitation - 1, 0);
        }

        if (!canonical) {
            detail::assign_Fia_terms(residual, max_excitation, &out, true, max_excitation - 1, 0);
        }

        detail::assign_Wabcd_terms(residual, max_excitation - 1, &out, skip_t1, max_excitation - 1, 0);

        detail::assign_Wijkl_terms(residual, max_excitation - 1, &out, skip_t1, max_excitation - 1, 0);

        detail::assign_Waijb_terms(residual, max_excitation - 1, &out, skip_t1, max_excitation - 1, 0);

        detail::assign_Wiabc_terms(residual, max_excitation - 2, &out, skip_t1, max_excitation - 1, 0);

        detail::assign_Wijak_terms(residual, max_excitation - 2, &out, skip_t1, max_excitation - 1, 0);

        detail::assign_Wabic_terms(residual, max_excitation - 1, &out, skip_t1, max_excitation - 1, 0);

        detail::assign_Waijk_terms(residual, max_excitation - 1, &out, skip_t1, max_excitation - 1, 0);

        detail::assign_Wijab_terms(residual, max_excitation - 1, &out, skip_t1, max_excitation - 1, 0);

        if (!canonical) {
            detail::assign_Fai_terms(residual, &out);
        }

        detail::assign_Wabij_terms(residual, &out);

        return out;
    }
}

std::list<Diagram> compute_ccp_residual(int residual, int max_excitation, bool skip_t1, bool canonical) {
    if (residual < max_excitation) {
        if (max_excitation == 2) {
            return compute_cc_residual(residual, max_excitation, skip_t1, canonical);
        }
        return compute_cc_residual(residual, max_excitation - 1, skip_t1, canonical);
    } else {
        std::list<Diagram> out;

        detail::assign_Fab_terms(residual, max_excitation, &out, false, max_excitation - 1, 2);
        detail::assign_Fij_terms(residual, max_excitation, &out, false, max_excitation - 1, 2);

        detail::assign_Wabcd_terms(residual, max_excitation - 1, &out, skip_t1, max_excitation - 1, 2);

        detail::assign_Wijkl_terms(residual, max_excitation - 1, &out, skip_t1, max_excitation - 1, 2);

        detail::assign_Waijb_terms(residual, max_excitation - 1, &out, skip_t1, max_excitation - 1, 2);

        detail::assign_Wiabc_terms(residual, max_excitation - 2, &out, skip_t1, max_excitation - 1, 2);

        detail::assign_Wijak_terms(residual, max_excitation - 2, &out, skip_t1, max_excitation - 1, 2);

        detail::assign_Wabic_terms(residual, max_excitation - 1, &out, skip_t1, max_excitation - 1, 2);

        detail::assign_Waijk_terms(residual, max_excitation - 1, &out, skip_t1, max_excitation - 1, 2);

        detail::assign_Wijab_terms(residual, max_excitation - 1, &out, skip_t1, max_excitation - 1, 2);

        if (!canonical) {
            detail::assign_Fai_terms(residual, &out);
        }

        detail::assign_Wabij_terms(residual, &out);

        return out;
    }
}

std::list<Diagram> compute_lcc_residual(int residual, int max_excitation, bool skip_t1, bool canonical) {
    std::list<Diagram> out;
    if (residual == 0 && max_excitation == 0) {
        return out;
    }
    if (residual < 0 || max_excitation < 1 || residual > max_excitation) {
        throw std::runtime_error("Invalid choice for the residual and "
                "the excitation level!");
    }

    detail::assign_Fab_terms(residual, max_excitation, &out, skip_t1);
    detail::assign_Fij_terms(residual, max_excitation, &out, skip_t1);

    if (!canonical) {
        detail::assign_Fia_t_terms(residual, max_excitation, &out, skip_t1);
    }

    detail::assign_Wabcd_t_terms(residual, max_excitation, &out, skip_t1);
    detail::assign_Wijkl_t_terms(residual, max_excitation, &out, skip_t1);
    detail::assign_Waijb_t_terms(residual, max_excitation, &out, skip_t1);
    detail::assign_Wiabc_t_terms(residual, max_excitation, &out, skip_t1);
    detail::assign_Wijak_t_terms(residual, max_excitation, &out, skip_t1);
    detail::assign_Wabic_terms(residual, max_excitation, &out, skip_t1);
    detail::assign_Waijk_terms(residual, max_excitation, &out, skip_t1);
    detail::assign_Wijab_t_terms(residual, max_excitation, &out, skip_t1);

    if (!canonical) {
        detail::assign_Fai_terms(residual, &out);
    }

    detail::assign_Wabij_terms(residual, &out);

    return out;
}
}

