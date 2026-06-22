/*
 * @BEGIN LICENSE
 *
 * general_dlpno by Psi4 Developer, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*
 * Author: Connor Briggs
 */

#ifndef __PSI_GENERAL_DLPNO_HPP_
#define __PSI_GENERAL_DLPNO_HPP_

#include <cstdint>
#include <vector>
#include <string>

#include "psi4/pragma.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/dlpno/dlpno.h"

namespace psi {
namespace general_dlpno {

class PSI_API Theory {
public:

    enum TheoryType {
        FULL, // This is full treatment of the given excitations.
        PARENTHESIS,// This is full treatment of all excitations but the highest one. The highest is treated perturbatively through the Wigner 2n + 1 rule. This is the parenthesis formalism.
        // When the order is 2, this gives MP2. That is, CC(D) is essentially MP2.
        PERTURB// Treat all excitations but the highest at the last perturbative order that only includes the considered excitations. Treat the highest excitations at the first non-vanishing perturbative order.
        // This is the CCn formalism. Singles are always treated exactly. CC2 treats doubles at 1st order. CC3 treats doubles at 3rd order, triples at 2nd order. CC4 treats doubles and triples at 4th order,
        // quadruples at 3rd order.
    };

    Theory(std::string const &theory_str);

    // 0 is hf. 1 is ccs if full, otherwise, it doesn't make sense. 2 is ccsd/mp2/cc2 depending on the type. 3 and up is what would be expected.
    // 3 would be ccsdt/ccsd(t)/cc3, 4 would be ccsdtq/ccsdt(q)/cc4, etc.
    constexpr Theory(unsigned int level, TheoryType type) : level_ {level}, type_ {type} {

    }

    constexpr unsigned int get_level() const {
        return level_;
    }
    constexpr TheoryType get_type() const {
        return type_;
    }

private:
    unsigned int level_;

    TheoryType type_;
};

class PSI_API GeneralDLPNO : public psi::dlpno::DLPNO {
public:
    GeneralDLPNO(psi::SharedWavefunction ref_wfn, psi::Options &options);

    double compute_energy() override;

protected:
    std::vector<double> cutoffs_;

};
}
}

#endif
