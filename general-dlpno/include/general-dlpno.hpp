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
