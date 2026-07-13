/*
 * cc_energy.cpp
 *
 *  Created on: Jul 13, 2026
 *      Author: connor
 */




#include <general_cc/general-cc.hpp>
#include <general_cc/tei_transform.hpp>
#include <general_cc/orbit-string.hpp>
#include <Einsums/Tensor/Tensor.hpp>
#include <Einsums/Tensor/DiskTensor.hpp>

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.hpp"

#include <list>

using namespace einsums;
using namespace std;

namespace psi::general_cc {

double GeneralCC::compute_iteration_energy() const {
    for(auto diagram : energy_expression_) {

    }
}





}
