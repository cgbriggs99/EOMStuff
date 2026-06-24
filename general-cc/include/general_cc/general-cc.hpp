/*
 * general-cc.hpp
 *
 *  Created on: Jun 23, 2026
 *      Author: connor
 */

#ifndef GENERAL_CC_INCLUDE_GENERAL_CC_HPP_
#define GENERAL_CC_INCLUDE_GENERAL_CC_HPP_

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.hpp"

#include <diagrams/diagram.hpp>
#include <diagrams/theory-input.hpp>

#include <Einsums/Tensor/Tensor.hpp>
#include <Einsums/Tensor/TiledTensor.hpp>
#include <Einsums/Tensor/BlockTensor.hpp>

#include <vector>

// My editor is complaining. This should be a no-op in real contexts.
#ifndef PSI_API
#define PSI_API
#endif

namespace psi {
namespace general_cc {


namespace detail {
struct PSI_API FullIntegrals {
    einsums::TiledTensor<>
};

}

class PSI_API GeneralCC : Wavefunction {
public:
    GeneralCC(diagram::Theory theory, SharedWavefunction ref_wfn, Options &options);

    virtual ~GeneralCC() = default;

protected:
    diagram::Theory theory;

    einsums::BlockTensor<double, 2> Ca_, Cb_, F_, H_;

    enum IntegralTypes {
        FULL, DENSITY_FITTED
    };

    // Get the Fock terms out of the way.
    einsums::BlockTensor<double, 2> W_b_a_x_, W_j_x_i_, W_ai_x_x_, W_x_a_i_;

};

}
}

#endif /* GENERAL_CC_INCLUDE_GENERAL_CC_HPP_ */
