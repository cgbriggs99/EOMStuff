/*
 * contract.cpp
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

void GeneralCC::resolvent_contract(int num_inds, DiskTensor<double, 2> *residual) const {
    const size_t num_A = OrbitString::number_of_strings(num_inds, nvirt_), num_I = OrbitString::number_of_strings(num_inds, nocc_);

    for (size_t A_block = 0; A_block < num_A / block_size; A_block++) {
        const OrbitString start_A_string { A_block * block_size, num_inds, nvirt_, OrbitString::ALPHA };
        for (size_t I_block = 0; I_block < num_I / block_size; I_block++) {
            auto res_view = (*residual)(Range { A_block * block_size, (A_block + 1) * block_size },
                    Range { I_block * block_size, (I_block + 1) * block_size });

            auto &res_tensor = res_view.get();

            const OrbitString start_I_string { I_block * block_size, num_inds, nocc_, OrbitString::ALPHA };
            OrbitString A_string = start_A_string;
            for (size_t A = 0; A < block_size; A++) {
                OrbitString I_string = start_A_string;

                double virt_denom = 0;
                for (int p = 0; p < num_inds; p++) {
                    virt_denom -= mo_energies_v_->subscript(A_string[p]);
                }

                for (size_t I = 0; I < block_size; I++) {
                    double denom = virt_denom;

                    for (int p = 0; p < num_inds; p++) {
                        denom += mo_energies_o_->subscript(I_string[p]);
                    }

                    res_tensor(A, I) /= denom;
                }
            }

            res_tensor.put();
        }
    }
}

void GeneralCC::contract_interaction(diagram::triplet const &amps_spec, diagram::Interaction interaction,
        einsums::DiskTensor<double, 2> &scaled_t, einsums::DiskTensor<double, 2> *out) const {

    switch(interaction) {
    case diagram::Fab :
        for(size_t )
    }
}

}

