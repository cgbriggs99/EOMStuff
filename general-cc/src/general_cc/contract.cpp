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
                OrbitString I_string = start_I_string;

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
                    if (I_string.has_next()) {
                        ++I_string;
                    }
                }
                if (A_string.has_next()) {
                    ++A_string;
                }
            }

            res_tensor.put();
        }

        // Handle remainder of I
        auto res_view = (*residual)(Range { A_block * block_size, (A_block + 1) * block_size },
                Range { (num_I / block_size) * block_size, (num_I / block_size) * block_size + num_I % block_size });

        auto &res_tensor = res_view.get();

        const OrbitString start_I_string { (num_I / block_size) * block_size, num_inds, nocc_, OrbitString::ALPHA };
        OrbitString A_string = start_A_string;
        for (size_t A = 0; A < block_size; A++) {
            OrbitString I_string = start_A_string;

            double virt_denom = 0;
            for (int p = 0; p < num_inds; p++) {
                virt_denom -= mo_energies_v_->subscript(A_string[p]);
            }

            for (size_t I = 0; I < num_I % block_size; I++) {
                double denom = virt_denom;

                for (int p = 0; p < num_inds; p++) {
                    denom += mo_energies_o_->subscript(I_string[p]);
                }

                res_tensor(A, I) /= denom;

                if (I_string.has_next()) {
                    ++I_string;
                }
            }

            if (A_string.has_next()) {
                ++A_string;
            }
        }

        res_tensor.put();
    }

    // Handle remainder of A
    const OrbitString start_A_string { (num_A / block_size) * block_size, num_inds, nvirt_, OrbitString::ALPHA };
    for (size_t I_block = 0; I_block < num_I / block_size; I_block++) {
        auto res_view = (*residual)(Range { (num_A / block_size) * block_size, (num_A / block_size) * block_size + num_I % block_size },
                Range { I_block * block_size, (I_block + 1) * block_size });

        auto &res_tensor = res_view.get();

        const OrbitString start_I_string { I_block * block_size, num_inds, nocc_, OrbitString::ALPHA };
        OrbitString A_string = start_A_string;
        for (size_t A = 0; A < num_A % block_size; A++) {
            OrbitString I_string = start_I_string;

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

                if (I_string.has_next()) {
                    ++I_string;
                }
            }

            if (A_string.has_next()) {
                ++A_string;
            }
        }

        res_tensor.put();
    }

    auto res_view = (*residual)(Range { (num_A / block_size) * block_size, (num_A / block_size) * block_size + num_A % block_size },
            Range { (num_I / block_size) * block_size, (num_I / block_size) * block_size + num_I % block_size });

    auto &res_tensor = res_view.get();

    const OrbitString start_I_string { (num_I / block_size) * block_size, num_inds, nocc_, OrbitString::ALPHA };
    OrbitString A_string = start_A_string;
    for (size_t A = 0; A < block_size; A++) {
        OrbitString I_string = start_I_string;

        double virt_denom = 0;
        for (int p = 0; p < num_inds; p++) {
            virt_denom -= mo_energies_v_->subscript(A_string[p]);
        }

        for (size_t I = 0; I < num_I % block_size; I++) {
            double denom = virt_denom;

            for (int p = 0; p < num_inds; p++) {
                denom += mo_energies_o_->subscript(I_string[p]);
            }

            res_tensor(A, I) /= denom;

            if (I_string.has_next()) {
                ++I_string;
            }
        }

        if (A_string.has_next()) {
            ++A_string;
        }
    }

    res_tensor.put();
}

void GeneralCC::contract_interaction(diagram::triplet const &amps_spec, diagram::Interaction interaction,
        einsums::DiskTensor<double, 4> *out) const {

    switch (interaction) {
    case diagram::Fab:
        contract_Fab(amps_spec, out);
        break;
    case diagram::Fji:
        contract_Fij(amps_spec, out);
        break;
    case diagram::Fia:
        contract_Fia(amps_spec, out);
        break;
    case diagram::Wadbc:
        contract_Wabcd(amps_spec, out);
        break;
    case diagram::Wjkil:
        contract_Wjkil(amps_spec, out);
        break;
    case diagram::Wjaib:
        contract_Wiajb(amps_spec, out);
        break;
    case diagram::Wicab:
        contract_Wiabc(amps_spec, out);
        break;
    case diagram::Wijak:
        contract_Wijak(amps_spec, out);
        break;
    case diagram::Wacib:
        contract_Wabic(amps_spec, out);
        break;
    case diagram::Wajik:
        contract_Waijk(amps_spec, out);
        break;
    case diagram::Wijab:
        contract_Wijab(amps_spec, out);
        break;
    case diagram::Fai:
        contract_Fai(out);
        break;
    case diagram::Wabij:
        contract_Wabij(out);
        break;
    }

}

void GeneralCC::contract_with_intermediate(diagram::triplet const &amps_spec, einsums::DiskTensor<double, 4> const &intermediate,
        einsums::DiskTensor<double, 4> *out, int Ao_inds, int Io_inds, int C_inds, int K_inds) const {
    int At_inds, It_inds, B_inds, J_inds;
    size_t num_At, num_Ao, num_It, num_Io, num_B, num_J, num_C, num_K;

    // Special case when the amplitude is the T1.
    if (std::get < 0 > (amps_spec) == 1) {

    }
}

}

