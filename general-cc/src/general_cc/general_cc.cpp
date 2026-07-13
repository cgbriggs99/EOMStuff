/*
 * general_cc.cpp
 *
 *  Created on: Jul 9, 2026
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

using namespace einsums;
using namespace std;

namespace psi::general_cc {

GeneralCC::GeneralCC(diagram::Theory theory, SharedWavefunction ref_wfn, Options &options) : Wavefunction(ref_wfn) {
    this->theory_ = theory;

    Ca_ = std::make_shared<Tensor<double, 2>>("Alpha orbitals", ref_wfn->nmo(), ref_wfn->nmo());
    Fa_ = std::make_shared<Tensor<double, 2>>("Alpha Fock matrix", ref_wfn->nmo(), ref_wfn->nmo());
    mo_energies_a_ = std::make_shared<Tensor<double, 1>>("Alpha MO energies", ref_wfn->nmo());

    TensorView<double, 2> psi_Ca_view { ref_wfn->Ca()->pointer(), Dim { nmo(), nmo() } };

    *Ca_ = psi_Ca_view;

    TensorView<double, 2> psi_Fa_view { ref_wfn->Fa()->pointer(), Dim { nmo(), nmo() } };

    *Fa_ = psi_Fa_view;

    TensorView<double, 1> psi_moa_view { ref_wfn->epsilon_a()->pointer(), Dim { nmo() } };

    *mo_energies_a_ = psi_moa_view;

    if (options.get_str_i("REFERENCE") == "RHF") {
        Cb_ = Ca_;
        Fb_ = Fa_;
        mo_energies_b_ = mo_energies_a_;
    } else {
        Cb_ = std::make_shared<Tensor<double, 2>>("Beta orbitals", ref_wfn->nmo(), ref_wfn->nmo());
        Fb_ = std::make_shared<Tensor<double, 2>>("Beta Fock matrix", ref_wfn->nmo(), ref_wfn->nmo());
        mo_energies_b_ = std::make_shared<Tensor<double, 2>>("Beta MO energies", ref_wfn->nmo());

        TensorView<double, 2> psi_Cb_view { ref_wfn->Cb()->pointer(), Dim { nmo(), nmo() } };

        *Cb_ = psi_Cb_view;

        TensorView<double, 2> psi_Fb_view { ref_wfn->Fb()->pointer(), Dim { nmo(), nmo() } };

        *Fb_ = psi_Fb_view;

        TensorView<double, 1> psi_mob_view { ref_wfn->epsilon_b()->pointer(), Dim { nmo() } };

        *mo_energies_b_ = psi_mob_view;
    }

    if (options.get_str_i("CC_TYPE") == "DF") {
        integral_type_ = DENSITY_FITTED;
    } else {
        integral_type_ = FULL;
    }

    nocc_ = nalpha() + nbeta();
    nvirt_ = 2 * nmo() - nalpha() - nbeta();
    navirt_ = nmo() - nalpha();
    nbvirt_ = nmo() - nbeta();

    mo_energies_o_ = std::make_shared<Tensor<double, 1>>("Occupied MO energies", nocc_);

    for (int i = 0; i < nalpha(); i++) {
        mo_energies_o_->subscript(i) = mo_energies_a_->subscript(i);
    }

    for (int i = 0; i < nbeta(); i++) {
        mo_energies_o_->subscript(i + nalpha()) = mo_energies_b_->subscript(i);
    }

    mo_energies_v_ = std::make_shared<Tensor<double, 1>>("Virtual MO energies", nvirt_);

    for (int i = 0; i < navirt_; i++) {
        mo_energies_v_->subscript(i) = mo_energies_a_->subscript(i + nalpha());
    }

    for (int i = 0; i < nbvirt_; i++) {
        mo_energies_v_->subscript(i + navirt_) = mo_energies_b_->subscript(i + nbeta());
    }

    if (theory.get_type() == diagram::Theory::LINEAR) {
        energy_expression_ = diagram::compute_lcc_residual(0, theory.get_level(), false, true);
    } else {
        energy_expression_ = diagram::compute_cc_residual(0, theory.get_level(), false, true);
    }

    t1_transformed_ = theory.get_type() == diagram::Theory::PERTURB || options.get_bool("T1_TRANSFORM");

    residual_expressions_ = std::vector < diagram::FactoredDiagram > (theory.get_level());

    for (int i = 0; i < theory.get_level(); i++) {
        if (i == 0 && theory.do_skip_singles() && !t1_transformed_) {
            continue;
        }
        switch (theory.get_type()) {
        case diagram::Theory::FULL:
            residual_expressions_[i] = diagram::FactoredDiagram(
                    diagram::compute_cc_residual(i + 1, theory.get_level(), theory.do_skip_singles() || t1_transformed_, t1_transformed_));
            break;
        case diagram::Theory::PERTURB:
            residual_expressions_[i] = diagram::FactoredDiagram(diagram::compute_ccn_residual(i + 1, theory.get_level(), true, false));
            break;
        case diagram::Theory::BRACKET:
        case diagram::Theory::PARENTHESIS:
            residual_expressions_[i] = diagram::FactoredDiagram(
                    diagram::compute_ccp_residual(i + 1, theory.get_level(), theory.do_skip_singles() || t1_transformed_, t1_transformed_));
            break;
        case diagram::Theory::LINEAR:
            residual_expressions_[i] = diagram::FactoredDiagram(
                    diagram::compute_lcc_residual(i + 1, theory.get_level(), theory.do_skip_singles() || t1_transformed_, t1_transformed_));
            break;
        }
    }

}

GeneralCC::~GeneralCC() {
    if (integrals_.full != nullptr) {
        if (integral_type_ == FULL) {
            delete integrals_.full;
        } else {
            delete integrals_.df;
        }
    }

    if (t1_integrals_.full != nullptr) {
        if (integral_type_ == FULL) {
            delete t1_integrals_.full;
        } else {
            delete t1_integrals_.df;
        }
    }
}

std::pair<einsums::Tensor<signed char, 2>, einsums::Tensor<ptrdiff_t, 2>> GeneralCC::get_product_table(int nleft, int nright, int norbs) {
    std::string base_name = fmt::format("table_{}_{}_{}", nleft, nright, norbs);

    if (std::filesystem::exists(basename + ".indices.tensors") && std::filesystem::exists(basename + ".mults.tensor")) {
        return read_product_table(base_name);
    } else {
        auto[mults, inds] = generate_product_table(norbs, nleft, nright);

        write_product_table(mults, inds, base_name);

        return {mults, inds};
    }
}

}

