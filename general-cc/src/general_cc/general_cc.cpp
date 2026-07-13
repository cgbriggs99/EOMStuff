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

#define MAKE_INDEX(x)                                                                                                                      \
    namespace einsums::index {                                                                                                             \
    struct x : public LabelBase {                                                                                                          \
        static constexpr const char *letter = #x;                                                                                          \
        constexpr x()                       = default;                                                                                     \
        size_t operator()(std::va_list args) const {                                                                                       \
            return va_arg(args, size_t);                                                                                                   \
        }                                                                                                                                  \
                                                                                                                                           \
        size_t operator()(size_t index) const {                                                                                            \
            return index;                                                                                                                  \
        }                                                                                                                                  \
                                                                                                                                           \
        template <typename T, typename Alloc>                                                                                              \
        size_t operator()(std::vector<T, Alloc> *args) const {                                                                             \
            size_t out = args->at(0);                                                                                                      \
            args->erase(args->begin());                                                                                                    \
                                                                                                                                           \
            return out;                                                                                                                    \
        }                                                                                                                                  \
    };                                                                                                                                     \
    static constexpr struct x x;                                                                                                           \
                                                                                                                                           \
    inline auto operator<<(std::ostream &os, const struct x &) -> std::ostream & {                                                         \
        os << x::letter;                                                                                                                   \
        return os;                                                                                                                         \
    }                                                                                                                                      \
    }                                                                                                                                      \
    template <>                                                                                                                            \
    struct fmt::formatter<struct ::einsums::index::x> : fmt::formatter<const char *> {                                                     \
        template <typename FormatContext>                                                                                                  \
        auto format(const struct ::einsums::index::x &, FormatContext &ctx) const {                                                        \
            return formatter<const char *>::format(::einsums::index::x::letter, ctx);                                                      \
        }                                                                                                                                  \
    };

MAKE_INDEX(ij)
MAKE_INDEX(kl)
MAKE_INDEX(ab)
MAKE_INDEX(cd)

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

void GeneralCC::fock_ao_to_so() {
    Tensor<double, 2> F_alpha { "alpha fock", nmo(), nmo() }, F_beta { beta_fock, nmo(), nmo() };

    F_alpha = *Fa_;
    F_beta = *Fb_;

    F_ij_ = std::make_shared<Tensor<double, 2>>("F oo", nocc_, nocc_);
    F_ab_ = std::make_shared<Tensor<double, 2>>("F vv", nvirt_, nvirt_);
    F_ia_ = std::make_shared<Tensor<double, 2>>("F ov", nocc_, nvirt_);
    F_ai_ = std::make_shared<Tensor<double, 2>>("F vo", nvirt_, nocc_);

    t1_F_ij_ = std::make_shared<Tensor<double, 2>>("F oo", nocc_, nocc_);
    t1_F_ab_ = std::make_shared<Tensor<double, 2>>("F vv", nvirt_, nvirt_);
    t1_F_ia_ = std::make_shared<Tensor<double, 2>>("F ov", nocc_, nvirt_);
    t1_F_ai_ = std::make_shared<Tensor<double, 2>>("F vo", nvirt_, nocc_);

    for (int i = 0; i < nocc_; i++) {
        for (int j = 0; j < nocc_; j++) {
            if (i < nalpha()) {
                if (j < nalpha()) {
                    (*F_ij_)(i, j) = F_alpha(i, j);
                } else {
                    (*F_ij_)(i, j) = 0.0;
                }
            } else {
                if (j < nalpha()) {
                    (*F_ij_)(i, j) = 0.0;
                } else {
                    (*F_ij_)(i, j) = F_beta(i - nalpha(), j - nalpha());
                }
            }
        }
    }

    for (int i = 0; i < nocc_; i++) {
        for (int a = 0; a < nvirt_; a++) {
            if (i < nalpha()) {
                if (a < navirt_) {
                    (*F_ia_)(i, a) = F_alpha(i, a + nalpha());
                } else {
                    (*F_ia_)(i, a) = 0.0;
                }
            } else {
                if (a < navirt_) {
                    (*F_ia_)(i, a) = 0.0;
                } else {
                    (*F_ia_)(i, a) = F_beta(i - nalpha(), a - navirt_ + nbeta());
                }
            }
        }
    }

    for (int a = 0; a < nvirt_; a++) {
        for (int i = 0; i < nocc_; i++) {
            if (i < nalpha()) {
                if (a < navirt_) {
                    (*F_ai_)(a, i) = F_alpha(a + nalpha(), i);
                } else {
                    (*F_ai_)(a, i) = 0.0;
                }
            } else {
                if (a < navirt_) {
                    (*F_ai_)(a, i) = 0.0;
                } else {
                    (*F_ai_)(a, i) = F_beta(a - navirt_ + nbeta(), i - nalpha());
                }
            }
        }
    }

    for (int a = 0; a < nvirt_; a++) {
        for (int b = 0; b < nvirt_; b++) {
            if (b < navirt_) {
                if (a < navirt_) {
                    (*F_ab_)(a, b) = F_alpha(a + nalpha(), b + nalpha());
                } else {
                    (*F_ab_)(a, b) = 0.0;
                }
            } else {
                if (a < nvirt_) {
                    (*F_ab_)(a, b) = 0.0;
                } else {
                    (*F_ab_)(a, b) = F_beta(a - navirt_ + nbeta(), b - navirt_ + nbeta());
                }
            }
        }
    }
}

void GeneralCC::tei_ao_to_antisym_so(std::shared_ptr<MintsHelper> mintshelper) {
    if (integral_type_ == FULL) {

        auto psi_matrix = mintshelper->ao_eri();

        const size_t nbf = mintshelper->nbf();

        auto out = std::make_shared<Tensor<double, 4>>("Two electron integrals", nbf, nbf, nbf, nbf);

        auto &unwrapped = *out;
        auto &unwrapped_psi = *psi_matrix;

        auto Ca_occ = (*Ca_)(All, Range { 0, nalpha() }), Cb_occ = (*Cb_)(All, Range { 0, nbeta() });
        auto Ca_virt = (*Ca_)(All, Range { nalpha(), nmo() }), Cb_virt = (*Cb_)(All, Range { nbeta(), nmo() });

        //EINSUMS_OMP_PRAGMA(parallel for collapse(4))
        for (int i = 0; i < nbf; i++) {
            for (int j = 0; j < nbf; j++) {
                for (int k = 0; k < nbf; k++) {
                    for (int l = 0; l < nbf; l++) {
                        unwrapped.subscript(i, j, k, l) = unwrapped_psi(i * nbf + j, k * nbf + l);
                    }
                }
            }
        }

        integrals_.full = new detail::FullIntegrals();

        ConvTEIMethod tei_transform { out };

        {
            auto TEI_aa_abij = tei_transform.contract_tei(Ca_virt, Ca_occ, Ca_virt, Ca_occ);
            auto TEI_ab_abij = tei_transform.contract_tei(Ca_virt, Ca_occ, Cb_virt, Cb_occ);
            auto TEI_bb_abij = tei_transform.contract_tei(Cb_virt, Cb_occ, Cb_virt, Cb_occ);

            integrals_.full->W_ab_ij_x_x = std::make_shared<Tensor<double, 2>>("Sorted <ab||ij>", OrbitString::number_of_strings(2, nvirt_),
                    OrbitString::number_of_strings(2, nocc_));

            integrals_.full->W_ab_ij_x_x->zero();

            for (int a = 0, ab_ind = 0; a < nvirt_; a++) {
                bool spin_a = a < navirt_;
                for (int b = 0; b < a; b++, ab_ind++) {
                    bool spin_b = b < navirt_;
                    for (int i = 0, ij_ind = 0; i < nocc_; i++) {
                        bool spin_i = i < nalpha();
                        for (int j = 0; j < i; j++, ij_ind++) {
                            bool spin_j = j < nalpha();

                            if (spin_a == spin_i && spin_b == spin_j) {
                                if (spin_a && spin_b) {
                                    integrals_.full->W_ab_ij_x_x->subscript(ab_ind, ij_ind) += TEI_aa_abij(a, b, i, j);
                                } else if (spin_a) {
                                    integrals_.full->W_ab_ij_x_x->subscript(ab_ind, ij_ind) += TEI_ab_abij(a, b, i, j);
                                } else if (spin_b) {
                                    integrals_.full->W_ab_ij_x_x->subscript(ab_ind, ij_ind) += TEI_ab_abij(b, a, j, i);
                                } else {
                                    integrals_.full->W_ab_ij_x_x->subscript(ab_ind, ij_ind) += TEI_bb_abij(a, b, i, j);
                                }
                            }

                            if (spin_a == spin_j && spin_b == spin_i) {
                                if (spin_a && spin_b) {
                                    integrals_.full->W_ab_ij_x_x->subscript(ab_ind, ij_ind) -= TEI_aa_abij(a, b, j, i);
                                } else if (spin_a) {
                                    integrals_.full->W_ab_ij_x_x->subscript(ab_ind, ij_ind) -= TEI_ab_abij(a, b, j, i);
                                } else if (spin_b) {
                                    integrals_.full->W_ab_ij_x_x->subscript(ab_ind, ij_ind) -= TEI_ab_abij(b, a, i, j);
                                } else {
                                    integrals_.full->W_ab_ij_x_x->subscript(ab_ind, ij_ind) -= TEI_bb_abij(a, b, j, i);
                                }
                            }
                        }
                    }
                }
            }
        }

        {
            auto TEI_aa_abcd = tei_transform.contract_tei(Ca_virt, Ca_virt, Ca_virt, Ca_virt);
            auto TEI_ab_abcd = tei_transform.contract_tei(Ca_virt, Ca_virt, Cb_virt, Cb_virt);
            auto TEI_bb_abcd = tei_transform.contract_tei(Cb_virt, Cb_virt, Cb_virt, Cb_virt);

            integrals_.full->W_ab_x_cd_x = std::make_shared<Tensor<double, 2>>("Sorted <ab||cd>", OrbitString::number_of_strings(2, nvirt_),
                    OrbitString::number_of_strings(2, nvirt_));

            integrals_.full->W_ab_x_cd_x->zero();

            for (int a = 0, ab_ind = 0; a < nvirt_; a++) {
                bool spin_a = a < navirt_;
                for (int b = 0; b < a; b++, ab_ind++) {
                    bool spin_b = b < navirt_;
                    for (int c = 0, cd_ind = 0; c < nvirt_; c++) {
                        bool spin_c = c < navirt_;
                        for (int d = 0; d < c; d++, cd_ind++) {
                            bool spin_d = d < navirt_;

                            if (spin_a == spin_c && spin_b == spin_c) {
                                if (spin_a && spin_b) {
                                    integrals_.full->W_ab_x_cd_x->subscript(ab_ind, cd_ind) += TEI_aa_abcd(a, b, c, d);
                                } else if (spin_a) {
                                    integrals_.full->W_ab_x_cd_x->subscript(ab_ind, cd_ind) += TEI_ab_abcd(a, b, c, d);
                                } else if (spin_b) {
                                    integrals_.full->W_ab_x_cd_x->subscript(ab_ind, cd_ind) += TEI_ab_abcd(b, a, c, d);
                                } else {
                                    integrals_.full->W_ab_x_cd_x->subscript(ab_ind, cd_ind) += TEI_bb_abcd(a, b, c, d);
                                }
                            }

                            if (spin_a == spin_d && spin_b == spin_c) {
                                if (spin_a && spin_b) {
                                    integrals_.full->W_ab_x_cd_x->subscript(ab_ind, cd_ind) -= TEI_aa_abcd(a, b, d, c);
                                } else if (spin_a) {
                                    integrals_.full->W_ab_x_cd_x->subscript(ab_ind, cd_ind) -= TEI_ab_abcd(a, b, d, c);
                                } else if (spin_b) {
                                    integrals_.full->W_ab_x_cd_x->subscript(ab_ind, cd_ind) -= TEI_ab_abcd(b, a, c, d);
                                } else {
                                    integrals_.full->W_ab_x_cd_x->subscript(ab_ind, cd_ind) -= TEI_bb_abcd(a, b, d, c);
                                }
                            }
                        }
                    }
                }
            }
        }

        {
            auto TEI_aa_ijkl = tei_transform.contract_tei(Ca_occ, Ca_occ, Ca_occ, Ca_occ);
            auto TEI_ab_ijkl = tei_transform.contract_tei(Ca_occ, Ca_occ, Cb_occ, Cb_occ);
            auto TEI_bb_ijkl = tei_transform.contract_tei(Cb_occ, Cb_occ, Cb_occ, Cb_occ);

            integrals_.full->W_x_kl_x_ij = std::make_shared<Tensor<double, 2>>("Sorted <ij||kl>", OrbitString::number_of_strings(2, nocc_),
                    OrbitString::number_of_strings(2, nocc_));

            integrals_.full->W_x_kl_x_ij->zero();

            for (int i = 0, ij_ind = 0; i < nocc_; i++) {
                bool spin_i = i < nalpha();
                for (int j = 0; j < i; j++, ij_ind++) {
                    bool spin_j = j < nalpha();
                    for (int k = 0, kl_ind = 0; k < nocc_; k++) {
                        bool spin_k = k < nalpha();
                        for (int l = 0; l < k; l++, kl_ind++) {
                            bool spin_l = l < nalpha();

                            if (spin_i == spin_k && spin_j == spin_l) {
                                if (spin_i && spin_j) {
                                    integrals_.full->W_x_kl_x_ij->subscript(kl_ind, ij_ind) += TEI_aa_ijkl(i, j, k, l);
                                } else if (spin_i) {
                                    integrals_.full->W_x_kl_x_ij->subscript(kl_ind, ij_ind) += TEI_ab_ijkl(i, j, k, l);
                                } else if (spin_j) {
                                    integrals_.full->W_x_kl_x_ij->subscript(kl_ind, ij_ind) += TEI_ab_ijkl(j, i, l, k);
                                } else {
                                    integrals_.full->W_x_kl_x_ij->subscript(kl_ind, ij_ind) += TEI_bb_ijkl(i, j, k, l);
                                }
                            }

                            if (spin_i == spin_l && spin_j == spin_k) {
                                if (spin_i && spin_j) {
                                    integrals_.full->W_x_kl_x_ij->subscript(kl_ind, ij_ind) -= TEI_aa_ijkl(i, j, l, k);
                                } else if (spin_i) {
                                    integrals_.full->W_x_kl_x_ij->subscript(kl_ind, ij_ind) -= TEI_ab_ijkl(i, j, l, k);
                                } else if (spin_j) {
                                    integrals_.full->W_x_kl_x_ij->subscript(kl_ind, ij_ind) -= TEI_ab_ijkl(j, i, k, l);
                                } else {
                                    integrals_.full->W_x_kl_x_ij->subscript(kl_ind, ij_ind) -= TEI_bb_ijkl(i, j, l, k);
                                }
                            }
                        }
                    }
                }
            }
        }

        {
            auto TEI_aa_ijab = tei_transform.contract_tei(Ca_occ, Ca_virt, Ca_occ, Ca_virt);
            auto TEI_ab_ijab = tei_transform.contract_tei(Ca_occ, Ca_virt, Cb_occ, Cb_virt);
            auto TEI_bb_ijab = tei_transform.contract_tei(Cb_occ, Cb_virt, Cb_occ, Cb_virt);

            integrals_.full->W_x_x_ab_ij = std::make_shared<Tensor<double, 2>>("Sorted <ij||ab>", OrbitString::number_of_strings(2, nvirt_),
                    OrbitString::number_of_strings(2, nocc_));

            integrals_.full->W_x_x_ab_ij->zero();

            for (int i = 0, ij_ind = 0; i < nocc_; i++) {
                bool spin_i = i < nalpha();
                for (int j = 0; j < i; j++, ij_ind++) {
                    bool spin_j = j < nalpha();
                    for (int a = 0, ab_ind = 0; a < nvirt_; a++) {
                        bool spin_a = a < navirt_;
                        for (int b = 0; b < a; b++, ab_ind++) {
                            bool spin_b = b < navirt_;

                            if (spin_i == spin_a && spin_j == spin_b) {
                                if (spin_i && spin_j) {
                                    integrals_.full->W_x_x_ab_ij->subscript(ab_ind, ij_ind) += TEI_aa_ijab(i, j, a, b);
                                } else if (spin_i) {
                                    integrals_.full->W_x_x_ab_ij->subscript(ab_ind, ij_ind) += TEI_ab_ijab(i, j, a, b);
                                } else if (spin_j) {
                                    integrals_.full->W_x_x_ab_ij->subscript(ab_ind, ij_ind) += TEI_ab_ijab(j, i, b, a);
                                } else {
                                    integrals_.full->W_x_x_ab_ij->subscript(ab_ind, ij_ind) += TEI_bb_ijab(i, j, a, b);
                                }
                            }

                            if (spin_i == spin_b && spin_j == spin_a) {
                                if (spin_i && spin_j) {
                                    integrals_.full->W_x_x_ab_ij->subscript(ab_ind, ij_ind) -= TEI_aa_ijab(i, j, b, a);
                                } else if (spin_i) {
                                    integrals_.full->W_x_x_ab_ij->subscript(ab_ind, ij_ind) -= TEI_ab_ijab(i, j, b, a);
                                } else if (spin_j) {
                                    integrals_.full->W_x_x_ab_ij->subscript(ab_ind, ij_ind) -= TEI_ab_ijab(j, i, b, a);
                                } else {
                                    integrals_.full->W_x_x_ab_ij->subscript(ab_ind, ij_ind) -= TEI_bb_ijab(i, j, b, a);
                                }
                            }
                        }
                    }
                }
            }
        }

        {
            auto TEI_aa_abci = tei_transform.contract_tei(Ca_virt, Ca_virt, Ca_virt, Ca_occ);
            auto TEI_ab_abci = tei_transform.contract_tei(Ca_virt, Ca_virt, Cb_virt, Cb_occ);
            auto TEI_ba_abci = tei_transform.contract_tei(Cb_virt, Cb_virt, Ca_virt, Ca_occ);
            auto TEI_bb_abci = tei_transform.contract_tei(Cb_virt, Cb_virt, Cb_virt, Cb_occ);

            integrals_.full->W_ab_i_c_x = std::make_shared<Tensor<double, 3>>("Sorted <ab||ci>", OrbitString::number_of_strings(2, nvirt_),
                    nocc_, nvirt_);

            integrals_.full->W_ab_i_c_x->zero();

            for (int a = 0, ab_ind = 0; a < nvirt_; a++) {
                bool spin_a = a < navirt_;
                for (int b = 0; b < a; b++, ab_ind++) {
                    bool spin_b = b < navirt_;
                    for (int c = 0; c < nvirt_; c++) {
                        bool spin_c = c < navirt_;
                        for (int i = 0; i < nocc_; i++) {
                            bool spin_i = i < nalpha();

                            if (spin_a == spin_c && spin_b == spin_i) {
                                if (spin_a && spin_b) {
                                    integrals_.full->W_ab_i_c_x->subscript(ab_ind, i, c) += TEI_aa_abci(a, b, c, i);
                                } else if (spin_a) {
                                    integrals_.full->W_ab_i_c_x->subscript(ab_ind, i, c) += TEI_ab_abci(a, b, c, i);
                                } else if (spin_b) {
                                    integrals_.full->W_ab_i_c_x->subscript(ab_ind, i, c) += TEI_ba_abci(a, b, c, i);
                                } else {
                                    integrals_.full->W_ab_i_c_x->subscript(ab_ind, i, c) += TEI_bb_abci(a, b, c, i);
                                }
                            }

                            if (spin_a == spin_i && spin_b == spin_c) {
                                if (spin_a && spin_b) {
                                    integrals_.full->W_ab_i_c_x->subscript(ab_ind, i, c) -= TEI_aa_abci(b, a, c, i);
                                } else if (spin_a) {
                                    integrals_.full->W_ab_i_c_x->subscript(ab_ind, i, c) -= TEI_ab_abci(b, a, c, i);
                                } else if (spin_b) {
                                    integrals_.full->W_ab_i_c_x->subscript(ab_ind, i, c) -= TEI_ba_abci(b, a, c, i);
                                } else {
                                    integrals_.full->W_ab_i_c_x->subscript(ab_ind, i, c) -= TEI_bb_abci(b, a, c, i);
                                }
                            }
                        }
                    }
                }
            }
        }

        {
            auto TEI_aa_iajk = tei_transform.contract_tei(Ca_occ, Ca_occ, Ca_virt, Ca_occ);
            auto TEI_ab_iajk = tei_transform.contract_tei(Ca_occ, Ca_occ, Cb_virt, Cb_occ);
            auto TEI_ba_iajk = tei_transform.contract_tei(Cb_occ, Cb_occ, Ca_virt, Ca_occ);
            auto TEI_bb_iajk = tei_transform.contract_tei(Cb_occ, Cb_occ, Cb_virt, Cb_occ);

            integrals_.full->W_a_jk_x_i = std::make_shared<Tensor<double, 3>>("Sorted <ia||jk>", nvirt_,
                    OrbitString::number_of_strings(2, nocc_), nocc_);

            integrals_.full->W_a_jk_x_i->zero();

            for (int i = 0; i < nocc_; i++) {
                bool spin_i = i < nalpha();
                for (int a = 0; a < nvirt_; a++) {
                    bool spin_a = a < navirt_;
                    for (int j = 0, jk_ind = 0; j < nocc_; j++) {
                        bool spin_j = j < nalpha();
                        for (int k = 0; k < j; k++, jk_ind++) {
                            bool spin_k = k < nalpha();

                            if (spin_i == spin_j && spin_a == spin_k) {
                                if (spin_i && spin_a) {
                                    integrals_.full->W_a_jk_x_i->subscript(a, jk_ind, i) += TEI_aa_iajk(i, a, j, k);
                                } else if (spin_i) {
                                    integrals_.full->W_a_jk_x_i->subscript(a, jk_ind, i) += TEI_ab_iajk(i, a, j, k);
                                } else if (spin_a) {
                                    integrals_.full->W_a_jk_x_i->subscript(a, jk_ind, i) += TEI_ba_iajk(i, a, j, k);
                                } else {
                                    integrals_.full->W_a_jk_x_i->subscript(a, jk_ind, i) += TEI_bb_iajk(i, a, j, k);
                                }
                            }

                            if (spin_i == spin_k && spin_a == spin_j) {
                                if (spin_i && spin_a) {
                                    integrals_.full->W_a_jk_x_i->subscript(a, jk_ind, i) -= TEI_aa_iajk(i, a, k, j);
                                } else if (spin_i) {
                                    integrals_.full->W_a_jk_x_i->subscript(a, jk_ind, i) -= TEI_ab_iajk(i, a, k, j);
                                } else if (spin_a) {
                                    integrals_.full->W_a_jk_x_i->subscript(a, jk_ind, i) -= TEI_ba_iajk(i, a, k, j);
                                } else {
                                    integrals_.full->W_a_jk_x_i->subscript(a, jk_ind, i) -= TEI_bb_iajk(i, a, k, j);
                                }
                            }
                        }
                    }
                }
            }
        }

        {
            auto TEI_aa_aibc = tei_transform.contract_tei(Ca_virt, Ca_virt, Ca_occ, Ca_virt);
            auto TEI_ab_aibc = tei_transform.contract_tei(Ca_virt, Ca_virt, Cb_occ, Cb_virt);
            auto TEI_ba_aibc = tei_transform.contract_tei(Cb_virt, Cb_virt, Ca_occ, Ca_virt);
            auto TEI_bb_aibc = tei_transform.contract_tei(Cb_virt, Cb_virt, Cb_occ, Cb_virt);

            integrals_.full->W_a_x_bc_i = std::make_shared<Tensor<double, 3>>("Sorted <ai||bc>", nvirt_,
                    OrbitString::number_of_strings(2, nvirt_), nocc_);

            integrals_.full->W_a_x_bc_i->zero();

            for (int a = 0; a < nvirt_; a++) {
                bool spin_a = a < navirt_;
                for (int i = 0; i < nocc_; i++) {
                    bool spin_i = i < nalpha();
                    for (int b = 0, bc_ind = 0; b < nvirt_; b++) {
                        bool spin_b = b < navirt_;
                        for (int c = 0; c < b; c++, bc_ind++) {
                            bool spin_c = c < navirt_;

                            if (spin_a == spin_b && spin_i == spin_c) {
                                if (spin_a && spin_b) {
                                    integrals_.full->W_a_x_bc_i->subscript(a, bc_ind, i) += TEI_aa_aibc(a, i, b, c);
                                } else if (spin_a) {
                                    integrals_.full->W_a_x_bc_i->subscript(a, bc_ind, i) += TEI_ab_aibc(a, i, b, c);
                                } else if (spin_b) {
                                    integrals_.full->W_a_x_bc_i->subscript(a, bc_ind, i) += TEI_ba_aibc(a, i, b, c);
                                } else {
                                    integrals_.full->W_a_x_bc_i->subscript(a, bc_ind, i) += TEI_bb_aibc(a, i, b, c);
                                }
                            }

                            if (spin_a == spin_i && spin_b == spin_c) {
                                if (spin_a && spin_b) {
                                    integrals_.full->W_a_x_bc_i->subscript(a, bc_ind, i) -= TEI_aa_aibc(a, i, c, b);
                                } else if (spin_a) {
                                    integrals_.full->W_a_x_bc_i->subscript(a, bc_ind, i) -= TEI_ab_aibc(a, i, c, b);
                                } else if (spin_b) {
                                    integrals_.full->W_a_x_bc_i->subscript(a, bc_ind, i) -= TEI_ba_aibc(a, i, c, b);
                                } else {
                                    integrals_.full->W_a_x_bc_i->subscript(a, bc_ind, i) -= TEI_bb_aibc(a, i, c, b);
                                }
                            }
                        }
                    }
                }
            }
        }

        {
            auto TEI_aa_ijka = tei_transform.contract_tei(Ca_occ, Ca_occ, Ca_occ, Ca_virt);
            auto TEI_ab_ijka = tei_transform.contract_tei(Ca_occ, Ca_occ, Cb_occ, Cb_virt);
            auto TEI_ba_ijka = tei_transform.contract_tei(Cb_occ, Cb_occ, Ca_occ, Ca_virt);
            auto TEI_bb_ijka = tei_transform.contract_tei(Cb_occ, Cb_occ, Cb_occ, Cb_virt);

            integrals_.full->W_x_k_a_ij = std::make_shared<Tensor<double, 3>>("<ij||ka>", nocc_, nvirt_,
                    OrbitString::number_of_strings(2, nocc_));

            integrals_.full->W_x_k_a_ij->zero();

            for (int i = 0, ij_ind = 0; i < nocc_; i++) {
                bool spin_i = i < nalpha();
                for (int j = 0; j < i; j++, ij_ind++) {
                    bool spin_j = j < nalpha();
                    for (int k = 0; k < nocc_; k++) {
                        bool spin_k = k < nalpha();
                        for (int a = 0; a < nvirt_; a++) {
                            bool spin_a = a < navirt_;

                            if (spin_i == spin_k && spin_j == spin_a) {
                                if (spin_i && spin_j) {
                                    integrals_.full->W_x_k_a_ij->subscript(k, a, ij_ind) += TEI_aa_ijka(i, j, k, a);
                                } else if (spin_i) {
                                    integrals_.full->W_x_k_a_ij->subscript(k, a, ij_ind) += TEI_ab_ijka(i, j, k, a);
                                } else if (spin_j) {
                                    integrals_.full->W_x_k_a_ij->subscript(k, a, ij_ind) += TEI_ba_ijka(i, j, k, a);
                                } else {
                                    integrals_.full->W_x_k_a_ij->subscript(k, a, ij_ind) += TEI_bb_ijka(i, j, k, a);
                                }
                            }

                            if (spin_i == spin_a && spin_j == spin_k) {
                                if (spin_i && spin_j) {
                                    integrals_.full->W_x_k_a_ij->subscript(k, a, ij_ind) -= TEI_aa_ijka(j, i, k, a);
                                } else if (spin_i) {
                                    integrals_.full->W_x_k_a_ij->subscript(k, a, ij_ind) -= TEI_ab_ijka(j, i, k, a);
                                } else if (spin_j) {
                                    integrals_.full->W_x_k_a_ij->subscript(k, a, ij_ind) -= TEI_ba_ijka(j, i, k, a);
                                } else {
                                    integrals_.full->W_x_k_a_ij->subscript(k, a, ij_ind) -= TEI_bb_ijka(j, i, k, a);
                                }
                            }
                        }
                    }
                }
            }
        }

        {
            auto TEI_aa_aijb = tei_transform.contract_tei(Ca_virt, Ca_occ, Ca_occ, Ca_virt);
            auto TEI_ab_aijb = tei_transform.contract_tei(Ca_virt, Ca_occ, Cb_occ, Cb_virt);
            auto TEI_ba_aijb = tei_transform.contract_tei(Cb_virt, Cb_occ, Ca_occ, Ca_virt);
            auto TEI_bb_aijb = tei_transform.contract_tei(Cb_virt, Cb_occ, Cb_occ, Cb_virt);

            integrals_.full->W_a_j_b_i = std::make_shared<Tensor<double, 4>>("Sorted <ai||jb>", nvirt_, nocc_, nvirt_, nocc_);

            integrals_.full->W_a_j_b_i->zero();

            for (int a = 0; a < nvirt_; a++) {
                bool spin_a = a < navirt_;
                for (int j = 0; j < nocc_; j++) {
                    bool spin_j = j < nalpha();
                    for (int b = 0; b < nvirt_; b++) {
                        bool spin_b = b < navirt_;
                        for (int i = 0; i < nocc_; i++) {
                            bool spin_i = i < nalpha();

                            if (spin_a == spin_j && spin_i == spin_b) {
                                if (spin_a && spin_i) {
                                    integrals_.full->W_a_j_b_i->subscript(a, j, b, i) += TEI_aa_aijb(a, i, j, b);
                                } else if (spin_a) {
                                    integrals_.full->W_a_j_b_i->subscript(a, j, b, i) += TEI_ab_aijb(a, i, j, b);
                                } else if (spin_i) {
                                    integrals_.full->W_a_j_b_i->subscript(a, j, b, i) += TEI_ba_aijb(a, i, j, b);
                                } else {
                                    integrals_.full->W_a_j_b_i->subscript(a, j, b, i) += TEI_bb_aijb(a, i, j, b);
                                }
                            }
                        }
                    }
                }
            }
        }

        {
            auto TEI_aa_aibj = tei_transform.contract_tei(Ca_virt, Ca_virt, Ca_occ, Ca_occ);
            auto TEI_ab_aibj = tei_transform.contract_tei(Ca_virt, Ca_virt, Cb_occ, Cb_occ);
            auto TEI_ba_aibj = tei_transform.contract_tei(Cb_virt, Cb_virt, Ca_occ, Ca_occ);
            auto TEI_bb_aibj = tei_transform.contract_tei(Cb_virt, Cb_virt, Cb_occ, Cb_occ);

            for (int a = 0; a < nvirt_; a++) {
                bool spin_a = a < navirt_;
                for (int j = 0; j < nocc_; j++) {
                    bool spin_j = j < nalpha();
                    for (int b = 0; b < nvirt_; b++) {
                        bool spin_b = b < navirt_;
                        for (int i = 0; i < nocc_; i++) {
                            bool spin_i = i < nalpha();

                            if (spin_a == spin_b && spin_i == spin_j) {
                                if (spin_a && spin_i) {
                                    integrals_.full->W_a_j_b_i->subscript(a, j, b, i) -= TEI_aa_aibj(a, i, b, j);
                                } else if (spin_a) {
                                    integrals_.full->W_a_j_b_i->subscript(a, j, b, i) -= TEI_ab_aibj(a, i, b, j);
                                } else if (spin_i) {
                                    integrals_.full->W_a_j_b_i->subscript(a, j, b, i) -= TEI_ba_aibj(a, i, b, j);
                                } else {
                                    integrals_.full->W_a_j_b_i->subscript(a, j, b, i) -= TEI_bb_aibj(a, i, b, j);
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
        auto aux = mintshelper_->get_basisset("RIFIT");
        auto three_index = mintshelper_->ao_eri(aux, psi::BasisSet::zero_ao_basis_set(), mintshelper_->basisset(),
                mintshelper_->basisset());
        auto two_index = mintshelper_->ao_eri(aux, psi::BasisSet::zero_ao_basis_set(), aux, psi::BasisSet::zero_ao_basis_set());

        const size_t nbf = mintshelper_->nbf();

        auto out = std::make_shared<Tensor<double, 3>>("Density Fitting Tensor", aux->nbf(), mintshelper_->nbf(), mintshelper_->nbf());

        // Wrap the psi4 matrix in a tensor view so that we can just directly copy.
        TensorView<double, 3> psi_wrap_3_index(three_index->get_pointer(), Dim { aux->nbf(), nbf, nbf }, Stride { nbf * nbf, nbf, 1 });
        TensorView<double, 2> psi_wrap_2_index(two_index->get_pointer(), Dim { aux->nbf(), aux->nbf() });

        auto Jhalf = linear_algebra::pow(psi_wrap_2_index, -0.5);

        tensor_algebra::einsum(Indices { index::Q, index::mu, index::nu }, out.get(), Indices { index::P, index::mu, index::nu },
                psi_wrap_3_index, Indices { index::Q, index::P }, Jhalf);

        auto Ca_occ = (*Ca_)(All, Range { 0, nalpha() }), Cb_occ = (*Cb_)(All, Range { 0, nbeta() });
        auto Ca_virt = (*Ca_)(All, Range { nalpha(), nmo() }), Cb_virt = (*Cb_)(All, Range { nbeta(), nmo() });

        DFMethod tei_method { out };

        integrals_.df = new detail::DFIntegrals();
        integrals_.df.naux = aux->nbf();

        {
            auto B_a_Qab = tei_method.contract_B(Ca_virt, Ca_virt);
            auto B_b_Qab = tei_method.contract_B(Cb_virt, Cb_virt);

            integrals_.df->B_abQ = std::make_shared<Tensor<double, 3>>("B_abQ", nvirt_, nvirt_, naux);
            integrals_.df->B_abQ.zero();

            for (int a = 0; a < nvirt_; a++) {
                bool spin_a = a < navirt_;
                for (int b = 0; b < nvirt_; b++) {
                    bool spin_b = b < navirt_;
                    for (int Q = 0; Q < nbf; Q++) {
                        if (spin_a == spin_b) {
                            if (spin_a) {
                                integrals_.df->B_abQ->subscript(a, b, Q) = B_a_Qab(Q, a, b);
                            } else {
                                integrals_.df->B_abQ->subscript(a, b, Q) = B_b_Qab(Q, a, b);
                            }
                        }
                    }
                }
            }
        }

        {
            auto B_a_Qai = tei_method.contract_B(Ca_virt, Ca_occ);
            auto B_b_Qai = tei_method.contract_B(Cb_virt, Cb_occ);

            integrals_.df->B_aiQ = std::make_shared<Tensor<double, 3>>("B_aiQ", nvirt_, nocc_, naux);
            integrals_.df->B_aiQ.zero();

            for (int a = 0; a < nvirt_; a++) {
                bool spin_a = a < navirt_;
                for (int i = 0; i < nocc_; i++) {
                    bool spin_i = i < nalpha();
                    for (int Q = 0; Q < nbf; Q++) {
                        if (spin_a == spin_i) {
                            if (spin_a) {
                                integrals_.df->B_aiQ->subscript(a, i, Q) = B_a_Qai(Q, a, i);
                            } else {
                                integrals_.df->B_aiQ->subscript(a, i, Q) = B_b_Qai(Q, a, i);
                            }
                        }
                    }
                }
            }
        }

        {
            auto B_a_Qia = tei_method.contract_B(Ca_occ, Ca_virt);
            auto B_b_Qia = tei_method.contract_B(Cb_occ, Cb_virt);

            integrals_.df->B_iaQ = std::make_shared<Tensor<double, 3>>("B_iaQ", nocc_, nvirt_, naux);
            integrals_.df->B_iaQ.zero();

            for (int i = 0; i < nocc_; i++) {
                bool spin_i = i < nalpha();
                for (int a = 0; a < nvirt_; a++) {
                    bool spin_a = a < navirt_;
                    for (int Q = 0; Q < nbf; Q++) {
                        if (spin_a == spin_i) {
                            if (spin_a) {
                                integrals_.df->B_iaQ->subscript(i, a, Q) = B_a_Qia(Q, i, a);
                            } else {
                                integrals_.df->B_iaQ->subscript(i, a, Q) = B_b_Qia(Q, i, a);
                            }
                        }
                    }
                }
            }
        }

        {
            auto B_a_Qij = tei_method.contract_B(Ca_occ, Ca_occ);
            auto B_b_Qij = tei_method.contract_B(Cb_occ, Cb_occ);

            integrals_.df->B_ijQ = std::make_shared<Tensor<double, 3>>("B_ijQ", nocc_, nvirt_, naux);
            integrals_.df->B_ijQ.zero();

            for (int i = 0; i < nocc_; i++) {
                bool spin_i = i < nalpha();
                for (int j = 0; j < nocc_; j++) {
                    bool spin_j = j < nalpha();
                    for (int Q = 0; Q < nbf; Q++) {
                        if (spin_i == spin_j) {
                            if (spin_i) {
                                integrals_.df->B_ijQ->subscript(i, j, Q) = B_a_Qia(Q, i, j);
                            } else {
                                integrals_.df->B_ijQ->subscript(i, j, Q) = B_b_Qia(Q, i, j);
                            }
                        }
                    }
                }
            }
        }

    }
}

void GeneralCC::t1_transform() {
    *t1_F_ab_ = *F_ab_;
    *t1_F_ij_ = *F_ij_;
    *t1_F_ia_ = *F_ia_;
    *t1_F_ai_ = *F_ai_;

    for (int i = 0; i < nocc_; i++) {
        if (i < nalpha()) {
            t1_F_ij_->subscript(i, i) -= mo_energies_a_(i);
        } else {
            t1_F_ij_->subscript(i, i) -= mo_energies_b_(i - nalpha());
        }
    }

    for (int a = 0; a < nvirt_; a++) {
        if (a < navirt_) {
            t1_F_ab_->subscript(a, a) -= mo_energies_a_(a + nocc_);
        } else {
            t1_F_ab_->subscript(a, a) -= mo_energies_b_(a + nocc_ - navirt_);
        }
    }

    if (integral_type_ == FULL) {
        if (t1_integrals_.full == nullptr) {
            t1_integrals_.full = new detail::FullIntegrals();

            // First one is trivial.
            t1_integrals_.full->W_x_x_ab_ij = integrals_.full->W_x_x_ab_ij;

            t1_integrals_.full->W_ab_ij_x_x = std::make_shared<Tensor<double, 2>>("Sorted t1 <ab||ij>",
                    OrbitString::number_of_strings(2, nvirt_), OrbitString::number_of_strings(2, nocc_));

            t1_integrals_.full->W_ab_x_cd_x = std::make_shared<Tensor<double, 2>>("Sorted t1 <ab||cd>",
                    OrbitString::number_of_strings(2, nvirt_), OrbitString::number_of_strings(2, nvirt_));

            t1_integrals_.full->W_x_kl_x_ij = std::make_shared<Tensor<double, 2>>("Sorted t1 <ij||kl>",
                    OrbitString::number_of_strings(2, nocc_), OrbitString::number_of_strings(2, nocc_));

            t1_integrals_.full->W_ab_i_c_x = std::make_shared<Tensor<double, 3>>("Sorted t1 <ab||ci>",
                    OrbitString::number_of_strings(2, nvirt_), nocc_, nvirt_);

            t1_integrals_.full->W_a_jk_x_i = std::make_shared<Tensor<double, 3>>("Sorted t1 <ia||jk>", nvirt_,
                    OrbitString::number_of_strings(2, nocc_), nocc_);

            t1_integrals_.full->W_a_x_bc_i = std::make_shared<Tensor<double, 3>>("Sorted t1 <ai||bc>", nvirt_,
                    OrbitString::number_of_strings(2, nvirt_), nocc_);

            t1_integrals_.full->W_x_k_a_ij = std::make_shared<Tensor<double, 3>>("Sorted t1 <ij||ka>", nocc_, nvirt_,
                    OrbitString::number_of_strings(2, nocc_));
            t1_integrals_.full->W_a_j_b_i = std::make_shared<Tensor<double, 3>>("Sorted t1 <ai||jb>", nvirt_, nocc_, nvirt_, nocc_);

        }

        *t1_integrals_.full->W_ab_ij_x_x = *integrals_.full->W_ab_ij_x_x;
        *t1_integrals_.full->W_ab_x_cd_x = *integrals_.full->W_ab_x_cd_x;
        *t1_integrals_.full->W_x_kl_x_ij = *integrals_.full->W_x_kl_x_ij;
        *t1_integrals_.full->W_ab_i_c_x = *integrals_.full->W_ab_i_c_x;
        *t1_integrals_.full->W_a_jk_x_i = *integrals_.full->W_a_jk_x_i;
        *t1_integrals_.full->W_a_x_bc_i = *integrals_.full->W_a_x_bc_i;
        *t1_integrals_.full->W_x_k_a_ij = *integrals_.full->W_x_k_a_ij;
        *t1_integrals_.full->W_a_j_b_i = *integrals_.full->W_a_j_b_i;

        auto[oo_mults, oo_inds] = get_product_table(1, 1, nocc_);
        auto[vv_mults, vv_inds] = get_product_table(1, 1, nvirt_);

        // Contract.

        // Fia
        for (int k = 0; k < nocc_; k++) {
            for (int e = 0; e < nvirt_; e++) {
#pragma omp parallel for collapse(2)
                for (int i = 0; i < nocc_; i++) {
                    for (int a = 0; a < nvirt_; a++) {
                        t1_F_ia_->subscript(i, a) += oo_mults(k, i) * vv_mults(e, a) * t1_amps_->subscript(e, k)
                                * integrals_.full->W_x_x_ab_ij->subscript(vv_inds(e, a), oo_inds(i, a));
                    }
                }
            }
        }

        // Fij
        tensor_algebra::einsum(1.0, Indices { index::i, index::j }, t1_F_ij_.get(), 1.0, Indices { index::j, index::e }, *t1_amps_,
                Indices { index::i, index::e }, *t1_F_ia_);

        for (int i = 0; i < nocc_; i++) {
            auto F_view = (*t1_F_ij_)(i, All);
            for (int k = 0; k < nocc_; k++) {
                if (oo_mults(i, k) == 0) {
                    continue;
                }
                auto W_view = (*integrals_.full->W_x_k_a_ij)(All, oo_inds(i, k), All);
                auto t1_view = (*t1_amps_)(All, k);
                tensor_algebra::einsum(1.0, Indices { index::j }, &F_view, (double) oo_mults(i, k), Indices { index::j, index::e }, W_view,
                        Indices { index::e }, t1_view);
            }
        }

        // Fab
        tensor_algebra::einsum(1.0, Indices { index::a, index::b }, t1_F_ab_.get(), -1.0, Indices { index::m, index::b }, *t1_F_ia_,
                Indices { index::m, index::a }, *t1_amps_);

        for (int b = 0; b < nvirt_; b++) {
            auto F_view = (*t1_F_ab_)(All, b);

            for (int e = 0; e < nvirt_; e++) {
                if (vv_mults(b, e) == 0) {
                    continue;
                }

                auto W_view = (*integrals_.full->W_a_x_bc_i)(All, vv_inds(b, e), All);
                auto t1_view = (*t1_amps_)(e, All);
                tensor_algebra::einsum(1.0, Indices { index::a }, &F_view, (double) vv_mults(b, e), Indices { index::a, index::k }, W_view,
                        Indices { index::k }, t1_view);
            }
        }

        // Fai
        tensor_algebra::einsum(1.0, Indices { index::a, index::i }, t1_F_ai_.get(), 1.0, Indices { index::a, index::e }, *t1_F_ab_,
                Indices { index::e, index::i }, *t1_amps_);

        {
            Tensor<double, 2> temp_mi { "temp oo", nocc_, nocc_ };

            temp_mi = *t1_F_ij_;

            tensor_algebra::einsum(1.0, Indices { index::m, index::i }, &temp_mi, -1.0, Indices { index::m, index::e }, *t1_F_ij_, Indices {
                    index::e, index::i }, *t1_amps_);

            tensor_algebra::einsum(1.0, Indices { index::a, index::i }, t1_F_ai_.get(), -1.0, Indices { index::m, index::i }, temp_mi,
                    Indices { index::a, index::m }, *t1_amps_);
        }

        tensor_algebra::einsum(1.0, Indices { index::a, index::i }, t1_F_ai_.get(), 1.0, Indices { index::a, index::i, index::e, index::k },
                *integrals_.full->W_a_j_b_i, Indices { index::e, index::k }, *t1_amps_);

        // <ij||ka>
        for (int a = 0; a < nvirt_; a++) {
            auto out_view = (*t1_integrals_.full->W_x_k_a_ij)(All, a, All);

            for (int e = 0; e < nvirt_; e++) {
                if (vv_mults(e, a) == 0) {
                    continue;
                }

                auto W_view = (*integrals_.full->W_x_x_ab_ij)(vv_inds(e, a), All);
                auto t1_view = (*t1_amps_)(e, All);
                tensor_algebra::einsum(1.0, Indices { index::k, index::ij }, &out_view, (double) vv_mults(e, a), Indices { index::ij },
                        W_view, Indices { index::k }, t1_view);
            }
        }

        // <ai||bc>
        for (int i = 0; i < nocc_; i++) {
            auto out_view = (*t1_integrals_.full->W_a_x_bc_i)(All, All, i);
            for (int m = 0; m < nocc_; m++) {
                if (oo_mults(i, m) == 0) {
                    continue;
                }
                auto W_view = (*integrals_.full->W_x_x_ab_ij)(All, oo_inds(i, m));
                auto t1_view = (*t1_amps_)(All, m);

                tensor_algebra::einsum(1.0, Indices { index::a, index::bc }, &out_view, (double) oo_mults(i, m), Indices { index::bc },
                        W_view, Indices { index::a }, t1_view);

            }
        }

        // <ij||kl>
        for (int k = 0; k < nocc_; k++) {
            auto W_view = (*t1_integrals_.full->W_x_k_a_ij)(k, All, All);
            for (int l = 0; l < nocc_; l++) {
                if (k == l) {
                    continue;
                }
                auto out_view = (*t1_integrals_.full->W_x_kl_x_ij)(oo_inds(k, l), All);
                auto t1_view = (*t1_amps_)(All, l);

                tensor_algebra::einsum(1.0, Indices { index:ij }, &out_view, (double) oo_mults(k, l), Indices { index::e, index::ij },
                        W_view, Indices { index::e }, t1_view);
            }
        }

        // <ab||cd>
        for (int a = 0; a < nvirt_; a++) {
            auto W_view = (*t1_integrals_.full->W_a_x_bc_i)(a, All, All);
            for (int b = 0; b < nvirt_; b++) {
                if (vv_mults(a, b) == 0) {
                    continue;
                }
                auto out_view = (*t1_integrals_.full->W_ab_x_cd_x)(vv_inds(a, b), All);
                auto t1_view = (*t1_amps_)(b, All);
                tensor_algebra::einsum(1.0, Indices { index::cd }, &out_view, (double) vv_mults(a, b), Indices { index::cd, index::i },
                        W_view, Indices { index::i }, t1_view);
            }
        }

        // <ai||jb>
        // Contract with non-transformed <ai||be>
        for (int e = 0; e < nvirt_; e++) {
            auto t1_view = (*t1_amps_)(e, All);
            for (int b = 0; b < nvirt_; b++) {
                if (vv_mults(b, e) == 0) {
                    continue;
                }
                auto W_view = (*integrals_.full->W_a_x_bc_i)(All, vv_inds(b, e), All);
                auto out_view = (*t1_integrals_.full->W_a_j_b_i)(All, All, b, All);
                tensor_algebra::einsum(1.0, Indices { index::a, index::j, index::i }, &out_view, (double) vv_mults(b, e), Indices {
                        index::a, index::i }, W_view, Indices { index::j }, t1_view);
            }
        }

        // Contract with transformed <im||jb>
        for(int m = 0; m < nocc_; m++) {
            auto t1_view = (*t1_amps_)(All, m);
        for(int i = 0; i < nocc_; i++) {
                if(oo_mults(i, m) == 0) {
                    continue;
                }
                auto W_view = (*t1_integrals_.full->W_x_k_a_ij)(All, oo_inds(i, m), All);
                auto out_view = (*t1_integrals_.full->W_a_j_b_i)(All, All, All, i);

                tensor_algebra::einsum(1.0, Indices{index::a, index::j, index::b}, &out_view, (double) oo_mults(i, m), Indices{index::})
            }
        }
    }

}

}

