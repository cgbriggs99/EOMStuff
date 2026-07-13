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
#include <Einsums/Tensor/DiskTensor.hpp>

#include <vector>

// My editor is complaining. This should be a no-op in real contexts.
#ifndef PSI_API
#define PSI_API
#endif

namespace psi {
namespace general_cc {

namespace detail {
struct PSI_API FullIntegrals {
    // Optimal ordering seems to be virtual out, occupied out, virtual sum, occupied sum.
    // Pairs of letters have the restriction that the first is less than the second.
    std::shared_ptr<einsums::Tensor<double, 2>> W_ab_ij_x_x, W_ab_x_cd_x, W_x_kl_x_ij, W_x_x_ab_ij;
    std::shared_ptr<einsums::Tensor<double, 3>> W_ab_i_c_x, W_a_jk_x_i, W_a_x_bc_i, W_x_k_a_ij;
    std::shared_ptr<einsums::Tensor<double, 4>> W_a_j_b_i;

    // Key:
    // W_ab_ij_x_x = <ab||ij>
    // W_ab_x_cd_x = <ab||cd>
    // W_x_kl_x_ij = <ij||kl>
    // W_x_x_ab_ij = <ij||ab>
    // W_ab_i_c_x = <ab||ci>
    // W_a_jk_x_i = <ia||jk>
    // W_a_x_bc_i = <ai||bc>
    // W_x_k_a_ij = <ij||ka>
    // W_a_j_b_i = <ai||jb>

    // Indices as given are sorted for optimality, not correspondence to the theory.
};

struct PSI_API DFIntegrals {
    // No restrictions on orbitals. This is because one index comes from the bra, one from the ket,
    // and Q is the auxiliary. Maybe I'll add other sorting forms of these.
    std::shared_ptr<einsums::Tensor<double, 3>> B_abQ, B_aiQ, B_iaQ, B_ijQ;
    size_t naux;
};

}

class PSI_API GeneralCC : public Wavefunction {
public:
    GeneralCC(diagram::Theory theory, SharedWavefunction ref_wfn, Options &options);

    virtual ~GeneralCC();

    virtual size_t memory_estimate() const;

    virtual double compute_iteration_energy() const;

    virtual double compute_energy();

    [[nodiscard]] inline double get_cc_energy() const {
        return cc_energy_;
    }

protected:

    void t1_transform();

    void fock_ao_to_so();

    void tei_ao_to_antisym_so(std::shared_ptr<MintsHelper> mintshelper);

    static std::pair<einsums::Tensor<signed char, 2>, einsums::Tensor<ptrdiff_t, 2>> get_product_table(int nleft, int nright,
            int norbs) const;

    diagram::Theory theory_;

    std::shared_ptr<einsums::Tensor<double, 2>> Ca_, Cb_, Fa_, Fb_;

    enum IntegralTypes {
        FULL, DENSITY_FITTED
    } integral_type_;

    // Get the Fock terms out of the way.
    std::shared_ptr<einsums::Tensor<double, 2>> F_ab_, F_ij_, F_ia_, F_ai_;

    std::shared_ptr<einsums::Tensor<double, 2>> t1_F_ab_, t1_F_ij_, t1_F_ia_, t1_F_ai_;

    std::shared_ptr<einsums::Tensor<double, 1>> mo_energies_a_, mo_energies_b_;

    union {
        detail::FullIntegrals *full { nullptr };
        detail::DFIntegrals *df;
    } integrals_, t1_integrals_;

    std::list<diagram::Diagram> energy_expression_;

    std::vector<diagram::FactoredDiagram> residual_expressions_;

    std::shared_ptr<einsums::Tensor<double, 2>> t1_amps_; // This one is a regular tensor because it's small and important.

    std::vector<std::shared_ptr<einsums::DiskTensor<double, 2>>> tn_amps_; // It's rank2 since we index it with strings.

    int nocc_, nvirt_, navirt_, nbvirt_;

    double cc_energy_;

    bool t1_transformed_ { false };
};

class PSI_API GeneralCCbracket : public GeneralCC {
public:
    GeneralCCbracket(diagram::Theory theory, SharedWavefunction ref_wfn, Options &options);

    virtual size_t memory_estimate() const;

    virtual double compute_iteration_energy() const;
    virtual double compute_energy();

    [[nodiscard]] inline double get_bracket_energy() const {
        return cc_bracket_energy_;
    }
protected:
    double cc_bracket_energy_;
};

class PSI_API GeneralCCparen : public GeneralCCbracket {
public:
    GeneralCCparen(diagram::Theory theory, SharedWavefunction ref_wfn, Options &options);

    virtual size_t memory_estimate() const;

    virtual double compute_iteration_energy() const;

    [[nodiscard]] inline double get_paren_energy() const {
        return cc_paren_energy_;
    }
protected:
    double cc_paren_energy_;
};

class PSI_API GeneralCCn : public GeneralCC {
public:
    GeneralCCn(diagram::Theory theory, SharedWavefunction ref_wfn, Options &options);

    virtual size_t memory_estimate() const;

    virtual double compute_energy();
};

}
}

#endif /* GENERAL_CC_INCLUDE_GENERAL_CC_HPP_ */
