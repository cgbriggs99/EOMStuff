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
#include <list>

// My editor is complaining. This should be a no-op in real contexts.
#ifndef PSI_API
#define PSI_API
#endif

namespace psi {
namespace general_cc {

// How large the blocks should be when reading in disk tensors.
constexpr int block_size = 512;

namespace detail {
struct PSI_API FullIntegrals {
    // Optimal ordering seems to be virtual out, occupied out, virtual sum, occupied sum.
    // Pairs of letters have the restriction that the first is less than the second.
    std::shared_ptr<einsums::DiskTensor<double, 4>> W_ab_ij_x_x, W_ab_x_cd_x, W_x_kl_x_ij, W_x_x_ab_ij;
    std::shared_ptr<einsums::DiskTensor<double, 4>> W_ab_i_c_x, W_a_jk_x_i, W_a_x_bc_i, W_x_k_a_ij;
    std::shared_ptr<einsums::DiskTensor<double, 4>> W_a_j_b_i;

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
    std::shared_ptr<einsums::Tensor<double, 3>> B_Qab, B_Qai, B_Qia, B_Qij;
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

    void find_loop_parameters(diagram::triplet const &amps_spec, int Ao_inds, int Io_inds, int C_inds, int K_inds, int *At_inds,
            int *It_inds, int *B_inds, int *J_inds, size_t *num_At, size_t *num_Ao, size_t *num_It, size_t *num_Io, size_t *num_B,
            size_t *num_J, size_t *num_C, size_t *num_K) const;

    void t1_transform();

    void fock_ao_to_so();

    void tei_ao_to_antisym_so(std::shared_ptr<MintsHelper> mintshelper);

    void contract_Fab(diagram::triplet const &amps_spec, einsums::DiskTensor<double, 4> *out) const;

    void contract_Fij(diagram::triplet const &amps_spec, einsums::DiskTensor<double, 4> *out) const;

    void contract_Fia(diagram::triplet const &amps_spec, einsums::DiskTensor<double, 4> *out) const;

    void contract_Wabcd(diagram::triplet const &amps_spec, einsums::DiskTensor<double, 4> *out) const;

    void contract_Wijkl(diagram::triplet const &amps_spec, einsums::DiskTensor<double, 4> *out) const;

    void contract_Wiajb(diagram::triplet const &amps_spec, einsums::DiskTensor<double, 4> *out) const;

    void contract_Wiabc(diagram::triplet const &amps_spec, einsums::DiskTensor<double, 4> *out) const;

    void contract_Wijak(diagram::triplet const &amps_spec, einsums::DiskTensor<double, 4> *out) const;

    void contract_Wabic(diagram::triplet const &amps_spec, einsums::DiskTensor<double, 4> *out) const;

    void contract_Waijk(diagram::triplet const &amps_spec, einsums::DiskTensor<double, 4> *out) const;

    void contract_Wijab(diagram::triplet const &amps_spec, einsums::DiskTensor<double, 4> *out) const;

    void contract_Fai(einsums::DiskTensor<double, 4> *out) const;

    void contract_Wabij(einsums::DiskTensor<double, 4> *out) const;

    void contract_interaction(diagram::triplet const &amps_spec, diagram::Interaction interaction,
            einsums::DiskTensor<double, 4> *out) const;

    void contract_with_intermediate(diagram::triplet const &amps_spec, einsums::DiskTensor<double, 4> const &intermediate,
            einsums::DiskTensor<double, 4> *out, int Ao_inds, int Io_inds, int C_inds, int K_inds) const;

    void contract(diagram::FactoredDiagram const &diagrams, einsums::DiskTensor<double, 2> *out) const;

    double energy_contract(diagram::Diagram const &diagram) const;

    void resolvent_contract(int num_inds, einsums::DiskTensor<double, 2> *residual) const;

    static std::pair<einsums::Tensor<signed char, 2>, einsums::Tensor<ptrdiff_t, 2>> get_product_table(int nleft, int nright, int norbs);

    diagram::Theory theory_;

    std::shared_ptr<einsums::Tensor<double, 2>> Ca_, Cb_, Fa_, Fb_;

    enum IntegralTypes {
        FULL, DENSITY_FITTED
    } integral_type_;

    // Get the Fock terms out of the way.
    std::shared_ptr<einsums::Tensor<double, 2>> F_ab_, F_ij_, F_ia_, F_ai_;

    std::shared_ptr<einsums::Tensor<double, 2>> t1_F_ab_, t1_F_ij_, t1_F_ia_, t1_F_ai_;

    std::shared_ptr<einsums::Tensor<double, 1>> mo_energies_a_, mo_energies_b_, mo_energies_o_, mo_energies_v_;

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

    h5::fd_t temp_file_;

    std::vector<h5::fd_t> tn_files_;
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

namespace detail {
template <size_t Rank>
void cache_view(std::unique_ptr<einsums::DiskTensorView<double, Rank>> *current, einsums::Dim<Rank> *cache_offset,
        einsums::DiskTensor<double, 2> &parent, einsums::Dim<Rank> const &position) {
    einsums::Dim<Rank> begin, end;

#pragma unroll
    for (size_t i = 0; i < Rank; i++) {
        begin[i] = block_size * (position[i] / block_size);
        end[i] = std::min(parent.dim(i), begin[i] + position[i] % block_size);
    }

    bool update_cache = !*current;

    if (cache_offset != nullptr) {
#pragma unroll
        for (size_t i = 0; i < Rank; i++) {
            update_cache = update_cache || (*cache_offset)[i] != begin[i];
        }
    }

    if (update_cache) {
        std::array < einsums::Range, Rank > ranges;

#pragma unroll
        for (size_t i = 0; i < Rank; i++) {
            ranges[i] = Range { begin[i], end[i] };
        }

        *current = std::make_unique<einsums::DiskView<double, 2>>(std::move(std::apply(parent, ranges)));

    }
}

template <size_t Rank>
void cache_view(std::unique_ptr<einsums::DiskTensorView<double, Rank>> const *current, einsums::Dim<Rank> *cache_offset,
        einsums::DiskTensor<double, 2> const &parent, einsums::Dim<Rank> const &position) {
    einsums::Dim<Rank> begin, end;

#pragma unroll
    for (size_t i = 0; i < Rank; i++) {
        begin[i] = block_size * (position[i] / block_size);
        end[i] = std::min(parent.dim(i), begin[i] + position[i] % block_size);
    }

    bool update_cache = !*current;

    if (cache_offset != nullptr) {
#pragma unroll
        for (size_t i = 0; i < Rank; i++) {
            update_cache = update_cache || (*cache_offset)[i] != begin[i];
        }
    }

    if (update_cache) {
        std::array < einsums::Range, Rank > ranges;

#pragma unroll
        for (size_t i = 0; i < Rank; i++) {
            ranges[i] = Range { begin[i], end[i] };
        }

        *current = std::make_unique<einsums::DiskView<double, 2> const>(std::move(std::apply(parent, ranges)));

    }
}
}

}
}

#endif /* GENERAL_CC_INCLUDE_GENERAL_CC_HPP_ */
