/*
 * Fij.cpp
 *
 *  Created on: Jul 27, 2026
 *      Author: connor
 */

#include <general_cc/general-cc.hpp>

using namespace einsums;

namespace psi::general_cc {
void GeneralCC::contract_Fij(diagram::triplet const &amps_spec, einsums::DiskTensor<double, 4> *out) const {

    Tensor<double 2> const *F;

    // Select which F to use.
    if (t1_transformed_) {
        F = t1_F_ij_.get();
    } else {
        F = F_ij_.get();
    }

    // The only pattern we can use is X11.
    int const order = std::get < 0 > (amps_spec);
    if (std::get < 1 > (amps_spec) != 1 || std::get < 2 > (amps_spec) != 0) {
        EINSUMS_THROW_EXCEPTION(std::runtime_error, "The Fij term needs to only have one internal hole line.");
    }

    // Handle T1 separately.
    if (order == 1) {
        auto out_view = (*out)(All, All, 0, 0);

        auto &out_tensor = out_view.get();

        tensor_algebra::einsum(1.0, Indices { index::a, index::i }, &out_tensor, 1.0, Indices { index::m, index::i }, *F, Indices {
                index::a, index::m }, *t1_amps_);
    } else {

        Tensor<double, 1> temp_o { "Temp o", nocc_ };

        auto const[occ_mult, occ_inds] = get_product_table(1, order - 1, nocc_);

        size_t const num_It = OrbitString::number_of_strings(order - 1, nocc_);
        size_t const num_A = OrbitString::number_of_strings(order, nvirt_);

        std::unique_ptr<DiskView<double, 2>> t_cache { nullptr };
        std::unique_ptr<DiskView<double, 4>> out_cache { nullptr };
        Dim < 2 > cache_pos;
        Dim < 4 > out_cache_pos;

        for (size_t A = 0; A < num_A; A++) {
            size_t const view_A = A % block_size;
            for (size_t It = 0; It < num_It; It++) {
                temp_o.zero();

                for (size_t J = 0; J < nocc_; J++) {
                    if (occ_mult(J, It) == 0) {
                        continue;
                    }

                    size_t const L = occ_inds(J, It);
                    size_t const view_L = L % block_size;

                    // Cache the T tensor.
                    cache_view(&t_cache, &cache_pos, *tn_amps_[order - 2], Dim<2> { L, A });

                    auto const &t_tensor = t_cache->get();

                    // Compute the element.
                    temp_o(J) = virt_mult(J, It) * t_tensor(A, L);
                }

                for (size_t Io = 0; Io < nocc_; Io++) {
                    if (occ_mult(Io, It) == 0) {
                        continue;
                    }
                    size_t const In = occ_inds(Io, It);

                    size_t const In_view = In % block_size;

                    cache_view(&out_cache, &out_cache_pos, *out, Dim<4> { A, In, 0, 0 });

                    auto &out_tens = out_cache->get();

                    tensor_algebra::einsum(1.0, Indices { }, out_tens.data(A_view, In_view, 0, 0), static_cast<double>(occ_mult(Io, It)),
                            Indices { index::J }, temp_o, Indices { index::J }, (*F)(Io, All));

                }
            }
        }
    }
}
}

