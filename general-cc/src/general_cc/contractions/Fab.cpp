/*
 * Fab.cpp
 *
 *  Created on: Jul 27, 2026
 *      Author: connor
 */

#include <general_cc/general-cc.hpp>

using namespace einsums;

namespace psi::general_cc {
void GeneralCC::contract_Fab(diagram::triplet const &amps_spec, einsums::DiskTensor<double, 4> *out) const {

    Tensor<double 2> const *F;

    // Select which F to use.
    if (t1_transformed_) {
        F = t1_F_ab_.get();
    } else {
        F = F_ab_.get();
    }

    // The only pattern we can use is X11.
    int const order = std::get < 0 > (amps_spec);
    if (std::get < 1 > (amps_spec) != 1 || std::get < 2 > (amps_spec) != 1) {
        EINSUMS_THROW_EXCEPTION(std::runtime_error, "The Fab term needs to only have one internal particle line.");
    }

    // Handle T1 separately.
    if (order == 1) {
        auto out_view = (*out)(All, All, 0, 0);

        auto &out_tensor = out_view.get();

        tensor_algebra::einsum(1.0, Indices { index::a, index::i }, &out_tensor, 1.0, Indices { index::a, index::e }, *F, Indices {
                index::e, index::i }, *t1_amps_);
    } else {

        Tensor<double, 1> temp_v { "Temp v", nvirt_ };

        auto const[virt_mult, virt_inds] = get_product_table(1, order - 1, nvirt_);

        size_t const num_I = OrbitString::number_of_strings(order, nocc_);
        size_t const num_At = OrbitString::number_of_strings(order - 1, nvirt_);

        std::unique_ptr<DiskView<double, 2>> t_cache { nullptr };
        std::unique_ptr<DiskView<double, 4>> out_cache { nullptr };
        Dim < 2 > cache_pos;
        Dim < 4 > out_cache_pos;

        for (size_t I = 0; I < num_I; I++) {
            size_t view_I = I % block_size;
            for (size_t At = 0; At < num_At; At++) {
                temp_v.zero();

                for (size_t B = 0; B < nvirt_; B++) {
                    if (virt_mult(B, At) == 0) {
                        continue;
                    }

                    size_t D = virt_inds(B, At);
                    size_t view_D = D % block_size;

                    // Cache the T tensor.
                    cache_view(&t_cache, &cache_pos, *tn_amps_[order - 2], Dim<2> { I, D });

                    auto const &t_tensor = t_cache->get();

                    // Compute the element.
                    temp_v(B) = virt_mult(B, At) * t_tensor(D, I);
                }

                for (size_t Ao = 0; Ao < nvirt_; Ao++) {
                    if (virt_mult(Ao, At) == 0) {
                        continue;
                    }
                    size_t const An = virt_inds(Ao, At);

                    size_t const An_view = An % block_size;

                    cache_view(&out_cache, &out_cache_pos, *out, Dim<4> { An, I, 0, 0 });

                    auto &out_tens = out_cache->get();

                    tensor_algebra::einsum(1.0, Indices { }, out_tens.data(An_view, I_view, 0, 0), static_cast<double>(virt_mult(Ao, At)),
                            Indices { index::B }, temp_v, Indices { index::B }, (*F)(All, Ao));

                }
            }
        }
    }
}
}
