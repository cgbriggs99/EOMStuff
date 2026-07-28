/*
 * Wabcd.cpp
 *
 *  Created on: Jul 27, 2026
 *      Author: connor
 */

#include <general_cc/general-cc.hpp>

using namespace einsums;

namespace psi::general_cc {

void GeneralCC::contract_Wabcd(diagram::triplet const &amps_spec, DiskTensor<double, 4> *out) {

    if (integral_type_ == FULL) {

        if (t1_transformed_) {
            contract_with_intermediate(amps_spec, *t1_integrals_.full->W_ab_x_cd_x, out);
        } else {
            contract_with_intermediate(amps_spec, *integrals_.full->W_ab_x_cd_x, out);
        }
    } else {

        Tensor<double, 3> *B_Qab;

        size_t naux;
        int order = std::get < 0 > (amps_spec);

        if (t1_transformed_) {
            B_Qab = t1_integrals_.df->B_Qab.get();
            naux = t1_integrals_.df->naux;
        } else {
            B_Qab = integrals_.df->B_Qab.get();
            naux = integrals_.df->naux;
        }
        if (std::get < 0 > (amps_spec) == 1) {

            size_t const num_Ao = OrbitString::number_of_strings(2, nvirt_);

            Tensor<double, 2> W_cd { "W cd", nvirt_, nvirt_ };

            for (size_t a = 0, Ao = 0; a < nvirt_; a++) {
                auto B_view_1 = (*B_Qab)(All, a, All);
                for (size_t b = 0; b < a; b++, Ao++) {
                    auto B_view_2 = (*B_Qab)(All, b, All);

                    tensor_algebra::einsum(Indices { index::c, index::d }, &W_cd, Indices { index::Q, index::c }, B_view_1, Indices {
                            index::Q, index::d }, B_view_2);

                    auto out_view = (*out)(Ao, All, All, 0);
                    auto &out_tensor = out_view.get();

                    tensor_algebra::einsum(1.0, Indices { index::C, index::I }, &out_tensor, 1.0, Indices { index::B, index::C }, W_cd,
                            Indices { index::B, index::I }, *t1_amps_);
                    tensor_algebra::einsum(1.0, Indices { index::C, index::I }, &out_tensor, -1.0, Indices { index::C, index::B }, W_cd,
                            Indices { index::B, index::I }, *t1_amps_);
                }
            }
        } else {
            size_t num_At, num_I, num_B, num_Ao, num_C;

            find_loop_parameters(amps_spec, 2, 0, 2 - std::get < 1 > (amps_spec), 0, nullptr, nullptr, nullptr, nullptr, nullptr, &num_At,
                    &num_Ao, &num_I, nullptr, &num_B, nullptr, &num_C, nullptr);

            auto const[D_mult, D_inds] = get_product_table(std::get < 1 > (amps_spec), order - std::get < 1 > (amps_spec), nvirt_);
            auto const[An_mult, An_inds] = get_product_table(2, order - std::get < 1 > (amps_spec), nvirt_);

            Tensor<double, 2> W_cd { "W cd", nvirt_, nvirt_ };
            Tensor<double, 1> temp_v { "temp v", num_B };

            std::unique_ptr<DiskView<double, 2>> t_cache { nullptr };
            std::unique_ptr<DiskView<double, 4>> out_cache { nullptr };
            Dim < 2 > cache_pos;
            Dim < 4 > out_cache_pos;

            for (size_t At = 0; At < num_At; At++) {
                for (size_t I = 0; I < num_I; I++) {
                    temp_v.zero();
                    size_t I_view = I % block_size;

                    for (size_t B = 0; B < num_B; B++) {
                        if (D_mult(B, At) == 0) {
                            continue;
                        }

                        size_t D = D_inds(B, At);
                        size_t view_D = D % block_size;

                        cache_view(&t_cache, &cache_pos, *tn_amps_[order - 2], Dim<2> { D, I });

                        temp_v(B) = D_mult(B, At) * (*t_cache)(D_view, I_view);
                    }

                    size_t Ao = 0;
                    for (int a = 0; a < nvirt_; a++) {
                        auto a_view = (*B_Qab)(All, a, All);
                        for (int b = 0; b < a; b++, Ao++) {
                            if (An_mults(Ao, At) == 0) {
                                continue;
                            }

                            auto b_view = (*B_Qab)(All, b, All);
                            // Compute the block of the <ab||cd> tensor.
                            tensor_algebra::einsum(Indices { index::c, index::d }, &W_cd, Indices { index::Q, index::c }, a_view, Indices {
                                    index::Q, index::d }, b_view);
                            tensor_algebra::einsum(1.0, Indices { index::c, index::d }, &W_cd, -1.0, Indices { index::Q, index::d }, a_view,
                                    Indices { index::Q, index::c }, b_view);

                            size_t An = An_inds(Ao, At);
                            size_t An_view = An % block_size;

                            // Contract.



                        }
                    }
                }
            }
        }
    }
}

}

