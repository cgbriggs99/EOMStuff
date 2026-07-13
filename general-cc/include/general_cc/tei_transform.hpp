/*
 * tei_transform.hpp
 *
 *  Created on: Jul 9, 2026
 *      Author: connor
 */

#ifndef GENERAL_CC_INCLUDE_GENERAL_CC_TEI_TRANSFORM_HPP_
#define GENERAL_CC_INCLUDE_GENERAL_CC_TEI_TRANSFORM_HPP_

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.hpp"

#include <Einsums/Tensor/Tensor.hpp>
#include <Einsums/TensorAlgebra.hpp>

// My editor is complaining. This should be a no-op in real contexts.
#ifndef PSI_API
#define PSI_API
#endif

namespace psi::general_cc {
class PSI_API ConvTEIMethod {
public:
    ConvTEIMethod(std::shared_ptr<einsums::Tensor<double, 4>> TEI);

    template <einsums::MatrixConcept C1Type, einsums::MatrixConcept C2Type, einsums::MatrixConcept C3Type, einsums::MatrixConcept C4Type,
            einsums::TensorConcept OutType>
    void contract_tei(C1Type const &C1, C2Type const &C2, C3Type const &C3, C4Type const &C4, OutType *out,
            einsums::Tensor<double, 4> &TEI_temp, einsums::Tensor<double, 4> &TEI_temp2) const {
        //    int nbf = C1.dim(0);
        //
        //    auto TEI_temp = std::make_shared<einsums::Tensor<double, 4>>(
        //        "TEI temp", nbf, nbf, nbf, C4.dim(1));
        //    auto TEI_temp2 = std::make_shared<einsums::Tensor<double, 4>>(
        //        "TEI temp2", nbf, nbf, C4.dim(1), nbf);
        // Convert TEI to MO basis and antisymmetrize.

        // First transform.
        //    psi::outfile->Printf("Performing integral transformation...\n");
        //    psi::outfile->Printf("Performing first transformation...\n");
        auto view1 = TEI_temp(einsums::All, einsums::All, einsums::All, einsums::Range { 0, C4.dim(1) });
        einsums::tensor_algebra::einsum(
                einsums::Indices { einsums::index::mu, einsums::index::nu, einsums::index::lambda, einsums::index::s }, &view1,
                einsums::Indices { einsums::index::mu, einsums::index::nu, einsums::index::lambda, einsums::index::sigma }, *TEI_,
                einsums::Indices { einsums::index::sigma, einsums::index::s }, C4);

        //    psi::outfile->Printf("Performing first permutation...\n");

        auto view2 = TEI_temp2(einsums::All, einsums::All, einsums::Range { 0, C4.dim(1) }, einsums::All);
        einsums::tensor_algebra::permute(0.0, einsums::Indices { einsums::index::mu, einsums::index::nu, einsums::index::s,
                einsums::index::lambda }, &view2, 1.0, einsums::Indices { einsums::index::mu, einsums::index::nu, einsums::index::lambda,
                einsums::index::s }, view1);

        //    psi::outfile->Printf("Performing second transformation...\n");
        auto view3 = TEI_temp(einsums::All, einsums::All, einsums::Range { 0, C4.dim(1) }, einsums::Range { 0, C3.dim(1) });

        einsums::tensor_algebra::einsum(einsums::Indices { einsums::index::mu, einsums::index::nu, einsums::index::s, einsums::index::r },
                &view3, einsums::Indices { einsums::index::mu, einsums::index::nu, einsums::index::s, einsums::index::lambda }, view2,
                einsums::Indices { einsums::index::lambda, einsums::index::r }, C3);

        //    psi::outfile->Printf("Performing second permutation...\n");
        auto view4 = TEI_temp2(einsums::All, einsums::Range { 0, C3.dim(1) }, einsums::Range { 0, C4.dim(1) }, einsums::All);

        einsums::tensor_algebra::permute(einsums::Indices { einsums::index::mu, einsums::index::r, einsums::index::s, einsums::index::nu },
                &view4, einsums::Indices { einsums::index::mu, einsums::index::nu, einsums::index::s, einsums::index::r }, view3);

        //    psi::outfile->Printf("Performing third transformation...\n");
        auto view5 = TEI_temp(einsums::All, einsums::Range { 0, C3.dim(1) }, einsums::Range { 0, C4.dim(1) },
                einsums::Range { 0, C2.dim(1) });

        einsums::tensor_algebra::einsum(einsums::Indices { einsums::index::mu, einsums::index::r, einsums::index::s, einsums::index::q },
                &view5, einsums::Indices { einsums::index::mu, einsums::index::r, einsums::index::s, einsums::index::nu }, view4,
                einsums::Indices { einsums::index::nu, einsums::index::q }, C2);

        //    psi::outfile->Printf("Performing third permutation...\n");
        auto view6 = TEI_temp2(einsums::Range { 0, C4.dim(1) }, einsums::Range { 0, C3.dim(1) }, einsums::Range { 0, C2.dim(1) },
                einsums::All);

        einsums::tensor_algebra::permute(einsums::Indices { einsums::index::s, einsums::index::r, einsums::index::q, einsums::index::mu },
                &view6, einsums::Indices { einsums::index::mu, einsums::index::r, einsums::index::s, einsums::index::q }, view5);

        //    psi::outfile->Printf("Performing fourth transformation...\n");
        auto view7 = TEI_temp(einsums::Range { 0, C4.dim(1) }, einsums::Range { 0, C3.dim(1) }, einsums::Range { 0, C2.dim(1) },
                einsums::Range { 0, C1.dim(1) });

        einsums::tensor_algebra::einsum(einsums::Indices { einsums::index::s, einsums::index::r, einsums::index::q, einsums::index::p },
                &view7, einsums::Indices { einsums::index::s, einsums::index::r, einsums::index::q, einsums::index::mu }, view6,
                einsums::Indices { einsums::index::mu, einsums::index::p }, C1);

        //    psi::outfile->Printf("Performing fourth permutation...\n");

        einsums::tensor_algebra::permute(einsums::Indices { einsums::index::p, einsums::index::r, einsums::index::q, einsums::index::s },
                out, einsums::Indices { einsums::index::s, einsums::index::r, einsums::index::q, einsums::index::p }, view7);

        //    psi::outfile->Printf("Finished integral transformations.\n");
    }

    /**
     * Contracts the two-electron integrals from chemists' notation in the AO
     * basis into physicists' notation in the MO basis.
     *
     * @param C1 The coefficients for the first index.
     * @param C2 The coefficients for the second index in the input, third index
     * in the output.
     * @param C3 The coefficients for the third index in the input, second index
     * in the output.
     * @param C4 The coefficients for the fourth index.
     */
    template <einsums::MatrixConcept C1Type, einsums::MatrixConcept C2Type, einsums::MatrixConcept C3Type, einsums::MatrixConcept C4Type>
    std::shared_ptr<einsums::Tensor<double, 4>> contract_tei(C1Type const &C1, C2Type const &C2, C3Type const &C3, C4Type const &C4) const {
        int nbf = C1.dim(0);

        auto TEI_temp = std::make_shared<einsums::Tensor<double, 4>>("TEI temp", nbf, nbf, nbf, nbf);
        auto TEI_temp2 = std::make_shared<einsums::Tensor<double, 4>>("TEI temp2", nbf, nbf, nbf, nbf);
        auto out = std::make_shared<einsums::Tensor<double, 4>>("out", C1.dim(1), C3.dim(1), C2.dim(1), C4.dim(1));
        contract_tei(C1, C2, C3, C4, out.get(), *TEI_temp, *TEI_temp2);

        return out;
    }

    std::shared_ptr<einsums::Tensor<double, 4>> get_TEI();

    std::shared_ptr<einsums::Tensor<double, 4> const> get_TEI() const;

private:
    std::shared_ptr<einsums::Tensor<double, 4>> TEI_;
};

class PSI_API DFMethod {
public:
    DFMethod(std::shared_ptr<einsums::Tensor<double, 3>> B);

    template <einsums::MatrixConcept C1Type, einsums::MatrixConcept C2Type>
    std::shared_ptr<einsums::Tensor<double, 3>>
    contract_B(C1Type const &C1, C2Type const &C2) const {
        using namespace einsums;
        using namespace einsums::tensor_algebra;
        auto temp1 = std::make_shared<Tensor<double, 3>>("temp1", B_->dim(0),
                C1.dim(1), B_->dim(2));
        auto temp2 = std::make_shared<Tensor<double, 3>>(
                "MO Density Fitting Tensor", B_->dim(0), C1.dim(1), C2.dim(1));

        temp1->zero();
        temp2->zero();

        einsum(Indices {index::Q, index::p, index::nu}, temp1.get(),
                Indices {index::Q, index::mu, index::nu}, *B_,
                Indices {index::mu, index::p}, C1);
        einsum(Indices {index::Q, index::p, index::q}, temp2.get(),
                Indices {index::Q, index::p, index::nu}, *temp1,
                Indices {index::nu, index::q}, C2);

        return temp2;
    }

    std::shared_ptr<einsums::Tensor<double, 3>> get_B();

    std::shared_ptr<einsums::Tensor<double, 3> const> get_B() const;

private:
    std::shared_ptr<einsums::Tensor<double, 3>> B_;
    mutable einsums::Tensor<double, 3>
    temp3d_; // mutable so that the storage space can be reused without issue.
    mutable einsums::Tensor<double, 1> temp1d_;
};

}

#endif /* GENERAL_CC_INCLUDE_GENERAL_CC_TEI_TRANSFORM_HPP_ */
