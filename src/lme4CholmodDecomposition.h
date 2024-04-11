// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// predModule.h: predictor module using Eigen
//
// Copyright (C)       2011 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#ifndef LME4_CHOLMODDECOMPOSITION_H
#define LME4_CHOLMODDECOMPOSITION_H

#include <RcppEigen.h>
#include <RcppEigenCholmod.h>

namespace lme4 {

    /* A wrapper around Eigen's CholmodDecomposition class that provides some extra
       functionality and member access required by lme4.
    */
    template<typename _MatrixType, int _UpLo = Eigen::Lower>
    class lme4CholmodDecomposition : public Eigen::CholmodDecomposition<_MatrixType, _UpLo> {
      protected:
        typedef Eigen::CholmodDecomposition<_MatrixType, _UpLo> Base;
        using Base::m_factorizationIsOk;
        using Base::m_analysisIsOk;

      public:
        cholmod_common& cholmod() const {
            return const_cast<lme4CholmodDecomposition<_MatrixType, _UpLo>*>(this)->Base::cholmod();
        }

        cholmod_factor* factor() const { return Base::m_cholmodFactor; }

        void factorize_p(const typename Base::MatrixType& matrix, Eigen::ArrayXi fset, double beta=0.) {
            // FIXME: add forceRectangular flag to allow backward compatibility;
            // restore an appropriate version of the square/rectangular test
            eigen_assert(m_analysisIsOk && "You must first call analyzePattern()");
            cholmod_sparse    A = 
                // **SKIP square test because we only call this function
                //   in lme4 when we want to treat the input as a rectangular
                //   matrix
                // (!forceRectangularmatrix.rows() == matrix.cols()) ?
                // viewAsCholmod(matrix.template selfadjointView<_UpLo>()) :
                viewAsCholmod(matrix);
            double beta_[2]; beta_[0] = beta; beta_[1] = 0.0;

            R_MATRIX_CHOLMOD(factorize_p)(
                &A, beta_, fset.data(), fset.size(), factor(), &cholmod());

            this->m_info = Eigen::Success;
            m_factorizationIsOk = true;
        }

        template<typename OtherDerived>
        void solveInPlace(const Eigen::MatrixBase<OtherDerived>& _other, int type) const {
            OtherDerived& other = _other.const_cast_derived();
            eigen_assert(m_factorizationIsOk &&
                         "The decomposition is not in a valid state for solving, you must first call either compute() or symbolic()/numeric()");
            eigen_assert((Base::Index)(factor()->n) == other.rows());

            // note: cd stands for Cholmod Dense
            cholmod_dense b_cd = viewAsCholmod(other.const_cast_derived());
            // m_cholmodFactor
            cholmod_dense* x_cd = R_MATRIX_CHOLMOD(solve)(
                type, factor(), &b_cd, &cholmod());
            if(!x_cd) {
                this->m_info = Eigen::NumericalIssue;
            }
            typename Base::Scalar* xpt =
                reinterpret_cast<typename Base::Scalar*>(x_cd->x);
            std::copy(xpt, xpt + other.rows() * other.cols(), other.data());
            R_MATRIX_CHOLMOD(free_dense)(&x_cd, &cholmod());
        }
    };

    template<typename T>
    SEXP Eigen_cholmod_wrap(const lme4CholmodDecomposition<Eigen::SparseMatrix<T> >& obj) {
        return M_chm_factor_to_SEXP(obj.factor(), 0);
    }

} // namespace lme4

namespace Rcpp {

    template<typename T>
    SEXP wrap(const lme4::lme4CholmodDecomposition<Eigen::SparseMatrix<T> >& obj) {
        return ::lme4::Eigen_cholmod_wrap(obj);
    }

} // namespace Rcpp

#endif // LME4_CHOLMODDECOMPOSITION_H
