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
		// **SKIP square test because we only all this function
		//   in lme4 when we want to treat the input as a rectangular
		//   matrix
		// (!forceRectangularmatrix.rows() == matrix.cols()) ?
		// viewAsCholmod(matrix.template selfadjointView<_UpLo>()) :
		viewAsCholmod(matrix);

            cholmod_factorize_p(&A, &beta, fset.data(), fset.size(),
                                factor(), &cholmod());

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
            cholmod_dense* x_cd = cholmod_solve(type, factor(), &b_cd,
                                                &cholmod());
            if(!x_cd) {
                this->m_info = Eigen::NumericalIssue;
            }
	    typename Base::Scalar* xpt =
		reinterpret_cast<typename Base::Scalar*>(x_cd->x);
            std::copy(xpt, xpt + other.rows() * other.cols(), other.data());
            cholmod_free_dense(&x_cd, &cholmod());
        }
    };

    template<typename T>
    SEXP Eigen_cholmod_wrap(const lme4CholmodDecomposition<Eigen::SparseMatrix<T> >& obj) {
        typedef T* Tpt;
        const cholmod_factor* f = obj.factor();
        if (f->minor < f->n)
            throw std::runtime_error("CHOLMOD factorization was unsuccessful");

        //FIXME: Should extend this selection according to T
        ::Rcpp::S4 ans(std::string(f->is_super ? "dCHMsuper" : "dCHMsimpl"));
        ::Rcpp::IntegerVector  dd(2);
        dd[0] = dd[1] = f->n;
        ans.slot("Dim") = dd;
        ans.slot("perm") = ::Rcpp::wrap((int*)f->Perm, (int*)f->Perm + f->n);
        ans.slot("colcount") = ::Rcpp::wrap((int*)f->ColCount, (int*)f->ColCount + f->n);
        ::Rcpp::IntegerVector tt(f->is_super ? 6 : 4);
        tt[0] = f->ordering; tt[1] = f->is_ll;
        tt[2] = f->is_super; tt[3] = f->is_monotonic;
        ans.slot("type") = tt;
        if (f->is_super) {
            tt[4] = f->maxcsize; tt[5] = f->maxesize;
            ans.slot("super") = ::Rcpp::wrap((int*)f->super, ((int*)f->super) + f->nsuper + 1);
            ans.slot("pi")    = ::Rcpp::wrap((int*)f->pi, ((int*)f->pi) + f->nsuper + 1);
            ans.slot("px")    = ::Rcpp::wrap((int*)f->px, ((int*)f->px) + f->nsuper + 1);
            ans.slot("s")     = ::Rcpp::wrap((int*)f->s, ((int*)f->s) + f->ssize);
            ans.slot("x")     = ::Rcpp::wrap((Tpt)f->x, ((T*)f->x) + f->xsize);
        } else {
            ans.slot("i")     = ::Rcpp::wrap((int*)f->i, ((int*)f->i) + f->nzmax);
            ans.slot("p")     = ::Rcpp::wrap((int*)f->p, ((int*)f->p) + f->n + 1);
            ans.slot("x")     = ::Rcpp::wrap((Tpt)f->x, ((T*)f->x) + f->nzmax);
            ans.slot("nz")    = ::Rcpp::wrap((int*)f->nz, ((int*)f->nz) + f->n);
            ans.slot("nxt")   = ::Rcpp::wrap((int*)f->next, ((int*)f->next) + f->n + 2);
            ans.slot("prv")   = ::Rcpp::wrap((int*)f->prev, ((int*)f->prev) + f->n + 2);
        }
        return ::Rcpp::wrap(ans);
    }

} // namespace lme4

namespace Rcpp {

    template<typename T>
    SEXP wrap(const lme4::lme4CholmodDecomposition<Eigen::SparseMatrix<T> >& obj) {
        return ::lme4::Eigen_cholmod_wrap(obj);
    }

} // namespace Rcpp

#endif // LME4_CHOLMODDECOMPOSITION_H
