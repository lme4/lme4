//
// predModule.cpp: implementation of predictor module using Eigen
//
// Copyright (C)       2011 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4Eigen.

#include "predModule.h"

using namespace Rcpp;
using namespace std;

namespace lme4Eigen {
    merPredD::merPredD(S4 X, S4 Zt, S4 Lambdat, IntegerVector Lind,
		       NumericVector theta)             throw (invalid_argument,
							       runtime_error)
	: d_X(X),
          d_Zt(Zt),
	  d_theta(clone(theta)),
	  d_Lind(Lind),
	  d_n(d_X.rows()),
	  d_p(d_X.cols()),
	  d_q(d_Zt.rows()),
	  d_Lambdat(Lambdat),
	  d_RZX(d_q, d_p),
	  d_V(d_X),
	  d_VtV(d_p,d_p),
	  d_Vtr(d_p),
	  d_Utr(d_q),
	  d_delb(d_p),
	  d_delu(d_q),
	  d_beta0(d_p),
	  d_u0(d_q),
	  d_Ut(d_Zt),
	  d_isDiagLam(false),
	  d_RX(d_p)
    {				// Check consistency of dimensions
	if (d_n != d_Zt.cols())
	    throw invalid_argument("Z dimension mismatch");
	if (d_Lind.size() != d_Lambdat.nonZeros())
	    throw invalid_argument("size of Lind does not match nonzeros in Lambda");
	// checking of the range of Lind is now done in R code for reference class

				// initialize beta0, u0 and VtV
	d_beta0.setZero(); d_u0.setZero();
	d_VtV.setZero().selfadjointView<Eigen::Upper>().rankUpdate(d_V.adjoint());
				// check for diagonal Lambdat
	if (d_Lambdat.nonZeros() == d_q) {
	    d_isDiagLam = true;
	    for (int j = 0; j < d_Lambdat.outerSize(); ++j) {
		SpMatrixXd::InnerIterator   it(d_Lambdat, j);
		if (it.index() != j)
		    throw invalid_argument("Lambdat is missing a diagonal element?");
	    }
	}
	if (d_isDiagLam) cout << "Diagonal Lambda detected" << endl;
				// starting values into Lambda
	setTheta(d_theta);
				// perform symbolic analysis
	SpMatrixXd   LamtUt(d_Lambdat * d_Ut);
	d_nnz = LamtUt.nonZeros();
cout << "LamtUt nonzeros = " << d_nnz << endl;
        d_L.setMode(Eigen::CholmodAutoLLt);      
	d_L.analyzePattern(LamtUt);
        if (d_L.info() != Eigen::Success)
	    throw runtime_error("CholeskyDecomposition.analyzePattern failed");
cout << "size of factor = " <<	d_L.rows()
     << ", is_ll = " << d_L.is_ll()
     << ", is_super = " << d_L.is_super()
     << ", ordering = " << d_L.ordering()
     << ", nzmax = " << d_L.nonZeros()
     << endl;

	const cholmod_factor* ff(d_L.factor());
        MVectorXi        Perm((int*)(ff->Perm), ff->n);
	MVectorXi    ColCount((int*)(ff->ColCount), ff->n);
cout << "Permutation: " << Perm.adjoint() << endl;
cout << "Column Counts: " << ColCount.adjoint() << endl;
    }

    void merPredD::updateL()                           throw (runtime_error) {
	if (d_isDiagLam) {
	    SpMatrixXd       LamtUt(d_Ut);
	    const Scalar*      dptr(d_Lambdat._valuePtr());
	    for (int j = 0; j < LamtUt.outerSize(); ++j) {
		for (SpMatrixXd::InnerIterator it(LamtUt, j); it; ++it)
		    it.valueRef() *= dptr[it.index()];
	    }
	    d_L.factorize_p(LamtUt, ArrayXi(), 1.);
	} else {		// do a lot of dancing around to avoid pruning zeros
	    const cholmod_sparse Lamt=viewAsCholmod(d_Lambdat), Ut=viewAsCholmod(d_Ut);
	    cholmod_common c = d_L.cholmod();
	    CHM_SP LamtUt = ::M_cholmod_ssmult(&Lamt, &Ut, 0/*stype*/, 1/*values*/, 0/*sorted*/, &c);
	    if (::M_cholmod_nnz(LamtUt, &c) != d_nnz)
		throw runtime_error("Number of nonzeros in LamtUt has changed");
	    d_L.factorize_p(LamtUt, ArrayXi(), 1.);
	    ::M_cholmod_free_sparse(&LamtUt, &c);
	}
	d_ldL2 = ::M_chm_factor_ldetL2(d_L.factor());
    }

    void merPredD::setTheta(const NumericVector& theta) throw (invalid_argument,
							       runtime_error) {
	if (theta.size() != d_theta.size())
	    throw invalid_argument("theta size mismatch");
				// update theta
	copy(theta.begin(), theta.end(), d_theta.begin());
				// update Lambdat
	int    *lipt = d_Lind.begin();
	double *LamX = d_Lambdat._valuePtr(), *thpt = d_theta.begin();
	for (int i = 0; i < d_Lind.size(); ++i) {
	    LamX[i] = thpt[lipt[i] - 1];
	}
    }

    void merPredD::solve() {
	d_L.setSolveType(CHOLMOD_P);
	d_delu           = d_L.solve(d_Utr);
	d_L.setSolveType(CHOLMOD_L);
	d_delu           = d_L.solve(d_delu);
				// d_delu now contains cu
	d_delb           = d_RX.solve(d_Vtr - d_RZX.adjoint() * d_delu);
	d_delu          -= d_RZX * d_delb;
	d_L.setSolveType(CHOLMOD_Lt);
	d_delu           = d_L.solve(d_delu);
	d_L.setSolveType(CHOLMOD_Pt);
	d_delu           = d_L.solve(d_delu);
    }

    void merPredD::updateXwts(const VectorType& sqrtXwt) throw (invalid_argument) {
	if (d_X.rows() != sqrtXwt.size())
	    throw invalid_argument("updateXwts: dimension mismatch");
	if (d_V.rows() == d_X.rows()) {  //FIXME: Generalize this for nlmer
	    DiagType  W(sqrtXwt.asDiagonal());
	    d_V         = W * d_X;
	    d_Ut        = d_Zt * W;
	} else throw invalid_argument("updateRes: no provision for nlmer yet");
	d_VtV.setZero().selfadjointView<Eigen::Upper>().rankUpdate(d_V.adjoint());
    }

    void merPredD::updateDecomp() {
				// update L, RZX and RX
	updateL();
	d_L.setSolveType(CHOLMOD_P);
	d_RZX           = d_L.solve(d_Lambdat * (d_Ut * d_V));
	d_L.setSolveType(CHOLMOD_L);
	d_RZX           = d_L.solve(d_RZX);

	MatrixXd      VtVdown(d_VtV);
	d_RX.compute(VtVdown.selfadjointView<Eigen::Upper>().rankUpdate(d_RZX.adjoint(), -1));
	if (d_RX.info() != Eigen::Success)
	    ::Rf_error("Downdated VtV is not positive definite");
	d_ldRX2         = 2. * d_RX.matrixLLT().diagonal().array().abs().log().sum();
    }

    void merPredD::updateRes(const VectorType& wtres)    throw (invalid_argument) {
	if (d_V.rows() != wtres.size())
	    throw invalid_argument("updateRes: dimension mismatch");
	d_Vtr           = d_V.adjoint() * wtres;
	d_Utr           = d_Lambdat * (d_Ut * wtres);
    }

    void merPredD::setBeta0(const VectorType& nBeta)     throw (invalid_argument) {
	if (nBeta.size() != d_p) throw invalid_argument("setBeta0: dimension mismatch");
	copy(nBeta.data(), nBeta.data() + d_p, d_beta0.data());
    }

    void merPredD::setU0(const VectorType& newU0)        throw (invalid_argument) {
	if (newU0.size() != d_q) throw invalid_argument("setU0: dimension mismatch");
	copy(newU0.data(), newU0.data() + d_q, d_u0.data());
    }

    IntegerVector merPredD::Pvec() const {
	const cholmod_factor* cf = d_L.factor();
	int*                 ppt = (int*)cf->Perm;
	return IntegerVector(ppt, ppt + cf->n);
    }
}

extern "C" {
    SEXP merPredDCreate(SEXP Xs, SEXP Zts, SEXP Lambdats, SEXP Linds, SEXP thetas) {
	BEGIN_RCPP;
	S4 X(Xs), Zt(Zts), Lambdat(Lambdats);
	IntegerVector Lind(Linds);
	NumericVector theta(thetas);
	lme4Eigen::merPredD *ans = new lme4Eigen::merPredD(X, Zt, Lambdat, Lind, theta);
	return wrap(XPtr<lme4Eigen::merPredD>(ans, true));
	END_RCPP;
    }

    SEXP merPredDsetBeta0(SEXP ptr, SEXP beta0) {
	BEGIN_RCPP;
	XPtr<lme4Eigen::merPredD>(ptr)->setBeta0(as<VectorXd>(beta0));
	END_RCPP;
    }
    
    SEXP merPredDsetTheta(SEXP ptr, SEXP theta) {
	BEGIN_RCPP;
	XPtr<lme4Eigen::merPredD>(ptr)->setTheta(NumericVector(theta));
	END_RCPP;
    }
    
    SEXP merPredDsetU0(SEXP ptr, SEXP u0) {
	BEGIN_RCPP;
	XPtr<lme4Eigen::merPredD>(ptr)->setU0(as<VectorXd>(u0));
	END_RCPP;
    }
    
    // SEXP merPredDI(SEXP ptr) {
    // 	BEGIN_RCPP;
    // 	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->I());
    // 	END_RCPP;
    // }
    
    SEXP merPredDLambdat(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->Lambdat());
	END_RCPP;
    }
    
    // SEXP merPredDL(SEXP ptr) {
    // 	BEGIN_RCPP;
    // 	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->L());
    // 	END_RCPP;
    // }
    
    SEXP merPredDPvec(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->Pvec());
	END_RCPP;
    }
    
    SEXP merPredDRX(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->RX());
	END_RCPP;
    }
    
    SEXP merPredDRXdiag(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->RXdiag());
	END_RCPP;
    }
    
    SEXP merPredDRZX(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->RZX());
	END_RCPP;
    }
    
    SEXP merPredDZt(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->Zt());
	END_RCPP;
    }
    
    SEXP merPredDVtV(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->VtV());
	END_RCPP;
    }
    
    SEXP merPredDbeta0(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->beta0());
	END_RCPP;
    }
    
    SEXP merPredDdelb(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->delb());
	END_RCPP;
    }
    
    SEXP merPredDdelu(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->delu());
	END_RCPP;
    }
    
    SEXP merPredDu0(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->u0());
	END_RCPP;
    }
    
    SEXP merPredDunsc(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->unsc());
	END_RCPP;
    }
    
    SEXP merPredDldL2(SEXP ptr) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lme4Eigen::merPredD>(ptr)->ldL2());
	END_RCPP;
    }
    
    SEXP merPredDldRX2(SEXP ptr) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lme4Eigen::merPredD>(ptr)->ldRX2());
	END_RCPP;
    }
    
    SEXP merPredDsqrL(SEXP ptr, SEXP fac) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lme4Eigen::merPredD>(ptr)->sqrL(::Rf_asReal(fac)));
	END_RCPP;
    }
    
    SEXP merPredDtheta(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->theta());
	END_RCPP;
    }
    
    SEXP merPredDinstallPars(SEXP ptr, SEXP fac) {
	BEGIN_RCPP;
	XPtr<lme4Eigen::merPredD>(ptr)->installPars(::Rf_asReal(fac));
	END_RCPP;
    }
    
    SEXP merPredDlinPred(SEXP ptr, SEXP fac) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->linPred(::Rf_asReal(fac)));
	END_RCPP;
    }
	    
}
