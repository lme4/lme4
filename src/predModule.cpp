//
// predModule.cpp: implementation of predictor module using Eigen
//
// Copyright (C)       2011 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4Eigen.

#include "predModule.h"

namespace lme4Eigen {
    merPredD::merPredD(S4 X, S4 Zt, S4 Lambdat, IntegerVector Lind,
		       NumericVector theta)
	: d_X(X),
          d_Zt(Zt),
	  d_theta(clone(theta)),
	  d_Lind(Lind),
	  d_n(d_X.rows()),
	  d_nnz(-1),
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
	  d_RX(d_p),
	  d_LamtUtRestructure(false)
    {				// Check consistency of dimensions
	if (d_n != d_Zt.cols())
	    throw invalid_argument("Z dimension mismatch");
	if (d_Lind.size() != d_Lambdat.nonZeros())
	    throw invalid_argument("size of Lind does not match nonzeros in Lambda");
	// checking of the range of Lind is now done in R code for reference class
				// initialize beta0, u0, delb, delu and VtV
	d_beta0.setZero(); d_u0.setZero(); d_delu.setZero(); d_delb.setZero();
	d_VtV.setZero().selfadjointView<Eigen::Upper>().rankUpdate(d_V.adjoint());
	d_RX.compute(d_VtV);	// ensure d_RX is initialized even in the 0-column X case

	setTheta(d_theta);		// starting values into Lambda
        d_L.cholmod().final_ll = 1;	// force an LL' decomposition
	d_LamtUt = d_Lambdat * d_Ut;
	d_L.analyzePattern(d_LamtUt); // perform symbolic analysis
        if (d_L.info() != Eigen::Success)
	    throw runtime_error("CholeskyDecomposition.analyzePattern failed");
    }

    void merPredD::updateLamtUt() {
	if (d_LamtUtRestructure) {
	    d_LamtUt                           = d_Lambdat * d_Ut;
	    return;
	}
	// This complicated code bypasses problems caused by Eigen's
	// sparse/sparse matrix multiplication pruning zeros.  The
	// Cholesky decomposition croaks if the structure of d_LamtUt changes.
	std::fill(d_LamtUt._valuePtr(), d_LamtUt._valuePtr() + d_LamtUt.nonZeros(), Scalar());
	for (Index j = 0; j < d_Ut.cols(); ++j) {
	    for(SpMatrixd::InnerIterator rhsIt(d_Ut, j); rhsIt; ++rhsIt) {
		SpMatrixd::InnerIterator prdIt(d_LamtUt, j);
		Scalar                       y = rhsIt.value();
		for (SpMatrixd::InnerIterator lhsIt(d_Lambdat, rhsIt.index()); lhsIt; ++lhsIt) {
		    Index                    i = lhsIt.index();
		    while (prdIt.index() != i && prdIt) ++prdIt;
		    if (!prdIt) throw runtime_error("logic error in updateLamtUt");
		    prdIt.valueRef()          += lhsIt.value() * y;
		}
	    }
	}
    }

    VectorXd merPredD::b(const double& f) const {return d_Lambdat.adjoint() * u(f);}

    VectorXd merPredD::beta(const double& f) const {return d_beta0 + f * d_delb;}

    VectorXd merPredD::linPred(const double& f) const {
	return d_X * beta(f) + d_Zt.adjoint() * b(f);
    }

    VectorXd merPredD::u(const double& f) const {return d_u0 + f * d_delu;}

    merPredD::Scalar merPredD::sqrL(const double& f) const {return u(f).squaredNorm();}

    void merPredD::updateL() {
	updateLamtUt();
	d_L.factorize_p(d_LamtUt, ArrayXi(), 1.);
	d_ldL2 = ::M_chm_factor_ldetL2(d_L.factor());
    }

    void merPredD::setTheta(const NumericVector& theta) {
	if (theta.size() != d_theta.size())
	    throw invalid_argument("theta size mismatch");
				// update theta
	std::copy(theta.begin(), theta.end(), d_theta.begin());
				// update Lambdat
	int    *lipt = d_Lind.begin();
	double *LamX = d_Lambdat._valuePtr(), *thpt = d_theta.begin();
	for (int i = 0; i < d_Lind.size(); ++i) {
	    LamX[i] = thpt[lipt[i] - 1];
	}
    }

    merPredD::Scalar merPredD::solve() {
	d_delu          = d_Utr - d_u0;
	d_L.solveInPlace(d_delu, CHOLMOD_P);
	d_L.solveInPlace(d_delu, CHOLMOD_L);    // d_delu now contains cu
	d_CcNumer       = d_delu.squaredNorm(); // numerator of convergence criterion

	d_delb          = d_RX.matrixL().solve(d_Vtr - d_RZX.adjoint() * d_delu);
	d_CcNumer      += d_delb.squaredNorm(); // increment CcNumer
	d_RX.matrixU().solveInPlace(d_delb);

	d_delu         -= d_RZX * d_delb;
	d_L.solveInPlace(d_delu, CHOLMOD_Lt);
	d_L.solveInPlace(d_delu, CHOLMOD_Pt);
	return d_CcNumer;
    }

    merPredD::Scalar merPredD::solveU() {
	d_delb.setZero(); // in calculation of linPred delb should be zero after solveU
	d_delu          = d_Utr - d_u0;
	d_L.solveInPlace(d_delu, CHOLMOD_P);
	d_L.solveInPlace(d_delu, CHOLMOD_L);    // d_delu now contains cu
	d_CcNumer       = d_delu.squaredNorm(); // numerator of convergence criterion
	d_L.solveInPlace(d_delu, CHOLMOD_Lt);
	d_L.solveInPlace(d_delu, CHOLMOD_Pt);
	return d_CcNumer;
    }

    void merPredD::updateXwts(const MatrixXd& sqrtXwt) {
	if (d_X.rows() != sqrtXwt.size())
	    throw invalid_argument("updateXwts: dimension mismatch");
	if (sqrtXwt.cols() == 1) {
	    DiagonalMatrix<double, Dynamic> W(sqrtXwt.col(0).asDiagonal());
	    d_V              = W * d_X;
	    d_Ut             = d_Zt * W;
	} else {
	    int n            = sqrtXwt.rows();
	    SpMatrixd      W(n, sqrtXwt.size());
	    const double *pt = sqrtXwt.data();
	    W.reserve(sqrtXwt.size());
	    for (Index j = 0; j < W.cols(); ++j, ++pt) {
		W.startVec(j);
		W.insertBack(j % n, j) = *pt;
	    }
	    W.finalize();
	    d_V              = W * d_X;
	    d_Ut             = d_Zt * W.adjoint();
	}
	if (d_LamtUt.rows() != d_Ut.rows() || 
	    d_LamtUt.cols() != d_Ut.cols()) d_LamtUtRestructure = true;
	d_VtV.setZero().selfadjointView<Upper>().rankUpdate(d_V.adjoint());
	updateL();
    }

    void merPredD::updateDecomp() { // update L, RZX and RX
	updateL();
	d_RZX           = d_LamtUt * d_V;
	if (d_p > 0) {
	    d_L.solveInPlace(d_RZX, CHOLMOD_P);
	    d_L.solveInPlace(d_RZX, CHOLMOD_L);

	    MatrixXd      VtVdown(d_VtV);
	    d_RX.compute(VtVdown.selfadjointView<Eigen::Upper>().rankUpdate(d_RZX.adjoint(), -1));
	    if (d_RX.info() != Eigen::Success)
		::Rf_error("Downdated VtV is not positive definite");
	    d_ldRX2         = 2. * d_RX.matrixLLT().diagonal().array().abs().log().sum();
	}
    }

    void merPredD::updateRes(const VectorXd& wtres) {
	if (d_V.rows() != wtres.size())
	    throw invalid_argument("updateRes: dimension mismatch");
	d_Vtr           = d_V.adjoint() * wtres;
	d_Utr           = d_LamtUt * wtres;
    }

    void merPredD::installPars(const Scalar& f) {
	d_u0 = u(f);
	d_beta0 = beta(f);
	d_delb.setZero();
	d_delu.setZero();
    }

    void merPredD::setBeta0(const VectorXd& nBeta) {
	if (nBeta.size() != d_p) throw invalid_argument("setBeta0: dimension mismatch");
	std::copy(nBeta.data(), nBeta.data() + d_p, d_beta0.data());
    }

    void merPredD::setU0(const VectorXd& newU0) {
	if (newU0.size() != d_q) throw invalid_argument("setU0: dimension mismatch");
	std::copy(newU0.data(), newU0.data() + d_q, d_u0.data());
    }

    IntegerVector merPredD::Pvec() const {
	const cholmod_factor* cf = d_L.factor();
	int*                 ppt = (int*)cf->Perm;
	return IntegerVector(ppt, ppt + cf->n);
    }

    MatrixXd merPredD::RX() const {
	return d_RX.matrixU();
    }

    MatrixXd merPredD::RXi() const {
	return d_RX.matrixU().solve(MatrixXd::Identity(d_p,d_p));
    }

    MatrixXd merPredD::unsc() const {
	return MatrixXd(MatrixXd(d_p, d_p).setZero().
			selfadjointView<Lower>().rankUpdate(RXi().adjoint()));
    }

    VectorXd merPredD::RXdiag() const {
	return d_RX.matrixLLT().diagonal();
    }
}
