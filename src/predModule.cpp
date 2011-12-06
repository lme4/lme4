//
// predModule.cpp: implementation of predictor module using Eigen
//
// Copyright (C)       2011 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4Eigen.

#include "predModule.h"

namespace lme4Eigen {
    using std::invalid_argument;

    typedef   Map<MatrixXd>    MMat;
    typedef   Map<VectorXd>    MVec;
    typedef   Map<VectorXi>    MiVec;

    merPredD::merPredD(S4 X, SEXP Lambdat, SEXP LamtUt, SEXP Lind,
		       SEXP RZX, SEXP Ut, SEXP Utr, SEXP V, SEXP VtV,
		       SEXP Vtr, SEXP Zt, SEXP beta0, SEXP delb, SEXP delu,
		       SEXP theta, SEXP u0)
	: d_X(       X),
          d_Zt(      as<MSpMatrixd>(Zt)),
	  d_theta(   as<MVec>(theta)),
	  d_Lind(    as<MiVec>(Lind)),
	  d_n(       d_X.rows()),
	  d_nnz(     -1),
	  d_p(       d_X.cols()),
	  d_q(       d_Zt.rows()),
	  d_Lambdat( as<MSpMatrixd>(Lambdat)),
	  d_RZX(     as<MMat>(RZX)),
	  d_V(       as<MMat>(V)),
	  d_VtV(     as<MMat>(VtV)),
	  d_Vtr(     as<MVec>(Vtr)),
	  d_Utr(     as<MVec>(Utr)),
	  d_beta0(   as<MVec>(beta0)),
	  d_delb(    as<MVec>(delb)),
	  d_delu(    as<MVec>(delu)),
	  d_u0(      as<MVec>(u0)),
	  d_Ut(      as<MSpMatrixd>(Ut)),
	  d_LamtUt(  as<MSpMatrixd>(LamtUt)),
	  d_RX(      d_p)
	  //	  d_LamtUtRestructure(false)
    {				// Check consistency of dimensions
	if (d_n != d_Zt.cols())
	    throw invalid_argument("Z dimension mismatch");
	if (d_Lind.size() != d_Lambdat.nonZeros())
	    throw invalid_argument("size of Lind does not match nonzeros in Lambda");
	// checking of the range of Lind is now done in R code for reference class
				// initialize beta0, u0, delb, delu and VtV
//	d_beta0.setZero(); d_u0.setZero(); d_delu.setZero(); d_delb.setZero();
//	d_V       = d_X;
	d_VtV.setZero().selfadjointView<Eigen::Upper>().rankUpdate(d_V.adjoint());
	d_RX.compute(d_VtV);	// ensure d_RX is initialized even in the 0-column X case

	setTheta(d_theta);	    // starting values into Lambda
        d_L.cholmod().final_ll = 1; // force an LL' decomposition
	updateLamtUt();
	d_L.analyzePattern(d_LamtUt); // perform symbolic analysis
        if (d_L.info() != Eigen::Success)
	    throw runtime_error("CholeskyDecomposition.analyzePattern failed");
    }

    void merPredD::updateLamtUt() {
	// This complicated code bypasses problems caused by Eigen's
	// sparse/sparse matrix multiplication pruning zeros.  The
	// Cholesky decomposition croaks if the structure of d_LamtUt changes.
	std::fill(d_LamtUt._valuePtr(), d_LamtUt._valuePtr() + d_LamtUt.nonZeros(), Scalar());
	for (Index j = 0; j < d_Ut.cols(); ++j) {
	    for(MSpMatrixd::InnerIterator rhsIt(d_Ut, j); rhsIt; ++rhsIt) {
		Scalar                       y(rhsIt.value());
		Index                        k(rhsIt.index());
//		Rcpp::Rcout << "d_Ut row k = " << k << ", v = " << y;
//		Rcpp::Rcout << ", Lambdat: col has " << d_Lambdat.innerNonZeros(k)
//			    << " nonzeros" << std::endl;
		MSpMatrixd::InnerIterator prdIt(d_LamtUt, j);
		for (MSpMatrixd::InnerIterator lhsIt(d_Lambdat, k); lhsIt; ++lhsIt) {
		    Index                    i = lhsIt.index();
//		    Rcpp::Rcout << "lhsIt: i = " << i << ", v = " << lhsIt.value() << std::endl;
//		    Rcpp::Rcout << "prdIt: i=" << prdIt.index() << std::endl;
		    while (prdIt && prdIt.index() != i) {
			++prdIt;
//			if (!prdIt) Rcpp::Rcout << "End of product column" << std::endl;
//			else Rcpp::Rcout << "i=" << prdIt.index() << std::endl;
		    }
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
	d_L.factorize_p(d_LamtUt, Eigen::ArrayXi(), 1.);
	d_ldL2 = ::M_chm_factor_ldetL2(d_L.factor());
    }

    void merPredD::setTheta(const VectorXd& theta) {
	if (theta.size() != d_theta.size())
	    throw invalid_argument("theta size mismatch");
				// update theta
	std::copy(theta.data(), theta.data() + theta.size(), d_theta.data());
				// update Lambdat
	int    *lipt = d_Lind.data();
	double *LamX = d_Lambdat._valuePtr(), *thpt = d_theta.data();
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

    void merPredD::updateXwts(const VectorXd& sqrtXwt) {
//	Rcpp::Rcout << "X[" << d_X.rows() << ", " << d_X.cols() << "]"
//		    << "V[" << d_V.rows() << ", " << d_V.cols() << "]"
//		    << ", sqrtXwt.size() = " << sqrtXwt.size() << std::endl;

	if (d_X.rows() != sqrtXwt.size())
	    throw invalid_argument("updateXwts: dimension mismatch");
	if (sqrtXwt.size() == d_V.rows()) {
	    d_V              = sqrtXwt.asDiagonal() * d_X;
	    for (int j = 0; j < d_n; ++j)
		for (MSpMatrixd::InnerIterator Uit(d_Ut, j), Zit(d_Zt, j);
		     Uit && Zit; ++Uit, ++Zit)
		    Uit.valueRef() = Zit.value() * sqrtXwt.data()[j];
	} else {
	    int n            = d_X.rows()/d_V.rows();
	    SpMatrixd      W(n, sqrtXwt.size());
	    const double *pt = sqrtXwt.data();
	    W.reserve(sqrtXwt.size());
	    for (Index j = 0; j < W.cols(); ++j, ++pt) {
		W.startVec(j);
		W.insertBack(j % n, j) = *pt;
	    }
	    W.finalize();
	    d_V              = W * d_X;
//FIXME: work out the corresponding code for Ut and Zt
//	    d_Ut             = d_Zt * W.adjoint();
	}
//	if (d_LamtUt.rows() != d_Ut.rows() || 
//	    d_LamtUt.cols() != d_Ut.cols()) d_LamtUtRestructure = true;
	d_VtV.setZero().selfadjointView<Eigen::Upper>().rankUpdate(d_V.adjoint());
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
			selfadjointView<Lower>().rankUpdate(RXi()));
    }

    VectorXd merPredD::RXdiag() const {
	return d_RX.matrixLLT().diagonal();
    }
}
