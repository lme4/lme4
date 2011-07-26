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
    merPredD::merPredD(S4 X, S4 Z, S4 Lambda, IntegerVector Lind,
		       NumericVector theta)             throw (invalid_argument,
							       runtime_error)
	: d_X(X),
          d_Z(Z),
	  d_theta(clone(theta)),
	  d_Lind(Lind),
	  d_n(d_X.rows()),
	  d_p(d_X.cols()),
	  d_q(d_Z.cols()),
	  d_Lambda(Lambda),
	  d_RZX(d_q, d_p),
	  d_V(d_X),
	  d_Vtr(d_p),
	  d_Utr(d_q),
	  d_delb(d_p),
	  d_delu(d_q),
	  d_beta0(d_p),
	  d_u0(d_q),
	  d_U(d_Z),
	  d_I(d_q, d_q),
	  d_RX(d_p)
    {				// Check consistency of dimensions
	if (d_n != d_Z.rows())
	    throw invalid_argument("Z dimension mismatch");
	// cout << "dimensions: n = " << d_n
	//      << ", p = " << d_p
	//      << ", q = " << d_q
	//      << ", length(theta) = " << d_theta.size()
	//      << endl;
	    
	if (d_Lind.size() != d_Lambda.nonZeros())
	    throw invalid_argument("size of Lind does not match nonzeros in Lambda");
	// check the range of Lind (now done in R code for reference class)
	// int thsz = d_theta.size(), *Li = Lind.begin();
	// ArrayXi  used(thsz);
	// used.setZero();
	// for (int i = 0; i < Lind.size(); ++i) {
	//     int ind = Li[i];
	//     if (ind < 1 || ind > thsz)
	// 	throw invalid_argument("Lind size mismatch or range mismatch for theta");
	//     used[ind - 1] = 1; // 0-based index
	// }
	// if (!used.all()) throw invalid_argument("Some indices of theta do not occur in Lind");
				// initialize beta0 and u0
	d_beta0.setZero(); d_u0.setZero();
				// set d_I to be the d_q by d_q sparse Identity matrix
	d_I.reserve(d_q);
	for (Index j = 0; j < d_q; ++j) {
	    d_I.startVec(j);
	    d_I.insertBack(j, j) = 1.;
	}
	d_I.finalize();
				// finish with steps that are followed for each theta change
	setTheta(d_theta);
    }

    void merPredD::updateL(const SpMatrixXd& ULam) {
	// Create Lambda'U'U Lambda + I
	SpMatrixXd   UtU(d_I);
	UtU.selfadjointView<Eigen::Lower>().rankUpdate(ULam.adjoint());

	if (d_L.isInitialized()) d_L.factorize(UtU); else d_L.compute(UtU);
	if (d_L.info() != Eigen::Success) throw runtime_error("factorization failure");
	d_ldL2            = d_L.vectorD().array().log().sum();
    }

    void merPredD::setTheta(const NumericVector& theta) throw (invalid_argument,
							       runtime_error) {
	if (theta.size() != d_theta.size())
	    throw invalid_argument("theta size mismatch");
				// update theta
	copy(theta.begin(), theta.end(), d_theta.begin());
				// update Lambda and L
	int    *lipt = d_Lind.begin();
	double *LamX = d_Lambda._valuePtr(), *thpt = d_theta.begin();
	for (int i = 0; i < d_Lind.size(); ++i) {
	    LamX[i] = thpt[lipt[i] - 1];
	}
	updateL(d_U * d_Lambda);
    }

    void merPredD::solve() {
        DiagType  sqrtDi(d_L.vectorD().array().sqrt().inverse().matrix());
	d_delu           = d_L.permutationPinv() * d_Utr;
	d_L.matrixLDL().triangularView<Eigen::UnitLower>().solveInPlace(d_delu);
	d_delu           = sqrtDi * d_delu;
				// d_delu now contains cu
	d_delb           = d_RX.solve(d_Vtr - d_RZX.adjoint() * d_delu);
	d_delu           = sqrtDi * (d_delu - d_RZX * d_delb);
	d_L.matrixLDL().adjoint().triangularView<Eigen::UnitUpper>().solveInPlace(d_delu);
	d_delu           = d_L.permutationP() * d_delu;
    }
static bool nonFinite(const double& x) {
    return !::R_finite(x);
}

static bool allFinite(const double* bb, const double* ee) {
    if (find_if(bb, ee, &nonFinite) < ee) return false;
    return true;
}

static bool chkFinite(const VectorXd& x) {
    return allFinite(x.data(), x.data() + x.size());
}

static bool chkFinite(const MatrixXd& x) {
    return allFinite(x.data(), x.data() + x.size());
}

static bool chkFinite(const SpMatrixXd& x) {
    return allFinite(x._valuePtr(), x._valuePtr() + x.nonZeros());
}

    void merPredD::updateRes(const VectorType& wtres,
			     const VectorType& sqrtXwt) throw (invalid_argument) {
	if (!chkFinite(wtres)) throw invalid_argument("updateRes: nonfinite residual");
	if (!chkFinite(sqrtXwt)) throw invalid_argument("updateRes: nonfinite weights");
	if (d_V.rows() != wtres.size() || d_X.rows() != sqrtXwt.size())
	    throw invalid_argument("updateRes: dimension mismatch");
	if (d_V.rows() == d_X.rows()) {  //FIXME: Generalize this for nlmer
	    DiagType  W(sqrtXwt.asDiagonal());
	    d_V         = W * d_X;
	    d_U         = W * d_Z;
	} else throw invalid_argument("updateRes: no provision for nlmer yet");
	if (!chkFinite(d_V)) ::Rf_error("nonfinite d_V");
	if (!chkFinite(d_U)) ::Rf_error("nonfinite d_U");
	d_Vtr           = d_V.adjoint() * wtres;
//cout << "Vtr: " << d_Vtr.adjoint() << endl;
	SpMatrixXd ULam(d_U * d_Lambda);
	d_Utr           = ULam.adjoint() * wtres;
//cout << "Utr: " << d_Utr.adjoint() << endl;
				// update L, RZX and RX
	updateL(ULam);
	d_RZX           = d_L.permutationPinv() * (ULam.adjoint() * d_V);
	d_L.matrixL().solveInPlace(d_RZX);
	d_RZX           = d_L.vectorD().array().sqrt().inverse().matrix().asDiagonal() * d_RZX;

	MatrixXd     VtV(d_p, d_p);
	d_RX.compute(VtV.setZero().selfadjointView<Eigen::Upper>().
		     rankUpdate(d_V.adjoint()).
		     rankUpdate(d_RZX.adjoint(), -1));
	d_ldRX2         = 2. * d_RX.matrixLLT().diagonal().array().log().sum();
    }

}

extern "C" {
    SEXP merPredDCreate(SEXP Xs, SEXP Zs, SEXP Lambdas, SEXP Linds, SEXP thetas) {
	BEGIN_RCPP;
	S4 X(Xs), Z(Zs), Lambda(Lambdas);
	IntegerVector Lind(Linds);
	NumericVector theta(thetas);
	lme4Eigen::merPredD *ans = new lme4Eigen::merPredD(X, Z, Lambda, Lind, theta);
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
    
    SEXP merPredDLambda(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->Lambda());
	END_RCPP;
    }
    
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
    
    SEXP merPredDRZX(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->RZX());
	END_RCPP;
    }
    
    SEXP merPredDZ(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::merPredD>(ptr)->Z());
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
