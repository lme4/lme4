// lmer.cpp: linear mixed-effects models using Eigen
//
// Copyright (C)       2011 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#include "respModule.h"
#include "predModule.h"
#include "lmer.h"

extern "C" {
    using Rcpp::XPtr;
    using Rcpp::NumericVector;
    using std::cout;
    using std::endl;

    SEXP lmerDeviance(SEXP pptr_, SEXP rptr_, SEXP theta_) {
	// Assume that ppt->updateWts(rpt->sqrtXwt()) has been called once
	BEGIN_RCPP;
	Rprintf("In lmerDeviance\n"); 

	XPtr<lme4Eigen::lmerResp>   rpt(rptr_);
	XPtr<lme4Eigen::merPredD>   ppt(pptr_);
	ppt->setTheta(NumericVector(theta_));
	ppt->updateDecomp();
        rpt->updateMu(ppt->linPred(0.));
        ppt->updateRes(rpt->wtres());
	ppt->solve();
	Rprintf("Before extracting sqrL(1.)\n");
	cout << "u0: " << (ppt->u0()).adjoint() << endl;
	cout << "delu: " << (ppt->delu()).adjoint() << endl;
	cout << "u(1,): " << (ppt->u(1.)).adjoint() << endl;
	double sqrL = ppt->sqrL(1.);
	Rprintf("sqrL = %g\n", sqrL);
	double wrss = rpt->updateMu(ppt->linPred(1.));
	cout << "wrss =" << wrss << endl; 
	return ::Rf_ScalarReal(rpt->Laplace(ppt->ldL2(), ppt->ldRX2(), sqrL));

	END_RCPP;
    }
}

