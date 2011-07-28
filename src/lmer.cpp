// lmer.cpp: linear mixed-effects models using Eigen
//
// Copyright (C)       2011 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#include "respModule.h"
#include "predModule.h"
#include "lmer.h"

using namespace Rcpp;
using namespace std;
extern "C" {
    SEXP lmerDeviance(SEXP pptr_, SEXP rptr_, SEXP theta_) {
	// Assume that ppt->updateWts(rpt->sqrtXwt()) has been called once
	BEGIN_RCPP;
	XPtr<lme4Eigen::lmerResp>   rpt(rptr_);
	XPtr<lme4Eigen::merPredD>   ppt(pptr_);
	ppt->setTheta(NumericVector(theta_));
	ppt->updateDecomp();
        rpt->updateMu(ppt->linPred(0.));
        ppt->updateRes(rpt->wtres());
	ppt->solve();
	rpt->updateMu(ppt->linPred(1.));
	return ::Rf_ScalarReal(rpt->Laplace(ppt->ldL2(), ppt->ldRX2(), ppt->sqrL(1.)));
	END_RCPP;
    }
}

