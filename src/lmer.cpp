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
	BEGIN_RCPP;
	XPtr<lme4Eigen::lmerResp>   rpt(rptr_);
	XPtr<lme4Eigen::merPredD>   ppt(pptr_);
//cout << "In lmerDeviance, theta:" << as<VectorXd>(theta_).adjoint() << endl;
	ppt->setTheta(NumericVector(theta_));
//cout << "past setTheta, ldL2 = " << ppt->ldL2() << endl;
        rpt->updateMu(ppt->linPred(0.));
//cout << "past updateMu(0.), wrss = " << rpt->wrss() <<  endl;
        ppt->updateRes(rpt->wtres(), rpt->sqrtXwt());
//cout << "past updateRes, ldRX2 = " << ppt->ldRX2() <<  endl;
	ppt->solve();
//cout << "past solve, sqrL = " << ppt->sqrL(1.) << ", delb:" << (ppt->delb()).adjoint() <<  endl;
	rpt->updateMu(ppt->linPred(1.));
//cout << "past updateMu(1.), wrss = " << rpt->wrss() <<  endl;
	return ::Rf_ScalarReal(rpt->Laplace(ppt->ldL2(), ppt->ldRX2(), ppt->sqrL(1.)));
	END_RCPP;
    }
}

