// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// external.h: externally .Call'able functions in lme4Eigen
//
// Copyright (C)       2011 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#ifndef LME4_LMER_H
#define LME4_LMER_H
#define R_NO_REMAP		// to suppress declaring length as a macro
#include <Rinternals.h>		// for SEXP declaration

extern "C" {
    // generalized linear model (and generalized linear mixed model) response

    SEXP glm_Create(SEXP,SEXP);	// constructor

    SEXP glm_setN(SEXP,SEXP);	// setter

    SEXP glm_devResid(SEXP);	// getters
    SEXP glm_eta(SEXP);
    SEXP glm_family(SEXP);
    SEXP glm_link(SEXP);
    SEXP glm_muEta(SEXP);
    SEXP glm_n(SEXP);		
    SEXP glm_resDev(SEXP);
    SEXP glm_sqrtWrkWt(SEXP);
    SEXP glm_updateWts(SEXP);
    SEXP glm_variance(SEXP);
    SEXP glm_wrkResids(SEXP);
    SEXP glm_wrkResp(SEXP);

    SEXP glm_Laplace(SEXP,SEXP,SEXP,SEXP); // methods
    SEXP glm_updateMu(SEXP,SEXP);

    // linear model response (also the base class for other response classes)

    SEXP lm_Create(SEXP);	   // constructor

    SEXP lm_setOffset(SEXP,SEXP);  // setters
    SEXP lm_setWeights(SEXP,SEXP);

    SEXP lm_mu(SEXP);		   // getters
    SEXP lm_offset(SEXP);
    SEXP lm_sqrtXwt(SEXP);
    SEXP lm_sqrtrwt(SEXP);
    SEXP lm_weights(SEXP);
    SEXP lm_wrss(SEXP);
    SEXP lm_wtres(SEXP);
    SEXP lm_y(SEXP);

    SEXP lm_updateMu(SEXP,SEXP);   // methods

    // linear mixed model response

    SEXP lmer_Create(SEXP);	   // constructor

    SEXP lmer_setREML(SEXP, SEXP); // setter

    SEXP lmer_REML(SEXP);	   // getter

    SEXP lmer_Laplace(SEXP,SEXP,SEXP,SEXP); // method

    SEXP lmerDeviance(SEXP, SEXP, SEXP);

    // dense predictor module
    SEXP merPredDCreate(SEXP, SEXP, SEXP, SEXP, SEXP); // constructor (returns external pointer)
    
    SEXP merPredDsetTheta(SEXP, SEXP); // setters
    SEXP merPredDsetBeta0(SEXP, SEXP);
    SEXP merPredDsetU0(SEXP, SEXP);

    SEXP merPredDCcNumer(SEXP);	       // getters
    SEXP merPredDL(SEXP);	
    SEXP merPredDLambdat(SEXP);
    SEXP merPredDLamtUt(SEXP);
    SEXP merPredDPvec(SEXP);
    SEXP merPredDRX(SEXP);
    SEXP merPredDRXdiag(SEXP);
    SEXP merPredDRXi(SEXP);
    SEXP merPredDRZX(SEXP);
    SEXP merPredDUt(SEXP);
    SEXP merPredDUtr(SEXP);
    SEXP merPredDV(SEXP);
    SEXP merPredDVtV(SEXP);
    SEXP merPredDVtr(SEXP);
    SEXP merPredDZt(SEXP);
    SEXP merPredDbeta0(SEXP);
    SEXP merPredDdelb(SEXP);
    SEXP merPredDdelu(SEXP);
    SEXP merPredDldL2(SEXP);
    SEXP merPredDldRX2(SEXP);
    SEXP merPredDtheta(SEXP);
    SEXP merPredDu0(SEXP);
    SEXP merPredDunsc(SEXP);

    SEXP merPredDb(SEXP, SEXP);	       // methods
    SEXP merPredDbeta(SEXP, SEXP);
    SEXP merPredDlinPred(SEXP, SEXP);
    SEXP merPredDinstallPars(SEXP, SEXP);
    SEXP merPredDsolve(SEXP);
    SEXP merPredDsolveU(SEXP);
    SEXP merPredDsqrL(SEXP,SEXP);
    SEXP merPredDu(SEXP, SEXP);
    SEXP merPredDupdateDecomp(SEXP);
    SEXP merPredDupdateRes(SEXP,SEXP);
    SEXP merPredDupdateXwts(SEXP,SEXP);

    // nonlinear model (and nonlinear mixed model) class

    SEXP nlm_Create(SEXP,SEXP,SEXP,SEXP); // constructor

    SEXP nlm_Laplace(SEXP,SEXP,SEXP,SEXP); // methods
    SEXP nlm_updateMu(SEXP,SEXP);
}

#endif
