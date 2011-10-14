// external.cpp: externally .Call'able functions in lme4Eigen
//
// Copyright (C)       2011 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#include "predModule.h"
#include "respModule.h"

extern "C" {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    using Rcpp::CharacterVector;
    using Rcpp::Environment;
    using Rcpp::IntegerVector;
    using Rcpp::Language;
    using Rcpp::List;
    using Rcpp::NumericVector;
    using Rcpp::S4;
    using Rcpp::XPtr;
    using Rcpp::as;
    using Rcpp::wrap;

    using glm::glmFamily;

    using lme4Eigen::glmResp;
    using lme4Eigen::lmResp;
    using lme4Eigen::lmerResp;
    using lme4Eigen::merPredD;
    using lme4Eigen::nlsResp;

    using std::runtime_error;

    SEXP Eigen_SSE() {
	BEGIN_RCPP;
	return wrap(Eigen::SimdInstructionSetsInUse());
	END_RCPP;
    }

    // generalized linear model (and generalized linear mixed model) response

    SEXP glm_Create(SEXP fams, SEXP ys) {
	BEGIN_RCPP;
	glmResp *ans = new glmResp(List(fams), NumericVector(ys));
	return wrap(XPtr<glmResp>(ans, true));
	END_RCPP;
    }

    SEXP glm_setN(SEXP ptr_, SEXP n) {
	BEGIN_RCPP;
	XPtr<glmResp>(ptr_)->setN(as<VectorXd>(n));
	END_RCPP;
    }

    SEXP glm_devResid(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<glmResp>(ptr_)->devResid());
	END_RCPP;
    }

    SEXP glm_eta(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<glmResp>(ptr_)->eta());
	END_RCPP;
    }

    SEXP glm_family(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<glmResp>(ptr_)->family());
	END_RCPP;
    }

    SEXP glm_link(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<glmResp>(ptr_)->link());
	END_RCPP;
    }

    SEXP glm_muEta(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<glmResp>(ptr_)->muEta());
	END_RCPP;
    }

    SEXP glm_n(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<glmResp>(ptr_)->n());
	END_RCPP;
    }

    SEXP glm_resDev(SEXP ptr_) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<glmResp>(ptr_)->resDev());
	END_RCPP;
    }

    SEXP glm_sqrtWrkWt(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<glmResp>(ptr_)->sqrtWrkWt());
	END_RCPP;
    }

    SEXP glm_updateWts(SEXP ptr_) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<glmResp>(ptr_)->updateWts());
	END_RCPP;
    }

    SEXP glm_variance(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<glmResp>(ptr_)->variance());
	END_RCPP;
    }

    SEXP glm_wrkResids(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<glmResp>(ptr_)->wrkResids());
	END_RCPP;
    }

    SEXP glm_wrkResp(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<glmResp>(ptr_)->wrkResp());
	END_RCPP;
    }

    SEXP glm_Laplace(SEXP ptr_, SEXP ldL2, SEXP ldRX2, SEXP sqrL) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<glmResp>(ptr_)->Laplace(::Rf_asReal(ldL2),
									 ::Rf_asReal(ldRX2),
									 ::Rf_asReal(sqrL)));
	END_RCPP;
    }

    SEXP glm_updateMu(SEXP ptr_, SEXP gamma) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<glmResp>(ptr_)->updateMu(as<Map<VectorXd> >(gamma)));
	END_RCPP;
    }

    // glm family objects

    SEXP glmFamily_Create(SEXP fam_) {
	BEGIN_RCPP;
	glmFamily *ans = new glmFamily(List(fam_));
	return wrap(XPtr<glmFamily>(ans, true));
	END_RCPP;
    }

    SEXP glmFamily_link(SEXP ptr, SEXP mu) {
	BEGIN_RCPP;
	return wrap(XPtr<glmFamily>(ptr)->linkFun(as<Map<VectorXd> >(mu)));
	END_RCPP;
    }

    SEXP glmFamily_linkInv(SEXP ptr, SEXP eta) {
	BEGIN_RCPP;
	return wrap(XPtr<glmFamily>(ptr)->linkInv(as<Map<VectorXd> >(eta)));
	END_RCPP;
    }

    SEXP glmFamily_devResid(SEXP ptr, SEXP mu, SEXP weights, SEXP y) {
	BEGIN_RCPP;
	return wrap(XPtr<glmFamily>(ptr)->devResid(as<Map<VectorXd> >(mu),
						   as<Map<VectorXd> >(weights),
						   as<Map<VectorXd> >(y)));
	END_RCPP;
    }

    SEXP glmFamily_muEta(SEXP ptr, SEXP eta) {
	BEGIN_RCPP;
	return wrap(XPtr<glmFamily>(ptr)->muEta(as<Map<VectorXd> >(eta)));
	END_RCPP;
    }

    SEXP glmFamily_variance(SEXP ptr, SEXP mu) {
	BEGIN_RCPP;
	return wrap(XPtr<glmFamily>(ptr)->variance(as<VectorXd>(mu)));
	END_RCPP;
    }

    inline double pwrss(lmResp *rp, merPredD *pp, double fac) {
	return rp->wrss() + (fac ? pp->sqrL(fac) : pp->u0().squaredNorm());
    }

    void stepFac(glmResp *rp, merPredD *pp, int verb) {
	double pwrss0(pwrss(rp, pp, 0.));

	for (double fac = 1.; fac > 0.001; fac /= 2.) {
	    double pwrss1 = rp->updateMu(pp->linPred(fac)) + pp->sqrL(fac);
	    if (verb > 3)
		::Rprintf("pwrss0=%10g, diff=%10g, fac=%6.4f\n",
			  pwrss0, pwrss0 - pwrss1, fac);
	    if (pwrss1 < pwrss0) {
		pp->installPars(fac);
		return;
	    }
	}
	throw runtime_error("step factor reduced below 0.001 without reducing pwrss");
    }

    void pwrssUpdate(glmResp *rp, merPredD *pp, int verb, bool uOnly, double tol) {
				// FIXME: limit the number of iterations?
	do {
	    rp->updateMu(pp->linPred(0.));
	    rp->updateWts();
	    pp->updateXwts(rp->sqrtXwt());
	    pp->updateDecomp();
	    pp->updateRes(rp->wtres());
	    if ((uOnly ? pp->solveU() : pp->solve())/pwrss(rp, pp, 0.) < tol) break;
	    stepFac(rp, pp, verb);
	} while (true);
    }

    SEXP glmerPwrssUpdate(SEXP pp_, SEXP rp_, SEXP verb_, SEXP uOnly_, SEXP tol) {
	BEGIN_RCPP;
	pwrssUpdate(XPtr<glmResp>(rp_), XPtr<merPredD>(pp_),
		    ::Rf_asInteger(verb_), ::Rf_asLogical(uOnly_), ::Rf_asReal(tol));
	END_RCPP;
    }

    SEXP glmerWrkIter(SEXP pp_, SEXP rp_) {
	BEGIN_RCPP;

	XPtr<glmResp>    rp(rp_);
	XPtr<merPredD>   pp(pp_);
	const VectorXd   wt(rp->sqrtWrkWt());
	pp->updateXwts(wt);
	pp->updateDecomp();
	pp->updateRes(rp->wrkResp());
	pp->solve();
	rp->updateMu(pp->linPred(1.));
	pp->installPars(1.);

	return ::Rf_ScalarReal(rp->updateWts());

	END_RCPP;
    }

    SEXP glmerLaplace(SEXP pp_, SEXP rp_, SEXP theta_, SEXP u0_, SEXP beta0_,
		      SEXP verbose_, SEXP uOnly_, SEXP tol_) {
	BEGIN_RCPP;

	XPtr<glmResp>     rp(rp_);
	XPtr<merPredD>    pp(pp_);

	pp->setTheta(NumericVector(theta_));
	pp->setU0(as<VectorXd>(u0_));
	pp->setBeta0(as<VectorXd>(beta0_));
	pwrssUpdate(rp, pp, ::Rf_asInteger(verbose_), ::Rf_asLogical(uOnly_),
		    ::Rf_asReal(tol_));
	return ::Rf_ScalarReal(rp->Laplace(pp->ldL2(), pp->ldRX2(), pp->sqrL(1.)));

	END_RCPP;
    }

    SEXP isNullExtPtr(SEXP Ptr) {
	BEGIN_RCPP;
	return ::Rf_ScalarLogical(XPtr<lmResp>(Ptr) == (lmResp*)NULL);
	END_RCPP;
    }

    // linear model response (also the base class for other response classes)

    SEXP lm_Create(SEXP ys) {
	BEGIN_RCPP;
	lmResp *ans = new lmResp(NumericVector(ys));
	return wrap(XPtr<lmResp>(ans, true));
	END_RCPP;
    }

    SEXP lm_setOffset(SEXP ptr_, SEXP offset) {
	BEGIN_RCPP;
	XPtr<lmResp>(ptr_)->setOffset(as<VectorXd>(offset));
	return offset;
	END_RCPP;
    }

    SEXP lm_setWeights(SEXP ptr_, SEXP weights) {
	BEGIN_RCPP;
	XPtr<lmResp>(ptr_)->setWeights(as<VectorXd>(weights));
	return weights;
	END_RCPP;
    }

    SEXP lm_mu(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lmResp>(ptr_)->mu());
	END_RCPP;
    }

    SEXP lm_offset(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lmResp>(ptr_)->offset());
	END_RCPP;
    }

    SEXP lm_sqrtXwt(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lmResp>(ptr_)->sqrtXwt());
	END_RCPP;
    }

    SEXP lm_sqrtrwt(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lmResp>(ptr_)->sqrtrwt());
	END_RCPP;
    }

    SEXP lm_weights(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lmResp>(ptr_)->weights());
	END_RCPP;
    }

    SEXP lm_wrss(SEXP ptr_) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lmResp>(ptr_)->wrss());
	END_RCPP;
    }

    SEXP lm_wtres(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lmResp>(ptr_)->wtres());
	END_RCPP;
    }

    SEXP lm_y(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lmResp>(ptr_)->y());
	END_RCPP;
    }

    SEXP lm_updateMu(SEXP ptr_, SEXP gamma) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lmerResp>(ptr_)->updateMu(as<Map<VectorXd> >(gamma)));
	END_RCPP;
    }

    // linear mixed-effects model response

    SEXP lmer_Create(SEXP ys) {
	BEGIN_RCPP;
	lmerResp *ans = new lmerResp(NumericVector(ys));
	return wrap(XPtr<lmerResp>(ans, true));
	END_RCPP;
    }

    SEXP lmer_setREML(SEXP ptr_, SEXP REML) {
	BEGIN_RCPP;
	int reml = ::Rf_asInteger(REML);
	XPtr<lmerResp>(ptr_)->setReml(reml);
	return ::Rf_ScalarInteger(reml);
	END_RCPP;
    }

    SEXP lmer_REML(SEXP ptr_) {
	BEGIN_RCPP;
	return ::Rf_ScalarInteger(XPtr<lmerResp>(ptr_)->REML());
	END_RCPP;
    }

    SEXP lmer_Laplace(SEXP ptr_, SEXP ldL2, SEXP ldRX2, SEXP sqrL) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lmerResp>(ptr_)->Laplace(::Rf_asReal(ldL2),
									::Rf_asReal(ldRX2),
									::Rf_asReal(sqrL)));
	END_RCPP;
    }

    SEXP lmer_Deviance(SEXP pptr_, SEXP rptr_, SEXP theta_) {
	// Assume that ppt->updateWts(rpt->sqrtXwt()) has been called once
	BEGIN_RCPP;

	XPtr<lmerResp>   rpt(rptr_);
	XPtr<merPredD>   ppt(pptr_);
	ppt->setTheta(NumericVector(theta_));
	ppt->updateDecomp();
        rpt->updateMu(ppt->linPred(0.));
        ppt->updateRes(rpt->wtres());
	ppt->solve();
        rpt->updateMu(ppt->linPred(1.));
	return ::Rf_ScalarReal(rpt->Laplace(ppt->ldL2(), ppt->ldRX2(), ppt->sqrL(1.)));

	END_RCPP;
    }

    // dense predictor module for mixed-effects models

    SEXP merPredDCreate(SEXP Xs, SEXP Zts, SEXP Lambdats, SEXP Linds, SEXP thetas) {
	BEGIN_RCPP;
	S4 X(Xs), Zt(Zts), Lambdat(Lambdats);
	IntegerVector Lind(Linds);
	NumericVector theta(thetas);
	merPredD *ans = new merPredD(X, Zt, Lambdat, Lind, theta);
	return wrap(XPtr<merPredD>(ans, true));
	END_RCPP;
    }

				// setters
    SEXP merPredDsetBeta0(SEXP ptr, SEXP beta0) {
	BEGIN_RCPP;
	XPtr<merPredD>(ptr)->setBeta0(as<Map<VectorXd> >(beta0));
	END_RCPP;
    }

    SEXP merPredDsetTheta(SEXP ptr, SEXP theta) {
	BEGIN_RCPP;
	XPtr<merPredD>(ptr)->setTheta(NumericVector(theta));
	END_RCPP;
    }

    SEXP merPredDsetU0(SEXP ptr, SEXP u0) {
	BEGIN_RCPP;
	XPtr<merPredD>(ptr)->setU0(as<Map<VectorXd> >(u0));
	END_RCPP;
    }

				// getters
    SEXP merPredDCcNumer(SEXP ptr) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<merPredD>(ptr)->CcNumer());
	END_RCPP;
    }

    SEXP merPredDL(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->L());
	END_RCPP;
    }

    SEXP merPredDLambdat(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->Lambdat());
	END_RCPP;
    }

    SEXP merPredDLamtUt(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->LamtUt());
	END_RCPP;
    }

    SEXP merPredDPvec(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->Pvec());
	END_RCPP;
    }

    SEXP merPredDRX(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->RX());
	END_RCPP;
    }

    SEXP merPredDRXi(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->RXi());
	END_RCPP;
    }

    SEXP merPredDRXdiag(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->RXdiag());
	END_RCPP;
    }

    SEXP merPredDRZX(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->RZX());
	END_RCPP;
    }

    SEXP merPredDUt(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->Ut());
	END_RCPP;
    }

    SEXP merPredDUtr(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->Utr());
	END_RCPP;
    }

    SEXP merPredDV(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->V());
	END_RCPP;
    }

    SEXP merPredDVtV(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->VtV());
	END_RCPP;
    }

    SEXP merPredDVtr(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->Vtr());
	END_RCPP;
    }

    SEXP merPredDZt(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->Zt());
	END_RCPP;
    }

    SEXP merPredDbeta0(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->beta0());
	END_RCPP;
    }

    SEXP merPredDdelb(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->delb());
	END_RCPP;
    }

    SEXP merPredDdelu(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->delu());
	END_RCPP;
    }

    SEXP merPredDldL2(SEXP ptr) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<merPredD>(ptr)->ldL2());
	END_RCPP;
    }

    SEXP merPredDldRX2(SEXP ptr) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<merPredD>(ptr)->ldRX2());
	END_RCPP;
    }

    SEXP merPredDtheta(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->theta());
	END_RCPP;
    }

    SEXP merPredDu0(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->u0());
	END_RCPP;
    }

    SEXP merPredDunsc(SEXP ptr) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->unsc());
	END_RCPP;
    }

    // methods

    SEXP merPredDb(SEXP ptr, SEXP fac) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->b(::Rf_asReal(fac)));
	END_RCPP;
    }

    SEXP merPredDbeta(SEXP ptr, SEXP fac) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->beta(::Rf_asReal(fac)));
	END_RCPP;
    }

    SEXP merPredDinstallPars(SEXP ptr, SEXP fac) {
	BEGIN_RCPP;
	XPtr<merPredD>(ptr)->installPars(::Rf_asReal(fac));
	END_RCPP;
    }

    SEXP merPredDlinPred(SEXP ptr, SEXP fac) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->linPred(::Rf_asReal(fac)));
	END_RCPP;
    }

    SEXP merPredDsolve(SEXP ptr) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<merPredD>(ptr)->solve());
	END_RCPP;
    }

    SEXP merPredDsolveU(SEXP ptr) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<merPredD>(ptr)->solveU());
	END_RCPP;
    }

    SEXP merPredDsqrL(SEXP ptr, SEXP fac) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<merPredD>(ptr)->sqrL(::Rf_asReal(fac)));
	END_RCPP;
    }

    SEXP merPredDu(SEXP ptr, SEXP fac) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->u(::Rf_asReal(fac)));
	END_RCPP;
    }

				// methods
    SEXP merPredDupdateDecomp(SEXP ptr) {
	BEGIN_RCPP;
	XPtr<merPredD>(ptr)->updateDecomp();
	END_RCPP;
    }

    SEXP merPredDupdateRes(SEXP ptr, SEXP wtres) {
	BEGIN_RCPP;
	XPtr<merPredD>(ptr)->updateRes(as<Map<VectorXd> >(wtres));
	END_RCPP;
    }

    SEXP merPredDupdateXwts(SEXP ptr, SEXP wts) {
	BEGIN_RCPP;
	XPtr<merPredD>(ptr)->updateXwts(as<Map<MatrixXd> >(wts));
	END_RCPP;
    }

    SEXP nls_Create(SEXP ys, SEXP mods, SEXP envs, SEXP pnms, SEXP Ns) {
	BEGIN_RCPP;
	nlsResp *ans =
	    new nlsResp(NumericVector(ys), Language(mods),
			Environment(envs), CharacterVector(pnms),
			::Rf_asInteger(Ns));
	return wrap(XPtr<nlsResp>(ans, true));
	END_RCPP;
    }

    SEXP nls_Laplace(SEXP ptr_, SEXP ldL2, SEXP ldRX2, SEXP sqrL) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<nlsResp>(ptr_)->
			       Laplace(::Rf_asReal(ldL2),
				       ::Rf_asReal(ldRX2),
				       ::Rf_asReal(sqrL)));
	END_RCPP;
    }

    SEXP nls_updateMu(SEXP ptr_, SEXP gamma) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<nlsResp>(ptr_)->updateMu(as<Map<VectorXd> >(gamma)));
	END_RCPP;
    }


}

#include <R_ext/Rdynload.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {

    CALLDEF(Eigen_SSE, 0),

    CALLDEF(glm_Create, 2),	  // generate external pointer

    CALLDEF(glm_setN, 2),	  // setters

    CALLDEF(glm_devResid, 1),	  // getters
    CALLDEF(glm_eta, 1),
    CALLDEF(glm_family, 1),
    CALLDEF(glm_link, 1),
    CALLDEF(glm_muEta, 1),
    CALLDEF(glm_n, 1),
    CALLDEF(glm_resDev, 1),
    CALLDEF(glm_sqrtWrkWt, 1),
    CALLDEF(glm_variance, 1),
    CALLDEF(glm_wrkResids, 1),
    CALLDEF(glm_wrkResp, 1),

    CALLDEF(glm_Laplace, 4),	  // methods
    CALLDEF(glm_updateMu, 2),
    CALLDEF(glm_updateWts, 1),

    CALLDEF(glmFamily_Create, 1), // generate external pointer

    CALLDEF(glmFamily_link, 1),   // methods
    CALLDEF(glmFamily_linkInv, 1),
    CALLDEF(glmFamily_devResid, 3),
    CALLDEF(glmFamily_muEta, 1),
    CALLDEF(glmFamily_variance, 1),

    CALLDEF(glmerPwrssUpdate, 5),
    CALLDEF(glmerWrkIter, 2),
    CALLDEF(glmerLaplace, 8),

    CALLDEF(isNullExtPtr, 1),

    CALLDEF(lm_Create, 1),	  // generate external pointer

    CALLDEF(lm_setOffset, 2),	  // setters
    CALLDEF(lm_setWeights, 2),

    CALLDEF(lm_mu, 1),		  // getters
    CALLDEF(lm_offset, 1),
    CALLDEF(lm_sqrtXwt, 1),
    CALLDEF(lm_sqrtrwt, 1),
    CALLDEF(lm_weights, 1),
    CALLDEF(lm_wrss, 1),
    CALLDEF(lm_wtres, 1),
    CALLDEF(lm_y, 1),

    CALLDEF(lm_updateMu, 2),	  // method

    CALLDEF(lmer_Create, 1),	  // generate external pointer

    CALLDEF(lmer_setREML, 2),     // setter

    CALLDEF(lmer_REML, 1),	  // getter

    CALLDEF(lmer_Deviance, 3),    // method
    CALLDEF(lmer_Laplace, 4),

    CALLDEF(merPredDCreate, 5),	  // generate external pointer

    CALLDEF(merPredDsetTheta, 2), // setters
    CALLDEF(merPredDsetBeta0, 2),
    CALLDEF(merPredDsetU0, 2),

    CALLDEF(merPredDCcNumer, 1),  // getters
    CALLDEF(merPredDL, 1),
    CALLDEF(merPredDLambdat, 1),
    CALLDEF(merPredDLamtUt, 1),
    CALLDEF(merPredDPvec, 1),
    CALLDEF(merPredDRX, 1),
    CALLDEF(merPredDRXdiag, 1),
    CALLDEF(merPredDRZX, 1),
    CALLDEF(merPredDUt, 1),
    CALLDEF(merPredDUtr, 1),
    CALLDEF(merPredDV, 1),
    CALLDEF(merPredDVtV, 1),
    CALLDEF(merPredDVtr, 1),
    CALLDEF(merPredDZt, 1),
    CALLDEF(merPredDbeta0, 1),
    CALLDEF(merPredDdelb, 1),
    CALLDEF(merPredDdelu, 1),
    CALLDEF(merPredDldL2, 1),
    CALLDEF(merPredDldRX2, 1),
    CALLDEF(merPredDtheta, 1),
    CALLDEF(merPredDu0, 1),
    CALLDEF(merPredDunsc, 1),

    CALLDEF(merPredDb, 2),	// methods
    CALLDEF(merPredDbeta, 2),
    CALLDEF(merPredDlinPred, 2),
    CALLDEF(merPredDinstallPars, 2),
    CALLDEF(merPredDsolve, 1),
    CALLDEF(merPredDsolveU, 1),
    CALLDEF(merPredDsqrL, 2),
    CALLDEF(merPredDu, 2),
    CALLDEF(merPredDupdateDecomp, 1),
    CALLDEF(merPredDupdateRes, 2),
    CALLDEF(merPredDupdateXwts, 2),

    CALLDEF(nls_Create, 5),	  // generate external pointer

    CALLDEF(nls_Laplace, 4),	  // methods
    CALLDEF(nls_updateMu, 2),
    {NULL, NULL, 0}
};

/** Initializer for lme4Eigen, called upon loading the package.
 *
 *  Register routines that can be called directly from R.
 *  Initialize CHOLMOD and require the LL' form of the factorization.
 *  Install the symbols to be used by functions in the package.
 */
extern "C"
void R_init_lme4Eigen(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, (Rboolean)FALSE);
}

