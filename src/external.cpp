// external.cpp: externally .Call'able functions in lme4Eigen
//
// Copyright (C)       2011-2012 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#include "predModule.h"
#include "respModule.h"
#include "optimizer.h"

extern "C" {
    using     Eigen::ArrayXd;
    typedef   Eigen::Map<Eigen::MatrixXd>     MMat;
    typedef   Eigen::Map<Eigen::VectorXd>     MVec;
    typedef   Eigen::Map<Eigen::VectorXi>    MiVec;

    using      Rcpp::CharacterVector;
    using      Rcpp::Environment;
    using      Rcpp::IntegerVector;
    using      Rcpp::Language;
    using      Rcpp::List;
    using      Rcpp::Named;
    using      Rcpp::NumericVector;
    using      Rcpp::XPtr;
    using      Rcpp::as;
    using      Rcpp::wrap;

    using       glm::glmFamily;

    using lme4Eigen::glmResp;
    using lme4Eigen::lmResp;
    using lme4Eigen::lmerResp;
    using lme4Eigen::merPredD;
    using lme4Eigen::nlsResp;

    using optimizer::Golden;
    using optimizer::Nelder_Mead;
    using optimizer::nm_status;

    using      std::runtime_error;

    SEXP Eigen_SSE() {
	BEGIN_RCPP;
	return wrap(Eigen::SimdInstructionSetsInUse());
	END_RCPP;
    }

    // generalized linear model (and generalized linear mixed model) response

    SEXP glm_Create(SEXP fam, SEXP y, SEXP weights, SEXP offset, SEXP mu,
		    SEXP sqrtXwt, SEXP sqrtrwt, SEXP wtres, SEXP eta, SEXP n) {
	BEGIN_RCPP;
	glmResp *ans = new glmResp(List(fam), y, weights, offset, mu,
				   sqrtXwt, sqrtrwt, wtres, eta, n);
	return wrap(XPtr<glmResp>(ans, true));
	END_RCPP;
    }

    SEXP glm_setN(SEXP ptr_, SEXP n) {
	BEGIN_RCPP;
	XPtr<glmResp>(ptr_)->setN(as<MVec>(n));
	END_RCPP;
    }

    SEXP glm_devResid(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<glmResp>(ptr_)->devResid());
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
	return ::Rf_ScalarReal(XPtr<glmResp>(ptr_)->updateMu(as<MVec>(gamma)));
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
	return wrap(XPtr<glmFamily>(ptr)->linkFun(as<MVec>(mu)));
	END_RCPP;
    }

    SEXP glmFamily_linkInv(SEXP ptr, SEXP eta) {
	BEGIN_RCPP;
	return wrap(XPtr<glmFamily>(ptr)->linkInv(as<MVec>(eta)));
	END_RCPP;
    }

    SEXP glmFamily_devResid(SEXP ptr, SEXP mu, SEXP weights, SEXP y) {
	BEGIN_RCPP;
	return wrap(XPtr<glmFamily>(ptr)->devResid(as<MVec>(mu),
						   as<MVec>(weights),
						   as<MVec>(y)));
	END_RCPP;
    }

    SEXP glmFamily_muEta(SEXP ptr, SEXP eta) {
	BEGIN_RCPP;
	return wrap(XPtr<glmFamily>(ptr)->muEta(as<MVec>(eta)));
	END_RCPP;
    }

    SEXP glmFamily_variance(SEXP ptr, SEXP mu) {
	BEGIN_RCPP;
	return wrap(XPtr<glmFamily>(ptr)->variance(as<MVec>(mu)));
	END_RCPP;
    }

    static inline double pwrss(lmResp *rp, merPredD *pp, double fac) {
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

#define MAXITER 30
    static void pwrssUpdate(glmResp *rp, merPredD *pp, int verb, bool uOnly, double tol) {
	bool cvgd(false);
	for (int it=0; it < MAXITER; ++it) {
	    rp->updateMu(pp->linPred(0.));
	    rp->updateWts();
	    pp->updateXwts(rp->sqrtXwt());
	    pp->updateDecomp();
	    pp->updateRes(rp->wtres());
	    if ((uOnly ? pp->solveU() : pp->solve())/pwrss(rp, pp, 0.) < tol) {
		cvgd = true;
		break;
	    }
	    stepFac(rp, pp, verb);
	}
	if (!cvgd) throw runtime_error("pwrss failed to converge in 30 iterations");
    }

    static ArrayXd devcCol(const MiVec& fac, const ArrayXd& u, const ArrayXd& devRes) {
	ArrayXd  ans(u * u);
	for (int i = 0; i < devRes.size(); ++i) ans[fac[i] - 1] += devRes[i];
	return ans;
    }

    static double sqrt2pi = std::sqrt(2. * PI);

    SEXP glmerAGQ(SEXP pp_, SEXP rp_, SEXP fac_, SEXP GQmat_, SEXP theta_,
		  SEXP u0_, SEXP beta0_, SEXP tol_) {
	BEGIN_RCPP;
	
	XPtr<glmResp>     rp(rp_);
	XPtr<merPredD>    pp(pp_);
	const MiVec      fac(as<MiVec>(fac_));
	if (fac.size() != rp->mu().size())
	    throw std::invalid_argument("size of fac must match dimension of response vector");

	pp->setTheta(as<MVec>(theta_));
	pp->setU0(as<MVec>(u0_));
	pp->setBeta0(as<MVec>(beta0_));
	pwrssUpdate(rp, pp, 0, true, ::Rf_asReal(tol_)); // should be a no-op
	const ArrayXd     u0(pp->u0());
	const ArrayXd  devc0(devcCol(fac, u0, rp->devResid()));

	const unsigned int q(u0.size());
	if (pp->L().factor()->nzmax !=  q)
	    throw std::invalid_argument("AGQ only defined for a single scalar random-effects term");
	const ArrayXd     sd(Eigen::Map<ArrayXd>((double*)pp->L().factor()->x, q).inverse());

	const MMat     GQmat(as<MMat>(GQmat_));
	ArrayXd         mult(q);

	mult.setZero();
	for (int i = 0; i < GQmat.rows(); ++i) {
	    double            zknot(GQmat(i, 0));
	    if (zknot == 0)    mult += ArrayXd::Constant(q, GQmat(i, 1));
	    else {
		const ArrayXd   u0i(u0 + zknot * sd);
		pp->setU0(u0i);
		rp->updateMu(pp->linPred(0.));
		mult += (-0.5 * (devcCol(fac, u0i, rp->devResid()) - devc0) - GQmat(i, 2)).exp() *
		    GQmat(i, 1)/sqrt2pi;
	    }
	}
	pp->setU0(u0);		// restore settings from pwrssUpdate;
	rp->updateMu(pp->linPred(0.));
	return ::Rf_ScalarReal(devc0.sum() + pp->ldL2() - 2 * std::log(mult.prod()));
	END_RCPP;
    }

    SEXP glmerPwrssUpdate(SEXP pp_, SEXP rp_, SEXP verb_, SEXP uOnly_, SEXP tol) {
	BEGIN_RCPP;
	XPtr<glmResp> rp(rp_);
	XPtr<merPredD> pp(pp_);
	pwrssUpdate(rp, pp, ::Rf_asInteger(verb_), ::Rf_asLogical(uOnly_), ::Rf_asReal(tol));
	END_RCPP;
    }

    SEXP glmerWrkIter(SEXP pp_, SEXP rp_) {
	BEGIN_RCPP;

	XPtr<glmResp>    rp(rp_);
	XPtr<merPredD>   pp(pp_);
	const Eigen::VectorXd   wt(rp->sqrtWrkWt());
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

	pp->setTheta(as<MVec>(theta_));
	pp->setU0(as<MVec>(u0_));
	pp->setBeta0(as<MVec>(beta0_));
	pwrssUpdate(rp, pp, ::Rf_asInteger(verbose_), ::Rf_asLogical(uOnly_),
		    ::Rf_asReal(tol_));
	return ::Rf_ScalarReal(rp->Laplace(pp->ldL2(), pp->ldRX2(), pp->sqrL(1.)));

	END_RCPP;
    }

    static inline double prss(lmResp *rp, merPredD *pp, double fac) {
	return rp->wrss() + (fac ? pp->sqrL(fac) : pp->u0().squaredNorm());
    }

    void nstepFac(nlsResp *rp, merPredD *pp, int verb) {
	double prss0(prss(rp, pp, 0.));

	for (double fac = 1.; fac > 0.001; fac /= 2.) {
	    double prss1 = rp->updateMu(pp->linPred(fac)) + pp->sqrL(fac);
	    if (verb > 3)
		::Rprintf("pwrss0=%10g, diff=%10g, fac=%6.4f\n",
			  prss0, prss0 - prss1, fac);
	    if (prss1 < prss0) {
		pp->installPars(fac);
		return;
	    }
	}
	throw runtime_error("step factor reduced below 0.001 without reducing pwrss");
    }

#define MAXITER 30
    static void prssUpdate(nlsResp *rp, merPredD *pp, int verb, bool uOnly, double tol) {
	bool cvgd(false);
	for (int it=0; it < MAXITER; ++it) {
	    rp->updateMu(pp->linPred(0.));
	    pp->updateXwts(rp->sqrtXwt());
	    pp->updateDecomp();
	    pp->updateRes(rp->wtres());
	    if ((uOnly ? pp->solveU() : pp->solve())/pwrss(rp, pp, 0.) < tol) {
		cvgd = true;
		break;
	    }
	    nstepFac(rp, pp, verb);
	}
	if (!cvgd) throw runtime_error("pwrss failed to converge in 30 iterations");
    }

    SEXP nlmerLaplace(SEXP pp_, SEXP rp_, SEXP theta_, SEXP u0_, SEXP beta0_,
		      SEXP verbose_, SEXP uOnly_, SEXP tol_) {
	BEGIN_RCPP;

	XPtr<nlsResp>     rp(rp_);
	XPtr<merPredD>    pp(pp_);

	pp->setTheta(as<MVec>(theta_));
	pp->setU0(as<MVec>(u0_));
	pp->setBeta0(as<MVec>(beta0_));
	prssUpdate(rp, pp, ::Rf_asInteger(verbose_), ::Rf_asLogical(uOnly_),
		    ::Rf_asReal(tol_));
	return ::Rf_ScalarReal(rp->Laplace(pp->ldL2(), pp->ldRX2(), pp->sqrL(1.)));

	END_RCPP;
    }

    SEXP golden_Create(SEXP lower_, SEXP upper_) {
	BEGIN_RCPP;
	Golden *ans = new Golden(::Rf_asReal(lower_), ::Rf_asReal(upper_));
	return wrap(XPtr<Golden>(ans, true));
	END_RCPP;
    }

    SEXP golden_newf(SEXP ptr_, SEXP f_) {
	BEGIN_RCPP;
	XPtr<Golden>(ptr_)->newf(::Rf_asReal(f_));
	END_RCPP;
    }
	
    SEXP golden_xeval(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<Golden>(ptr_)->xeval());
	END_RCPP;
    }

    SEXP golden_value(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<Golden>(ptr_)->value());
	END_RCPP;
    }

    SEXP golden_xpos(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<Golden>(ptr_)->xpos());
	END_RCPP;
    }
	
    SEXP isNullExtPtr(SEXP Ptr) {
	BEGIN_RCPP;
	return ::Rf_ScalarLogical(XPtr<lmResp>(Ptr) == (lmResp*)NULL);
	END_RCPP;
    }

    // linear model response (also the base class for other response classes)

    SEXP lm_Create(SEXP y, SEXP weights, SEXP offset, SEXP mu,
		   SEXP sqrtXwt, SEXP sqrtrwt, SEXP wtres) {
	BEGIN_RCPP;
	lmResp *ans = new lmResp(y, weights, offset, mu, sqrtXwt, sqrtrwt, wtres);
	return wrap(XPtr<lmResp>(ans, true));
	END_RCPP;
    }

    SEXP lm_setOffset(SEXP ptr_, SEXP offset) {
	BEGIN_RCPP;
	XPtr<lmResp>(ptr_)->setOffset(as<MVec>(offset));
	END_RCPP;
    }

    SEXP lm_setWeights(SEXP ptr_, SEXP weights) {
	BEGIN_RCPP;
	XPtr<lmResp>(ptr_)->setWeights(as<MVec>(weights));
	END_RCPP;
    }

    SEXP lm_wrss(SEXP ptr_) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lmResp>(ptr_)->wrss());
	END_RCPP;
    }

    SEXP lm_updateMu(SEXP ptr_, SEXP gamma) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lmerResp>(ptr_)->updateMu(as<MVec>(gamma)));
	END_RCPP;
    }

    // linear mixed-effects model response

    SEXP lmer_Create(SEXP y, SEXP weights, SEXP offset, SEXP mu,
		     SEXP sqrtXwt, SEXP sqrtrwt, SEXP wtres) {
	BEGIN_RCPP;
	lmerResp *ans = new lmerResp(y, weights, offset, mu, sqrtXwt, sqrtrwt, wtres);
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

    SEXP lmer_Laplace(SEXP ptr_, SEXP ldL2, SEXP ldRX2, SEXP sqrL) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lmerResp>(ptr_)->Laplace(::Rf_asReal(ldL2),
							     ::Rf_asReal(ldRX2),
							     ::Rf_asReal(sqrL)));
	END_RCPP;
    }

    static double lmer_dev(XPtr<merPredD> ppt, XPtr<lmerResp> rpt, const Eigen::VectorXd& theta) {
	ppt->setTheta(theta);
	ppt->updateDecomp();
        rpt->updateMu(ppt->linPred(0.));
        ppt->updateRes(rpt->wtres());
	ppt->solve();
        rpt->updateMu(ppt->linPred(1.));
	return rpt->Laplace(ppt->ldL2(), ppt->ldRX2(), ppt->sqrL(1.));
    }

    SEXP lmer_Deviance(SEXP pptr_, SEXP rptr_, SEXP theta_) {
	BEGIN_RCPP;
	XPtr<lmerResp>   rpt(rptr_);
	XPtr<merPredD>   ppt(pptr_);
	return ::Rf_ScalarReal(lmer_dev(ppt, rpt, as<MVec>(theta_)));
	END_RCPP;
    }

    SEXP lmer_opt1(SEXP pptr_, SEXP rptr_, SEXP lower_, SEXP upper_) {
	BEGIN_RCPP;
	XPtr<lmerResp>     rpt(rptr_);
	XPtr<merPredD>     ppt(pptr_);
	Eigen::VectorXd     th(1);
	optimizer::Golden gold(::Rf_asReal(lower_), ::Rf_asReal(upper_));
	for (int i = 0; i < 30; ++i) {
	    th[0] = gold.xeval();
	    gold.newf(lmer_dev(ppt, rpt, th));
	}
	return List::create(Named("theta") = ::Rf_ScalarReal(gold.xpos()),
			    Named("objective") = ::Rf_ScalarReal(gold.value()));
	END_RCPP;
    }

    // dense predictor module for mixed-effects models

    SEXP merPredDCreate(SEXP Xs, SEXP Lambdat, SEXP LamtUt, SEXP Lind,
			SEXP RZX, SEXP Ut, SEXP Utr, SEXP V, SEXP VtV,
			SEXP Vtr, SEXP Xwts, SEXP Zt, SEXP beta0,
			SEXP delb, SEXP delu, SEXP theta, SEXP u0) {
	BEGIN_RCPP;
	merPredD *ans = new merPredD(Xs, Lambdat, LamtUt, Lind, RZX, Ut, Utr, V, VtV,
				     Vtr, Xwts, Zt, beta0, delb, delu, theta, u0);
	return wrap(XPtr<merPredD>(ans, true));
	END_RCPP;
    }

				// setters
    SEXP merPredDsetTheta(SEXP ptr, SEXP theta) {
	BEGIN_RCPP;
	XPtr<merPredD>(ptr)->setTheta(as<MVec>(theta));
	return theta;
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

    SEXP merPredDupdateDecomp(SEXP ptr) {
	BEGIN_RCPP;
	XPtr<merPredD>(ptr)->updateDecomp();
	END_RCPP;
    }

    SEXP merPredDupdateL(SEXP ptr) {
	BEGIN_RCPP;
	XPtr<merPredD>(ptr)->updateL();
	END_RCPP;
    }

    SEXP merPredDupdateLamtUt(SEXP ptr) {
	BEGIN_RCPP;
	XPtr<merPredD>(ptr)->updateLamtUt();
	END_RCPP;
    }

    SEXP merPredDupdateRes(SEXP ptr, SEXP wtres) {
	BEGIN_RCPP;
	XPtr<merPredD>(ptr)->updateRes(as<MVec>(wtres));
	END_RCPP;
    }

    SEXP merPredDupdateXwts(SEXP ptr, SEXP wts) {
	BEGIN_RCPP;
	XPtr<merPredD>(ptr)->updateXwts(as<MVec>(wts));
	END_RCPP;
    }

    SEXP NelderMead_Create(SEXP lb_, SEXP ub_, SEXP xstep0_, SEXP x_, SEXP xtol_) {
	BEGIN_RCPP;
	MVec  lb(as<MVec>(lb_)), ub(as<MVec>(ub_)), xstep0(as<MVec>(xstep0_)), x(as<MVec>(x_)), xtol(as<MVec>(xtol_));
	Nelder_Mead *ans =
	    new Nelder_Mead(lb, ub, xstep0, x, optimizer::nl_stop(as<MVec>(xtol_)));
	return wrap(XPtr<Nelder_Mead>(ans, true));
	END_RCPP;
    }

    SEXP NelderMead_newf(SEXP ptr_, SEXP f_) {
	BEGIN_RCPP;
	switch (XPtr<Nelder_Mead>(ptr_)->newf(::Rf_asReal(f_))) {
	case optimizer::nm_evals:         return ::Rf_ScalarInteger(-4);
	case optimizer::nm_forced:        return ::Rf_ScalarInteger(-3);
	case optimizer::nm_nofeasible:    return ::Rf_ScalarInteger(-2);
	case optimizer::nm_x0notfeasible: return ::Rf_ScalarInteger(-1);
	case optimizer::nm_active:        return ::Rf_ScalarInteger(0);
	case optimizer::nm_minf_max:      return ::Rf_ScalarInteger(1);
	case optimizer::nm_fcvg:          return ::Rf_ScalarInteger(2);
	case optimizer::nm_xcvg:          return ::Rf_ScalarInteger(3);
	}
	END_RCPP;
    }

    SEXP NelderMead_setForce_stop(SEXP ptr_, SEXP stp_) {
	BEGIN_RCPP;
	XPtr<Nelder_Mead>(ptr_)->setForce_stop(::Rf_asLogical(stp_));
	END_RCPP;
    }

    SEXP NelderMead_setFtol_abs(SEXP ptr_, SEXP fta_) {
	BEGIN_RCPP;
	XPtr<Nelder_Mead>(ptr_)->setFtol_rel(::Rf_asReal(fta_));
	END_RCPP;
    }

    SEXP NelderMead_setFtol_rel(SEXP ptr_, SEXP ftr_) {
	BEGIN_RCPP;
	XPtr<Nelder_Mead>(ptr_)->setFtol_rel(::Rf_asReal(ftr_));
	END_RCPP;
    }

    SEXP NelderMead_setMaxeval(SEXP ptr_, SEXP mm_) {
	BEGIN_RCPP;
	XPtr<Nelder_Mead>(ptr_)->set_Maxeval(::Rf_asInteger(mm_));
	END_RCPP;
    }
	
    SEXP NelderMead_setMinf_max(SEXP ptr_, SEXP mm_) {
	BEGIN_RCPP;
	XPtr<Nelder_Mead>(ptr_)->setMinf_max(::Rf_asReal(mm_));
	END_RCPP;
    }
	
    SEXP NelderMead_xeval(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<Nelder_Mead>(ptr_)->xeval());
	END_RCPP;
    }

    SEXP NelderMead_value(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<Nelder_Mead>(ptr_)->value());
	END_RCPP;
    }

    SEXP NelderMead_xpos(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<Nelder_Mead>(ptr_)->xpos());
	END_RCPP;
    }

				// return the number of function evaluations performed
    SEXP NelderMead_evals(SEXP ptr_) { 
	BEGIN_RCPP;
	return wrap(XPtr<Nelder_Mead>(ptr_)->evals());
	END_RCPP;
    }
	
    // nonlinear model response (also the base class for other response classes)

    SEXP nls_Create(SEXP y, SEXP weights, SEXP offset, SEXP mu, SEXP sqrtXwt,
		    SEXP sqrtrwt, SEXP wtres, SEXP gamma, SEXP mod, SEXP env, SEXP pnms) {
	BEGIN_RCPP;
	nlsResp *ans =
	    new nlsResp(y, weights, offset, mu, sqrtXwt, sqrtrwt, wtres, gamma, mod, env, pnms);
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
	return ::Rf_ScalarReal(XPtr<nlsResp>(ptr_)->updateMu(as<MVec>(gamma)));
	END_RCPP;
    }
}

#include <R_ext/Rdynload.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {

    CALLDEF(Eigen_SSE, 0),

    CALLDEF(glm_Create, 10),	// generate external pointer

    CALLDEF(glm_setN, 2),	// setters

    CALLDEF(glm_devResid, 1),	// getters
    CALLDEF(glm_family, 1),
    CALLDEF(glm_link, 1),
    CALLDEF(glm_muEta, 1),
    CALLDEF(glm_resDev, 1),
    CALLDEF(glm_sqrtWrkWt, 1),
    CALLDEF(glm_variance, 1),
    CALLDEF(glm_wrkResids, 1),
    CALLDEF(glm_wrkResp, 1),

    CALLDEF(glm_Laplace, 4),	// methods
    CALLDEF(glm_updateMu, 2),
    CALLDEF(glm_updateWts, 1),

    CALLDEF(glmFamily_Create, 1), // generate external pointer

    CALLDEF(glmFamily_link, 2),	// methods
    CALLDEF(glmFamily_linkInv, 2),
    CALLDEF(glmFamily_devResid, 4),
    CALLDEF(glmFamily_muEta, 2),
    CALLDEF(glmFamily_variance, 2),

    CALLDEF(glmerAGQ, 8),
    CALLDEF(glmerPwrssUpdate, 5),
    CALLDEF(glmerWrkIter, 2),
    CALLDEF(glmerLaplace, 8),

    CALLDEF(golden_Create, 2),
    CALLDEF(golden_newf, 2),
    CALLDEF(golden_value, 1),
    CALLDEF(golden_xeval, 1),
    CALLDEF(golden_xpos, 1),

    CALLDEF(isNullExtPtr, 1),

    CALLDEF(lm_Create, 7),	// generate external pointer

    CALLDEF(lm_setOffset, 2),	// setters
    CALLDEF(lm_setWeights, 2),

    CALLDEF(lm_wrss, 1),	// getter

    CALLDEF(lm_updateMu, 2),	// method

    CALLDEF(lmer_Create, 7),	// generate external pointer

    CALLDEF(lmer_setREML, 2),	// setter

    CALLDEF(lmer_Deviance, 3),	// methods
    CALLDEF(lmer_Laplace, 4),
    CALLDEF(lmer_opt1, 4),

    CALLDEF(merPredDCreate, 17), // generate external pointer

    CALLDEF(merPredDsetTheta, 2), // setters

    CALLDEF(merPredDCcNumer, 1), // getters
    CALLDEF(merPredDL, 1),
    CALLDEF(merPredDPvec, 1),
    CALLDEF(merPredDRX, 1),
    CALLDEF(merPredDRXi, 1),
    CALLDEF(merPredDRXdiag, 1),
    CALLDEF(merPredDldL2, 1),
    CALLDEF(merPredDldRX2, 1),
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
    CALLDEF(merPredDupdateL, 1),
    CALLDEF(merPredDupdateLamtUt, 1),
    CALLDEF(merPredDupdateRes, 2),
    CALLDEF(merPredDupdateXwts, 2),

    CALLDEF(NelderMead_Create, 5),
    CALLDEF(NelderMead_newf, 2),
    CALLDEF(NelderMead_setForce_stop, 2),
    CALLDEF(NelderMead_setFtol_abs, 2),
    CALLDEF(NelderMead_setFtol_rel, 2),
    CALLDEF(NelderMead_setMaxeval, 2),
    CALLDEF(NelderMead_setMinf_max, 2),
    CALLDEF(NelderMead_value, 1),
    CALLDEF(NelderMead_xeval, 1),
    CALLDEF(NelderMead_xpos, 1),

    CALLDEF(nlmerLaplace, 8),

    CALLDEF(nls_Create, 11),	// generate external pointer

    CALLDEF(nls_Laplace, 4),	// methods
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

