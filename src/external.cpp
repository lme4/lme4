// external.cpp: externally .Call'able functions in lme4
//
// Copyright (C)       2011-2012 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#include "predModule.h"
#include "respModule.h"
#include "optimizer.h"

extern "C" {
    typedef   Eigen::VectorXi        iVec;
    typedef   Eigen::Map<iVec>      MiVec;
    typedef   Eigen::MatrixXd         Mat;
    typedef   Eigen::Map<Mat>        MMat;
    typedef   Eigen::VectorXd         Vec;
    typedef   Eigen::Map<Vec>        MVec;
    typedef   Eigen::ArrayXd          Ar1;
    typedef   Eigen::Map<Ar1>        MAr1;
    typedef   Eigen::ArrayXXd         Ar2;
    typedef   Eigen::Map<Ar2>        MAr2;

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

    using      lme4::glmResp;
    using      lme4::lmResp;
    using      lme4::lmerResp;
    using      lme4::merPredD;
    using      lme4::nlsResp;

    using optimizer::Golden;
    using optimizer::Nelder_Mead;
    using optimizer::nm_status;

    using      std::runtime_error;

    // utilities

    SEXP allPerm_int(SEXP v_) {
	BEGIN_RCPP;
	iVec     v(as<iVec>(v_));   // forces a copy
	int     sz(v.size());
	std::vector<iVec> vec;
	
	std::sort(v.data(), v.data() + sz);
	do {
	    vec.push_back(iVec(v));
	} while (std::next_permutation(v.data(), v.data() + sz));
	
	int  nperm(vec.size());
	List allPerm(nperm);
	for (int j = 0; j < nperm; ++j) allPerm[j] = wrap(vec[j]);
	return allPerm;
	END_RCPP;
    }

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

    SEXP glm_aic(SEXP ptr_) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<glmResp>(ptr_)->aic());
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

    SEXP glm_setTheta(SEXP ptr, SEXP newtheta) {
	BEGIN_RCPP;
	XPtr<glmResp>(ptr)->setTheta(::Rf_asReal(newtheta));
	END_RCPP;
    }

    SEXP glm_sqrtWrkWt(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<glmResp>(ptr_)->sqrtWrkWt());
	END_RCPP;
    }

    SEXP glm_theta(SEXP ptr) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<glmResp>(ptr)->theta());
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

    SEXP glm_wtWrkResp(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<glmResp>(ptr_)->wtWrkResp());
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

    SEXP glmFamily_devResid(SEXP ptr, SEXP y, SEXP mu, SEXP wt) {
	BEGIN_RCPP;
	return wrap(XPtr<glmFamily>(ptr)->devResid(as<MVec>(y),
						   as<MVec>(mu),
						   as<MVec>(wt)));
	END_RCPP;
    }

    SEXP glmFamily_aic(SEXP ptr, SEXP y, SEXP n, SEXP mu, SEXP wt, SEXP dev) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<glmFamily>(ptr)->aic(as<MVec>(y),
							 as<MVec>(n),
							 as<MVec>(mu),
							 as<MVec>(wt),
							 ::Rf_asReal(dev)));
	END_RCPP;
    }

    SEXP glmFamily_muEta(SEXP ptr, SEXP eta) {
	BEGIN_RCPP;
	return wrap(XPtr<glmFamily>(ptr)->muEta(as<MVec>(eta)));
	END_RCPP;
    }

    SEXP glmFamily_setTheta(SEXP ptr, SEXP ntheta) {
	BEGIN_RCPP;
	XPtr<glmFamily>(ptr)->setTheta(::Rf_asReal(ntheta));
	END_RCPP;
    }

    SEXP glmFamily_theta(SEXP ptr) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<glmFamily>(ptr)->theta());
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

    static double internal_glmerWrkIter(merPredD *pp, glmResp *rp, bool uOnly) {
	pp->updateXwts(rp->sqrtWrkWt());
	pp->updateDecomp();
	pp->updateRes(rp->wtWrkResp());
	if (uOnly) pp->solveU();
	else pp->solve();
	rp->updateMu(pp->linPred(1.));
	return rp->resDev() + pp->sqrL(1.);
    }

    static void pwrssUpdate(glmResp *rp, merPredD *pp, bool uOnly, double tol, int verbose) {
	double oldpdev=std::numeric_limits<double>::max();
	bool   cvgd = false, verb = verbose > 2;
	for (int i = 0; i < 30; i++) {
	    Vec   olddelu(pp->delu()), olddelb(pp->delb());
	    double pdev=internal_glmerWrkIter(pp, rp, uOnly);
	    if (verb) Rcpp::Rcout << i << ": " << pdev << std::endl;
	    if (std::abs((oldpdev - pdev) / pdev) < tol) {cvgd = true; break;}
	    if (pdev > oldpdev) { // try step halving
		if (verb) Rcpp::Rcout << "Step halving: oldpdev = "
				      << oldpdev << ", pdev = " << pdev
				      << std::endl;
		for (int k = 0; k < 10 && pdev > oldpdev; k++) {
		    pp->setDelu((olddelu + pp->delu())/2.);
		    if (!uOnly) pp->setDelb((olddelb + pp->delb())/2.);
		    pdev = internal_glmerWrkIter(pp, rp, uOnly);
		    if (verb) Rcpp::Rcout << "k = " << k << ", pdev = "
					  << pdev << std::endl;
		}
		if (pdev > oldpdev) throw runtime_error("PIRLS step failed");
	    }
	    oldpdev = pdev;
	}
	if (!cvgd)
	    throw runtime_error("pwrssUpdate did not converge in 30 iterations");
    }

    SEXP glmerLaplace(SEXP pp_, SEXP rp_, SEXP nAGQ_, SEXP tol_, SEXP verbose_) {
	BEGIN_RCPP;
	XPtr<glmResp>  rp(rp_);
	XPtr<merPredD> pp(pp_);
	pwrssUpdate(rp, pp, ::Rf_asInteger(nAGQ_), ::Rf_asReal(tol_), ::Rf_asInteger(verbose_));
	return ::Rf_ScalarReal(rp->Laplace(pp->ldL2(), pp->ldRX2(), pp->sqrL(1.)));
	END_RCPP;
    }

    static Ar1 devcCol(const MiVec& fac, const Ar1& u, const Ar1& devRes) {
	Ar1  ans(u.square());
	for (int i = 0; i < devRes.size(); ++i) ans[fac[i] - 1] += devRes[i];
	return ans;
    }

    static double sqrt2pi = std::sqrt(2. * PI);

    SEXP glmerAGQ(SEXP pp_, SEXP rp_, SEXP tol_, SEXP GQmat_, SEXP fac_, SEXP verbose_) {
	BEGIN_RCPP;
	
	XPtr<glmResp>     rp(rp_);
	XPtr<merPredD>    pp(pp_);
	const MiVec      fac(as<MiVec>(fac_));
	double           tol(::Rf_asReal(tol_));
	double          verb(::Rf_asReal(verbose_));
	if (fac.size() != rp->mu().size())
	    throw std::invalid_argument("size of fac must match dimension of response vector");

	pwrssUpdate(rp, pp, true, tol, verb); // should be a no-op
	const Ar1      devc0(devcCol(fac, pp->u(1.), rp->devResid()));

	const unsigned int q(pp->u0().size());
	if (pp->L().factor()->nzmax !=  q)
	    throw std::invalid_argument("AGQ only defined for a single scalar random-effects term");
	const Ar1         sd(MAr1((double*)pp->L().factor()->x, q).inverse());

	const MMat     GQmat(as<MMat>(GQmat_));
	Ar1             mult(q);

	mult.setZero();
	for (int i = 0; i < GQmat.rows(); ++i) {
	    double     zknot(GQmat(i, 0));
	    if (zknot == 0)
		mult += Ar1::Constant(q, GQmat(i, 1));
	    else {
		pp->setU0(zknot * sd); // to be added to current delu 
		rp->updateMu(pp->linPred(1.));
		mult += (-0.5 * (devcCol(fac, pp->u(1.), rp->devResid()) - devc0) -
			 GQmat(i, 2)).exp() * GQmat(i, 1)/sqrt2pi;
	    }
	}
	pp->setU0(Vec::Zero(q)); // restore settings from pwrssUpdate;
	rp->updateMu(pp->linPred(1.));
	return ::Rf_ScalarReal(devc0.sum() + pp->ldL2() - 2 * std::log(mult.prod()));
	END_RCPP;
    }

    void nstepFac(nlsResp *rp, merPredD *pp, int verb) {
	double prss0(pwrss(rp, pp, 0.));

	for (double fac = 1.; fac > 0.001; fac /= 2.) {
	    double prss1 = rp->updateMu(pp->linPred(fac)) + pp->sqrL(fac);
	    if (verb > 3)
		::Rprintf("prss0=%10g, diff=%10g, fac=%6.4f\n",
			  prss0, prss0 - prss1, fac);
	    if (prss1 < prss0) {
		pp->installPars(fac);
		return;
	    }
	}
	throw runtime_error("step factor reduced below 0.001 without reducing pwrss");
    }

#define NMAXITER 300
    static void prssUpdate(nlsResp *rp, merPredD *pp, int verb, bool uOnly, double tol) {
	bool cvgd(false);
	for (int it=0; it < NMAXITER; ++it) {
	    rp->updateMu(pp->linPred(0.));
	    pp->updateXwts(rp->sqrtXwt());
	    pp->updateDecomp();
	    pp->updateRes(rp->wtres());
	    double ccrit((uOnly ? pp->solveU() : pp->solve())/pwrss(rp, pp, 0.));
	    if (verb > 3) 
		::Rprintf("ccrit=%10g, tol=%10g\n", ccrit, tol);
	    if (ccrit < tol) {
		cvgd = true;
		break;
	    }
	    nstepFac(rp, pp, verb);
	}
	if (!cvgd) throw runtime_error("prss failed to converge in 300 iterations");
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

    SEXP lm_setResp(SEXP ptr_, SEXP resp) {
	BEGIN_RCPP;
	XPtr<lmResp>(ptr_)->setResp(as<MVec>(resp));
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

    SEXP merPredDsetBeta0(SEXP ptr, SEXP beta0) {
	BEGIN_RCPP;
	XPtr<merPredD>(ptr)->setBeta0(as<MVec>(beta0));
	END_RCPP;
    }


    SEXP merPredDsetDelu(SEXP ptr, SEXP delu) {
	BEGIN_RCPP;
	XPtr<merPredD>(ptr)->setDelu(as<MVec>(delu));
	END_RCPP;
    }


    SEXP merPredDsetDelb(SEXP ptr, SEXP delb) {
	BEGIN_RCPP;
	XPtr<merPredD>(ptr)->setDelb(as<MVec>(delb));
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

    SEXP merPredDcondVar(SEXP ptr, SEXP rho) {
	BEGIN_RCPP;
	return wrap(XPtr<merPredD>(ptr)->condVar(Rcpp::Environment(rho)));
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

    SEXP NelderMead_setIprint(SEXP ptr_, SEXP ip_) {
	BEGIN_RCPP;
	XPtr<Nelder_Mead>(ptr_)->set_Iprint(::Rf_asInteger(ip_));
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
	return wrap(int(XPtr<Nelder_Mead>(ptr_)->evals()));
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

    CALLDEF(allPerm_int, 1),

    CALLDEF(glm_Create, 10),	// generate external pointer

    CALLDEF(glm_setN, 2),	// setters

    CALLDEF(glm_aic,            1), // getters
    CALLDEF(glm_devResid,       1),	
    CALLDEF(glm_family,         1),
    CALLDEF(glm_link,           1),
    CALLDEF(glm_muEta,          1),
    CALLDEF(glm_resDev,         1),
    CALLDEF(glm_setTheta,       2),
    CALLDEF(glm_sqrtWrkWt,      1),
    CALLDEF(glm_theta,          1),
    CALLDEF(glm_variance,       1),
    CALLDEF(glm_wtWrkResp,      1),
    CALLDEF(glm_wrkResids,      1),
    CALLDEF(glm_wrkResp,        1),

    CALLDEF(glm_Laplace,        4), // methods
    CALLDEF(glm_updateMu,       2),
    CALLDEF(glm_updateWts,      1),

    CALLDEF(glmFamily_Create,   1), // generate external pointer

    CALLDEF(glmFamily_aic,      6), // methods
    CALLDEF(glmFamily_link,     2),
    CALLDEF(glmFamily_linkInv,  2),
    CALLDEF(glmFamily_devResid, 4),
    CALLDEF(glmFamily_muEta,    2),
    CALLDEF(glmFamily_setTheta, 2),
    CALLDEF(glmFamily_theta,    1),
    CALLDEF(glmFamily_variance, 2),

    CALLDEF(glmerAGQ,           6),
    CALLDEF(glmerLaplace,       5),

    CALLDEF(golden_Create,      2),
    CALLDEF(golden_newf,        2),
    CALLDEF(golden_value,       1),
    CALLDEF(golden_xeval,       1),
    CALLDEF(golden_xpos,        1),

    CALLDEF(isNullExtPtr, 1),

    CALLDEF(lm_Create,          7), // generate external pointer

    CALLDEF(lm_setOffset,       2), // setters
    CALLDEF(lm_setResp,         2),
    CALLDEF(lm_setWeights,      2),

    CALLDEF(lm_wrss,            1), // getter

    CALLDEF(lm_updateMu,        2), // method

    CALLDEF(lmer_Create,        7), // generate external pointer

    CALLDEF(lmer_setREML,       2), // setter

    CALLDEF(lmer_Deviance,      3), // methods
    CALLDEF(lmer_Laplace,       4),
    CALLDEF(lmer_opt1,          4),

    CALLDEF(merPredDCreate,    17), // generate external pointer

    CALLDEF(merPredDsetTheta,   2), // setters
    CALLDEF(merPredDsetBeta0,   2), 

    CALLDEF(merPredDsetDelu,    2), // setters
    CALLDEF(merPredDsetDelb,    2), 

    CALLDEF(merPredDCcNumer,    1), // getters
    CALLDEF(merPredDL,          1),
    CALLDEF(merPredDPvec,       1),
    CALLDEF(merPredDRX,         1),
    CALLDEF(merPredDRXdiag,     1),
    CALLDEF(merPredDRXi,        1),
    CALLDEF(merPredDldL2,       1),
    CALLDEF(merPredDldRX2,      1),
    CALLDEF(merPredDunsc,       1),

    CALLDEF(merPredDb,          2), // methods
    CALLDEF(merPredDbeta,       2),
    CALLDEF(merPredDcondVar,    2),
    CALLDEF(merPredDlinPred,    2),
    CALLDEF(merPredDinstallPars,2),
    CALLDEF(merPredDsolve,      1),
    CALLDEF(merPredDsolveU,     1),
    CALLDEF(merPredDsqrL,       2),
    CALLDEF(merPredDu,          2),
    CALLDEF(merPredDupdateDecomp,1),
    CALLDEF(merPredDupdateL,    1),
    CALLDEF(merPredDupdateLamtUt,1),
    CALLDEF(merPredDupdateRes,  2),
    CALLDEF(merPredDupdateXwts, 2),

    CALLDEF(NelderMead_Create,  5),
    CALLDEF(NelderMead_newf,    2),
    CALLDEF(NelderMead_setForce_stop, 2),
    CALLDEF(NelderMead_setFtol_abs, 2),
    CALLDEF(NelderMead_setFtol_rel, 2),
    CALLDEF(NelderMead_setIprint, 2),
    CALLDEF(NelderMead_setMaxeval, 2),
    CALLDEF(NelderMead_setMinf_max, 2),
    CALLDEF(NelderMead_value,   1),
    CALLDEF(NelderMead_xeval,   1),
    CALLDEF(NelderMead_xpos,    1),

    CALLDEF(nlmerLaplace,       8),

    CALLDEF(nls_Create,        11), // generate external pointer

    CALLDEF(nls_Laplace,        4), // methods
    CALLDEF(nls_updateMu,       2),

    {NULL, NULL, 0}
};

/** Initializer for lme4, called upon loading the package.
 *
 *  Register routines that can be called directly from R.
 *  Initialize CHOLMOD and require the LL' form of the factorization.
 *  Install the symbols to be used by functions in the package.
 */
extern "C"
void R_init_lme4(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, (Rboolean)FALSE);
}

