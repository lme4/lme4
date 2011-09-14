// respModule.cpp: response modules using Eigen
//
// Copyright (C)       2011 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#include "respModule.h"
#include <cmath>

namespace lme4Eigen {
    using Rcpp::List;
    using Rcpp::NumericMatrix;
    using std::copy;
    using std::string;
    using std::invalid_argument;

    lmResp::lmResp(NumericVector y)
	: d_yR(y),
	  d_y(d_yR.begin(), d_yR.size()),
	  d_weights(VectorXd::Constant(y.size(), 1.0)),
	  d_offset( VectorXd::Zero(y.size())),
	  d_mu(     y.size()),
	  d_sqrtXwt(y.size()),
	  d_sqrtrwt(y.size()),
	  d_wtres(  y.size()) {
	d_sqrtrwt = d_weights.cwiseSqrt();
	if (d_mu.size() == d_offset.size()) {	
	    d_mu = d_offset;
	    copy(d_sqrtrwt.begin(), d_sqrtrwt.end(), d_sqrtXwt.begin());
	}
	updateWrss();
    }

    double lmResp::updateMu(const VectorXd& gamma) {
	if (gamma.size() != d_offset.size())
	    throw invalid_argument("updateMu: Size mismatch");
	d_mu = d_offset + gamma;
	return updateWrss();
    }

    /** 
     * Update the wtres vector and return its sum of squares
     *   wtres <- sqrtrwt * (y - mu)
     *   return(wrss <- sum(wtres^2))
     *
     * @return Updated weighted residual sum of squares
     */
    double lmResp::updateWrss() {
	d_wtres = d_sqrtrwt.cwiseProduct(d_y - d_mu);
	d_wrss  = d_wtres.squaredNorm();
	return d_wrss;
    }

    void lmResp::setOffset(const VectorXd& oo) {
	if (oo.size() != d_offset.size())
	    throw invalid_argument("setOffset: Size mismatch");
	d_offset = oo;
    }

    void lmResp::setWeights(const VectorXd& ww) {
	if (ww.size() != d_weights.size())
	    throw invalid_argument("setWeights: Size mismatch");
	d_weights = ww;
    }

    lmerResp::lmerResp(NumericVector y)
	: lmResp(y),
	  d_reml(0) {
    }

    double lmerResp::Laplace(double ldL2, double ldRX2, double sqrL) const {
	double lnum = std::log(2.* M_PI * (d_wrss + sqrL));
	if (d_reml == 0) return ldL2 + d_y.size() * (1. + lnum - std::log(d_y.size()));
	double nmp = d_y.size() - d_reml;
	return ldL2 + ldRX2 + nmp * (1. + lnum - std::log(nmp));
    }

    void lmerResp::setReml(int rr) {
	if (rr < 0) throw invalid_argument("setReml: negative value for REML not meaningful");
	d_reml = rr;
    }
    
    glmResp::glmResp(List fam, NumericVector y)
	: lmResp(y),
	  d_fam(fam),
	  d_eta(y.size()),
	  d_n(VectorXd::Constant(y.size(), 1.)) {
    }

    VectorXd glmResp::devResid() const {
	return d_fam.devResid(d_mu, d_weights, d_y);
    }

    VectorXd glmResp::muEta() const {
	return d_fam.muEta(d_eta);
    }

    VectorXd glmResp::variance() const {
	return d_fam.variance(d_mu);
    }

    VectorXd glmResp::wrkResids() const {
	return (d_y - d_mu).cwiseQuotient(muEta());
    }

    VectorXd glmResp::wrkResp() const {
	return (d_eta - d_offset) + wrkResids();
    }

    VectorXd glmResp::sqrtWrkWt() const {
	VectorXd me = muEta();
	return d_weights.cwiseProduct(me).cwiseProduct(me).cwiseQuotient(variance()).cwiseSqrt();
    }

    double glmResp::Laplace(double ldL2, double ldRX2, double sqrL) const {
	return ldL2 + sqrL + resDev();
    }

    double glmResp::resDev() const {
	return devResid().sum();
    }

    double glmResp::updateMu(const VectorXd& gamma) {
	d_eta = d_offset + gamma; // lengths are checked here
	d_mu  = d_fam.linkInv(d_eta);
	return updateWrss();
    }

    double glmResp::updateWts() {
	d_sqrtrwt = d_weights.cwiseQuotient(variance()).cwiseSqrt();
	d_sqrtXwt = muEta().cwiseProduct(d_sqrtrwt);
	return updateWrss();
    }

    void glmResp::setN(const VectorXd& n) {
	if (n.size() != d_n.size())
	    throw invalid_argument("n size mismatch");
	d_n = n;
    }

    nlmerResp::nlmerResp(NumericVector ys, Language mm, Environment ee, CharacterVector pp)
	: lmResp(ys), d_nlenv(ee), d_nlmod(mm), d_pnames(pp) {
    }

    double nlmerResp::Laplace(double ldL2, double ldRX2, double sqrL) const {
	double lnum = 2.* PI * (d_wrss + sqrL), n = d_y.size();
	return ldL2 + n * (1 + log(lnum / n));
    }

    double nlmerResp::updateMu(VectorXd const &gamma) {
	int             n = d_y.size();
	VectorXd      gam = gamma + d_offset;
	const double  *gg = gam.begin();

	for (int p = 0; p < d_pnames.size(); p++) {
	    string pn(d_pnames[p]);
	    NumericVector pp = d_nlenv.get(pn);
	    copy(gg + n * p, gg + n * (p + 1), pp.begin());
	}
	NumericVector  rr = d_nlmod.eval(SEXP(d_nlenv));
	if (rr.size() != n)
	    throw invalid_argument("dimension mismatch");
	copy(rr.begin(), rr.end(), d_mu.begin());
	NumericMatrix rrg = rr.attr("gradient");
	copy(rrg.begin(), rrg.end(), d_sqrtXwt.begin());
	return updateWrss();
    }

}

extern "C" {
    using Eigen::VectorXd;

    using Rcpp::List;
    using Rcpp::NumericVector;
    using Rcpp::XPtr;
    using Rcpp::as;
    using Rcpp::wrap;

    // generalized linear model (and generalized linear mixed model) response

    SEXP glm_Create(SEXP fams, SEXP ys) {
	BEGIN_RCPP;
	lme4Eigen::glmResp *ans = new lme4Eigen::glmResp(List(fams), NumericVector(ys));
	return wrap(XPtr<lme4Eigen::glmResp>(ans, true));
	END_RCPP;
    }
    
    SEXP glm_setN(SEXP ptr_, SEXP n) {
	BEGIN_RCPP;
	XPtr<lme4Eigen::glmResp>(ptr_)->setN(as<VectorXd>(n));
	END_RCPP;
    }

    SEXP glm_devResid(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::glmResp>(ptr_)->devResid());
	END_RCPP;
    }

    SEXP glm_eta(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::glmResp>(ptr_)->eta());
	END_RCPP;
    }

    SEXP glm_family(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::glmResp>(ptr_)->family());
	END_RCPP;
    }

    SEXP glm_link(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::glmResp>(ptr_)->link());
	END_RCPP;
    }

    SEXP glm_muEta(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::glmResp>(ptr_)->muEta());
	END_RCPP;
    }

    SEXP glm_n(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::glmResp>(ptr_)->n());
	END_RCPP;
    }

    SEXP glm_resDev(SEXP ptr_) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lme4Eigen::glmResp>(ptr_)->resDev());
	END_RCPP;
    }

    SEXP glm_sqrtWrkWt(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::glmResp>(ptr_)->sqrtWrkWt());
	END_RCPP;
    }

    SEXP glm_updateWts(SEXP ptr_) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lme4Eigen::glmResp>(ptr_)->updateWts());
	END_RCPP;
    }

    SEXP glm_variance(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::glmResp>(ptr_)->variance());
	END_RCPP;
    }

    SEXP glm_wrkResids(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::glmResp>(ptr_)->wrkResids());
	END_RCPP;
    }

    SEXP glm_wrkResp(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::glmResp>(ptr_)->wrkResp());
	END_RCPP;
    }

    SEXP glm_Laplace(SEXP ptr_, SEXP ldL2, SEXP ldRX2, SEXP sqrL) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lme4Eigen::glmResp>(ptr_)->Laplace(::Rf_asReal(ldL2),
									 ::Rf_asReal(ldRX2),
									 ::Rf_asReal(sqrL)));
	END_RCPP;
    }

    SEXP glm_updateMu(SEXP ptr_, SEXP gamma) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lme4Eigen::glmResp>(ptr_)->updateMu(as<VectorXd>(gamma)));
	END_RCPP;
    }

    // linear model response (also the base class for other response classes)

    SEXP lm_Create(SEXP ys) {
	BEGIN_RCPP;
	lme4Eigen::lmResp *ans = new lme4Eigen::lmResp(NumericVector(ys));
	return wrap(XPtr<lme4Eigen::lmResp>(ans, true));
	END_RCPP;
    }
    
    SEXP lm_setOffset(SEXP ptr_, SEXP offset) {
	BEGIN_RCPP;
	XPtr<lme4Eigen::lmResp>(ptr_)->setOffset(as<VectorXd>(offset));
	END_RCPP;
    }
    
    SEXP lm_setWeights(SEXP ptr_, SEXP weights) {
	BEGIN_RCPP;
	XPtr<lme4Eigen::lmResp>(ptr_)->setWeights(as<VectorXd>(weights));
	END_RCPP;
    }
    
    SEXP lm_mu(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::lmResp>(ptr_)->mu());
	END_RCPP;
    }
    
    SEXP lm_offset(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::lmResp>(ptr_)->offset());
	END_RCPP;
    }

    SEXP lm_sqrtXwt(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::lmResp>(ptr_)->sqrtXwt());
	END_RCPP;
    }

    SEXP lm_sqrtrwt(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::lmResp>(ptr_)->sqrtrwt());
	END_RCPP;
    }

    SEXP lm_weights(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::lmResp>(ptr_)->weights());
	END_RCPP;
    }

    SEXP lm_wrss(SEXP ptr_) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lme4Eigen::lmResp>(ptr_)->wrss());
	END_RCPP;
    }

    SEXP lm_wtres(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::lmResp>(ptr_)->wtres());
	END_RCPP;
    }
    
    SEXP lm_y(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::lmResp>(ptr_)->y());
	END_RCPP;
    }

    SEXP lm_updateMu(SEXP ptr_, SEXP gamma) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lme4Eigen::lmerResp>(ptr_)->updateMu(as<VectorXd>(gamma)));
	END_RCPP;
    }
    
    SEXP lmer_Create(SEXP ys) {
	BEGIN_RCPP;
	lme4Eigen::lmerResp *ans = new lme4Eigen::lmerResp(NumericVector(ys));
	return wrap(XPtr<lme4Eigen::lmerResp>(ans, true));
	END_RCPP;
    }
    
    SEXP lmer_setREML(SEXP ptr_, SEXP REML) {
	BEGIN_RCPP;
	XPtr<lme4Eigen::lmerResp>(ptr_)->setReml(::Rf_asInteger(REML));
	END_RCPP;
    }

    SEXP lmer_REML(SEXP ptr_) {
	BEGIN_RCPP;
	return ::Rf_ScalarInteger(XPtr<lme4Eigen::lmerResp>(ptr_)->REML());
	END_RCPP;
    }
    
    SEXP lmer_Laplace(SEXP ptr_, SEXP ldL2, SEXP ldRX2, SEXP sqrL) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lme4Eigen::lmerResp>(ptr_)->Laplace(::Rf_asReal(ldL2),
									::Rf_asReal(ldRX2),
									::Rf_asReal(sqrL)));
	END_RCPP;
    }

}
