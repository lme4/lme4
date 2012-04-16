// respModule.cpp: response modules using Eigen
//
// Copyright (C) 2011-2012 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#include "respModule.h"
#include <cmath>

namespace lme4 {
    using Eigen::ArrayXd;
    using Eigen::VectorXd;

    using Rcpp::List;
    using Rcpp::NumericMatrix;
    using Rcpp::as;

    using std::copy;
    using std::invalid_argument;

    typedef Eigen::Map<VectorXd>  MVec;

    lmResp::lmResp(SEXP y, SEXP weights, SEXP offset, SEXP mu, SEXP sqrtXwt,
		   SEXP sqrtrwt, SEXP wtres)
	: d_y(      as<MVec>(y)),
	  d_weights(as<MVec>(weights)),
	  d_offset( as<MVec>(offset)),
	  d_mu(     as<MVec>(mu)),
	  d_sqrtXwt(as<MVec>(sqrtXwt)),
	  d_sqrtrwt(as<MVec>(sqrtrwt)),
	  d_wtres(  as<MVec>(wtres)) {
	updateWrss();
    }

    /** 
     * Update the (conditional) mean response and weighted residuals
     * 
     * @param gamma New value of the linear predictor
     * @return updated sum of squared, weighted residuals
     */
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

    /** 
     * Set a new value of the offset.
     *
     * The values are copied into the d_offset member because that member is mapped. 
     * @param oo New value of the offset
     */
    void lmResp::setOffset(const VectorXd& oo) {
	if (oo.size() != d_offset.size())
	    throw invalid_argument("setOffset: Size mismatch");
	d_offset = oo;		// this copies the values
    }

    void lmResp::setResp(const VectorXd& yy) {
	if (yy.size() != d_y.size())
	    throw invalid_argument("setResp: Size mismatch");
	d_y = yy;
    }

    void lmResp::setWeights(const VectorXd& ww) {
	if (ww.size() != d_weights.size())
	    throw invalid_argument("setWeights: Size mismatch");
	d_weights = ww;
    }

    lmerResp::lmerResp(SEXP y, SEXP weights, SEXP offset, SEXP mu,
		       SEXP sqrtXwt, SEXP sqrtrwt, SEXP wtres)
	: lmResp(y, weights, offset, mu, sqrtXwt, sqrtrwt, wtres),
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
    
    glmResp::glmResp(List fam, SEXP y, SEXP weights, SEXP offset,
		     SEXP mu, SEXP sqrtXwt, SEXP sqrtrwt, SEXP wtres, SEXP eta, SEXP n)
	: lmResp(y, weights, offset, mu, sqrtXwt, sqrtrwt, wtres),
	  d_fam(fam),
	  d_eta(as<MVec>(eta)),
	  d_n(as<MVec>(n)) {
    }

    double  glmResp::aic() const {
	return d_fam.aic(d_y, d_n, d_mu, d_weights, resDev());
    }

    ArrayXd glmResp::devResid() const {
	return d_fam.devResid(d_y, d_mu, d_weights);
    }

    ArrayXd glmResp::muEta() const {
	return d_fam.muEta(d_eta);
    }

    ArrayXd glmResp::variance() const {
	return d_fam.variance(d_mu);
    }

    ArrayXd glmResp::wrkResids() const {
	return (d_y - d_mu).array() / muEta();
    }

    ArrayXd glmResp::wrkResp() const {
	return (d_eta - d_offset).array() + wrkResids();
    }

    ArrayXd glmResp::wtWrkResp() const {
	return wrkResp() * sqrtWrkWt();
    }

    ArrayXd glmResp::sqrtWrkWt() const {
	return muEta() * (d_weights.array() / variance()).sqrt();
    }

    double glmResp::Laplace(double ldL2, double ldRX2, double sqrL) const {
	return ldL2 + sqrL + aic();
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
	d_sqrtrwt = d_weights.array() / variance().sqrt();
	d_sqrtXwt = muEta() * d_sqrtrwt.array();
	return updateWrss();
    }

    void glmResp::setN(const VectorXd& n) {
	if (n.size() != d_n.size())
	    throw invalid_argument("n size mismatch");
	d_n = n;
    }

    nlsResp::nlsResp(SEXP y, SEXP weights, SEXP offset, SEXP mu, SEXP sqrtXwt,
		     SEXP sqrtrwt, SEXP wtres, SEXP gamma, SEXP mm, SEXP ee,
		     SEXP pp)
	: lmResp(y, weights, offset, mu, sqrtXwt, sqrtrwt, wtres),
	  d_gamma(as<MVec>(gamma)),
	  d_nlenv(as<Environment>(ee)),
	  d_nlmod(as<Language>(mm)),
	  d_pnames(as<CharacterVector>(pp)) {
    }

    double nlsResp::Laplace(double ldL2, double ldRX2, double sqrL) const {
	double lnum = 2.* PI * (d_wrss + sqrL), n = d_y.size();
	return ldL2 + n * (1 + std::log(lnum / n));
    }

    double nlsResp::updateMu(const VectorXd& gamma) {
	int             n = d_y.size();
	if (gamma.size() != d_gamma.size())
	    throw invalid_argument("size mismatch in updateMu");
	std::copy(gamma.data(), gamma.data() + gamma.size(), d_gamma.data());
	const VectorXd lp(d_gamma + d_offset); // linear predictor
	const double  *gg = lp.data();

	for (int p = 0; p < d_pnames.size(); ++p) {
	    std::string pn(d_pnames[p]);
	    NumericVector pp = d_nlenv.get(pn);
	    std::copy(gg + n * p, gg + n * (p + 1), pp.begin());
	}
	NumericVector  rr = d_nlmod.eval(SEXP(d_nlenv));
	if (rr.size() != n)
	    throw invalid_argument("dimension mismatch");
	std::copy(rr.begin(), rr.end(), d_mu.data());
	NumericMatrix  gr = rr.attr("gradient");
	std::copy(gr.begin(), gr.end(), d_sqrtXwt.data());
	return updateWrss();
    }

}
