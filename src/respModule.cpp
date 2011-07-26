// respModule.cpp: response modules using Eigen
//
// Copyright (C)       2011 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#include "respModule.h"

using namespace Rcpp;
using namespace std;

namespace lme4Eigen {

#if 0				// no longer needed
    modResp::modResp(S4 xp)                                       throw (invalid_argument)
	: d_yR(xp.slot("y")),
	  d_y(d_yR.begin(), d_yR.size()),
	  d_weights(d_y.size()),
	  d_offset( NumericVector(xp.slot("offset")).size()),
	  d_mu(     d_y.size()),
	  d_sqrtXwt(d_offset.size()),
	  d_sqrtrwt(d_y.size()),
	  d_wtres(  d_y.size()) {
	NumericVector weights(xp.slot("weights")), offset(xp.slot("offset"));
	int n = d_y.size(), os = d_offset.size();
	if (weights.size() != n) throw invalid_argument("size mismatch of y and weights");
	copy(weights.begin(), weights.end(), d_weights.begin());
	if (os < 1 || os % n)
	    throw invalid_argument("length(offset) must be a positive multiple of length(y)");
	copy(offset.begin(), offset.end(), d_offset.begin());
	init();
    }
#endif

    modResp::modResp(NumericVector y)                             throw (invalid_argument)
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

#if 0    			// no longer needed
    modResp::modResp(NumericVector y, 
		     NumericVector weights)                       throw (invalid_argument)
	: d_yR(y),
	  d_y(d_yR.begin(), d_yR.size()),
	  d_weights(y.size()),
	  d_offset( VectorXd::Zero(y.size())),
	  d_mu(     y.size()),
	  d_sqrtXwt(y.size()),
	  d_sqrtrwt(y.size()),
	  d_wtres(  y.size()) {
	if (weights.size() != y.size())
	    throw invalid_argument("lengths of y and wts must agree");
	copy(weights.begin(), weights.end(), d_weights.begin());
	init();
    }

    modResp::modResp(NumericVector y,
		     NumericVector weights,
		     NumericVector offset)                        throw (invalid_argument)
	: d_yR(y),
	  d_y(d_yR.begin(), d_yR.size()),
	  d_weights(y.size()),
	  d_offset( offset.size()),
	  d_mu(     y.size()),
	  d_sqrtXwt(offset.size()),
	  d_sqrtrwt(y.size()),
	  d_wtres(  y.size()) {
	int n = y.size(), os = offset.size();
	if (weights.size() != n) throw invalid_argument("lengths of y and wts must agree");
	copy(weights.begin(), weights.end(), d_weights.begin());
	if (os < 1 || os % n)
	    throw invalid_argument("length(offset) must be a positive multiple of length(y)");
	copy(offset.begin(), offset.end(), d_offset.begin());
	init();
    }

    void modResp::init() {
	d_sqrtrwt = d_weights.cwiseSqrt();
	if (d_mu.size() == d_offset.size()) {	
	    d_mu = d_offset;
	    copy(d_sqrtrwt.begin(), d_sqrtrwt.end(), d_sqrtXwt.begin());
	}
	updateWrss();
    }
#endif
    /** 
     * Update the wtres vector and return its sum of squares
     *   wtres <- sqrtrwt * (y - mu)
     *   return(wrss <- sum(wtres^2))
     *
     * @return Updated weighted residual sum of squares
     */
    double modResp::updateWrss() {
	d_wtres = d_sqrtrwt.cwiseProduct(d_y - d_mu);
	d_wrss  = d_wtres.squaredNorm();
	return d_wrss;
    }

    void modResp::setOffset(const NumericVector& oo)              throw (invalid_argument) {
	if (oo.size() != d_offset.size())
	    throw invalid_argument("setOffset: Size mismatch");
	copy(oo.begin(), oo.end(), d_offset.data());
    }

    void modResp::setWeights(const NumericVector& ww)             throw (invalid_argument) {
	if (ww.size() != d_weights.size())
	    throw invalid_argument("setWeights: Size mismatch");
	copy(ww.begin(), ww.end(), d_weights.data());
    }
#if 0 				// no longer needed
    lmerResp::lmerResp(S4 xp)                                     throw (invalid_argument)
	: modResp(xp),
	  d_reml( ::Rf_asInteger(xp.slot("REML"))) {
    }
#endif    
    lmerResp::lmerResp(NumericVector y)                   throw (invalid_argument)
	: modResp(y),
	  d_reml(0) {
    }
#if 0    			// no longer needed
    lmerResp::lmerResp(int rr, NumericVector y, 
		       NumericVector weights)                     throw (invalid_argument)
	: modResp(y, weights),
	  d_reml(rr) {
    }

    lmerResp::lmerResp(int rr, NumericVector y,
		       NumericVector weights,
		       NumericVector offset)               throw (invalid_argument)
	: modResp(y, weights, offset),
	  d_reml(rr) {
    }
#endif    
    double lmerResp::Laplace(double ldL2, double ldRX2, double sqrL) const {
	double lnum = 2.* PI * (d_wrss + sqrL), n = d_y.size();
	if (d_reml == 0) return ldL2 + n * (1. + log(lnum / n));
	double nmp = n - d_reml;
	return ldL2 + ldRX2 + nmp * (1. + log(lnum / nmp));
    }
#if 0    
    void lmerResp::setWeights(const VectorXd& ww)                 throw (invalid_argument) {
	if (ww.size() != d_weights.size())
	    throw invalid_argument("setWeights: Size mismatch");
	copy(ww.begin(), ww.end(), d_weights.begin());
    }
    
    void lmerResp::setOffset(const VectorXd& oo)                  throw (invalid_argument) {
	if (oo.size() != d_offset.size())
	    throw invalid_argument("setOffset: Size mismatch");
	copy(oo.begin(), oo.end(), d_offset.begin());
    }
#endif    
    void lmerResp::setReml(int rr)                                throw (invalid_argument) {
	if (rr < 0) throw invalid_argument("setReml: negative rr");
	d_reml = rr;
    }
    
    double lmerResp::updateMu(const VectorXd& gamma) {
	d_mu = d_offset + gamma;
	return updateWrss();
    }
#if 0    
    glmerResp::glmerResp(S4 xp)                                   throw (invalid_argument)
	: modResp(xp),
	  d_fam(  SEXP(xp.slot("family"))),
	  d_eta(  d_y.size()),
	  d_n(    d_y.size()) {
	int n = d_y.size();
	NumericVector nn(xp.slot("n")), eta(xp.slot("eta"));
	if (nn.size() != n || eta.size() != n)
	    throw invalid_argument("y, n and eta must all be the same length");
	updateWts();
    }
#endif
    glmerResp::glmerResp(List fam, NumericVector y)               throw (invalid_argument)
	: modResp(y),
	  d_fam(fam),
	  d_eta(y.size()),
	  d_n(VectorXd::Constant(y.size(), 1.)) {
    }
#if 0
    glmerResp::glmerResp(List fam, NumericVector y,
			 NumericVector weights)                   throw (invalid_argument)
	: modResp(y, weights),
	  d_fam(fam),
	  d_eta(y.size()),
	  d_n(VectorXd::Constant(y.size(), 1.)) {
    }
    
    glmerResp::glmerResp(List fam, NumericVector y,
			 NumericVector weights,
			 NumericVector offset)                    throw (invalid_argument) 
	: modResp(y, weights, offset),
	  d_fam(  fam),
	  d_eta(  y.size()),
	  d_n(    VectorXd::Constant(y.size(), 1.)) {
    }

    glmerResp::glmerResp(List fam, NumericVector y,
			 NumericVector weights,
			 NumericVector offset,
			 NumericVector n)                         throw (invalid_argument)
 	: modResp(y, weights, offset), 
	  d_fam(  fam),
	  d_eta(  y.size()),
    	  d_n(    y.size()) {
	if (n.size() != y.size())
	    throw invalid_argument("lengths of y and n must agree");
	copy(n.begin(), n.end(), d_n.begin());
    }

    glmerResp::glmerResp(List fam, NumericVector y,
			 NumericVector weights,
			 NumericVector offset,
			 NumericVector n, 
			 NumericVector eta)                       throw (invalid_argument) 
	: modResp(y, weights, offset),
	  d_fam(  fam),
	  d_eta(  y.size()),
	  d_n(    y.size()) {
	int nn = y.size();
	if (n.size() != nn || eta.size() != nn )
	    throw invalid_argument("lengths of y, n and eta must agree");
	copy(n.begin(), n.end(), d_n.begin());
	copy(eta.begin(), eta.end(), d_eta.begin());
    }
#endif 
    double glmerResp::updateWts() {
	d_sqrtrwt = d_weights.cwiseQuotient(variance()).cwiseSqrt();
	d_sqrtXwt = muEta() * d_sqrtrwt;
	return updateWrss();
    }

    VectorXd glmerResp::sqrtWrkWt() const {
	VectorXd me = muEta();
	return d_weights.cwiseProduct(me).cwiseProduct(me).cwiseQuotient(variance()).cwiseSqrt();
    }

    double glmerResp::Laplace(double ldL2, double ldRX2, double sqrL) const {
	return ldL2 + sqrL + resDev();
    }

    double glmerResp::updateMu(const VectorXd& gamma) {
	d_eta = d_offset + gamma; // lengths are checked here
	d_mu  = d_fam.linkInv(d_eta);
	return updateWrss();
    }
#if 0
    void glmerResp::setWeights(const VectorXd& ww)                throw (invalid_argument) {
	if (ww.size() != d_weights.size())
	    throw invalid_argument("setWeights: Size mismatch");
	copy(ww.begin(), ww.end(), d_weights.begin());
    }
    
    void glmerResp::setOffset(const VectorXd& oo)                 throw (invalid_argument) {
	if (oo.size() != d_offset.size())
	    throw invalid_argument("setOffset: Size mismatch");
	copy(oo.begin(), oo.end(), d_offset.begin());
    }

    nlmerResp::nlmerResp(S4 xp)                                   throw (invalid_argument)
	: modResp(xp),
	  d_nlenv(SEXP(xp.slot("nlenv"))),
	  d_nlmod(SEXP(xp.slot("nlmod"))),
	  d_pnames(SEXP(xp.slot("pnames"))) {
    }
#endif    

    double nlmerResp::Laplace(double ldL2, double ldRX2, double sqrL) const {
	double lnum = 2.* PI * (d_wrss + sqrL), n = d_y.size();
	return ldL2 + n * (1 + log(lnum / n));
    }

    double nlmerResp::updateMu(VectorXd const &gamma)             throw (invalid_argument) {
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
    SEXP lmerRespCreate(SEXP ys) {
	BEGIN_RCPP;
	lme4Eigen::lmerResp *ans = new lme4Eigen::lmerResp(NumericVector(ys));
	return wrap(XPtr<lme4Eigen::lmerResp>(ans, true));
	END_RCPP;
    }
    
    SEXP lmerRespLaplace(SEXP ptr_, SEXP ldL2, SEXP ldRX2, SEXP sqrL) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lme4Eigen::lmerResp>(ptr_)->Laplace(::Rf_asReal(ldL2),
									::Rf_asReal(ldRX2),
									::Rf_asReal(sqrL)));
	END_RCPP;
    }

    SEXP lmerRespsetREML(SEXP ptr_, SEXP REML) {
	BEGIN_RCPP;
	XPtr<lme4Eigen::lmerResp>(ptr_)->setReml(::Rf_asInteger(REML));
	END_RCPP;
    }

    SEXP lmerRespREML(SEXP ptr_) {
	BEGIN_RCPP;
	return ::Rf_ScalarInteger(XPtr<lme4Eigen::lmerResp>(ptr_)->REML());
	END_RCPP;
    }

    SEXP lmerRespupdateMu(SEXP ptr_, SEXP gamma) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lme4Eigen::lmerResp>(ptr_)->updateMu(as<VectorXd>(gamma)));
	END_RCPP;
    }

    SEXP modRespwrss(SEXP ptr_) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lme4Eigen::modResp>(ptr_)->wrss());
	END_RCPP;
    }
    
    SEXP modRespwtres(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::modResp>(ptr_)->wtres());
	END_RCPP;
    }
    
    SEXP modRespsetOffset(SEXP ptr_, SEXP offset) {
	BEGIN_RCPP;
	XPtr<lme4Eigen::modResp>(ptr_)->setOffset(NumericVector(offset));
	END_RCPP;
    }
    
    SEXP modRespsetWeights(SEXP ptr_, SEXP weights) {
	BEGIN_RCPP;
	XPtr<lme4Eigen::modResp>(ptr_)->setWeights(NumericVector(weights));
	END_RCPP;
    }
    
    SEXP modRespy(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::modResp>(ptr_)->y());
	END_RCPP;
    }
    
    SEXP modRespmu(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::modResp>(ptr_)->mu());
	END_RCPP;
    }
    
    SEXP modRespoffset(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::modResp>(ptr_)->offset());
	END_RCPP;
    }

    SEXP modRespweights(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::modResp>(ptr_)->weights());
	END_RCPP;
    }

}
