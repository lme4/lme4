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

    modResp::modResp(NumericVector y)
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

    void modResp::setOffset(const NumericVector& oo) {
	if (oo.size() != d_offset.size())
	    throw invalid_argument("setOffset: Size mismatch");
	copy(oo.begin(), oo.end(), d_offset.data());
    }

    void modResp::setWeights(const NumericVector& ww) {
	if (ww.size() != d_weights.size())
	    throw invalid_argument("setWeights: Size mismatch");
	copy(ww.begin(), ww.end(), d_weights.data());
    }

    lmerResp::lmerResp(NumericVector y)
	: modResp(y),
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
    
    double lmerResp::updateMu(const VectorXd& gamma) {
	d_mu = d_offset + gamma;
	return updateWrss();
    }

    glmerResp::glmerResp(List fam, NumericVector y)
	: modResp(y),
	  d_fam(fam),
	  d_eta(y.size()),
	  d_n(VectorXd::Constant(y.size(), 1.)) {
    }

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
    using Rcpp::List;
    using Rcpp::NumericVector;
    using Rcpp::XPtr;
    using Rcpp::as;
    using Rcpp::wrap;

    SEXP glmerRespCreate(SEXP fams, SEXP ys) {
	BEGIN_RCPP;
	lme4Eigen::glmerResp *ans = new lme4Eigen::glmerResp(List(fams), NumericVector(ys));
	return wrap(XPtr<lme4Eigen::glmerResp>(ans, true));
	END_RCPP;
    }
    
    SEXP glmerRespLaplace(SEXP ptr_, SEXP ldL2, SEXP ldRX2, SEXP sqrL) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lme4Eigen::glmerResp>(ptr_)->Laplace(::Rf_asReal(ldL2),
									 ::Rf_asReal(ldRX2),
									 ::Rf_asReal(sqrL)));
	END_RCPP;
    }

    SEXP glmerRespupdateMu(SEXP ptr_, SEXP gamma) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lme4Eigen::glmerResp>(ptr_)->updateMu(as<Eigen::VectorXd>(gamma)));
	END_RCPP;
    }

    SEXP glmerRespdevResid(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::glmerResp>(ptr_)->devResid());
	END_RCPP;
    }

    SEXP glmerRespeta(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::glmerResp>(ptr_)->eta());
	END_RCPP;
    }

    SEXP glmerRespfamily(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::glmerResp>(ptr_)->family());
	END_RCPP;
    }

    SEXP glmerResplink(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::glmerResp>(ptr_)->link());
	END_RCPP;
    }

    SEXP glmerRespresDev(SEXP ptr_) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<lme4Eigen::glmerResp>(ptr_)->resDev());
	END_RCPP;
    }

    SEXP glmerRespsqrtWrkWt(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::glmerResp>(ptr_)->sqrtWrkWt());
	END_RCPP;
    }

    SEXP glmerRespupdateWts(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::glmerResp>(ptr_)->updateWts());
	END_RCPP;
    }

    SEXP glmerRespvariance(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::glmerResp>(ptr_)->variance());
	END_RCPP;
    }

    SEXP glmerRespwrkResids(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::glmerResp>(ptr_)->wrkResids());
	END_RCPP;
    }

    SEXP glmerRespwrkResp(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<lme4Eigen::glmerResp>(ptr_)->wrkResp());
	END_RCPP;
    }

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
	return ::Rf_ScalarReal(XPtr<lme4Eigen::lmerResp>(ptr_)->updateMu(as<Eigen::VectorXd>(gamma)));
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
