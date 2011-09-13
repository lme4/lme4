// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// respModule.h: response modules using Eigen
//
// Copyright (C)       2011 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#ifndef LME4_RESPMODULE_H
#define LME4_RESPMODULE_H

#include "glmFamily.h"

namespace lme4Eigen {
    using Eigen::VectorXd;
    using Eigen::Map;
    using Rcpp::NumericVector;

    class modResp {
    protected:
	double                     d_wrss;
	const NumericVector        d_yR;
	const Map<VectorXd>        d_y;
	VectorXd                   d_weights, d_offset, d_mu, d_sqrtXwt, d_sqrtrwt, d_wtres;
    public:
	modResp(NumericVector);

	const VectorXd&           mu() const {return d_mu;}
	const VectorXd&       offset() const {return d_offset;}
	const VectorXd&      sqrtXwt() const {return d_sqrtXwt;}
	const VectorXd&      sqrtrwt() const {return d_sqrtrwt;}
	const VectorXd&      weights() const {return d_weights;}
	const VectorXd&        wtres() const {return d_wtres;}
	const Map<VectorXd>&       y() const {return d_y;}
	double                  wrss() const {return d_wrss;}
	double             updateWts()       {return updateWrss();}
	double            updateWrss();
	void               setOffset(const NumericVector&);
	void              setWeights(const NumericVector&);
    };

    class lmerResp : public modResp {
    private:
	int d_reml;
    public:
	lmerResp(NumericVector);

	double               Laplace(double,double,double)const;
	double              updateMu(const VectorXd&);
	int                     REML() const {return d_reml;}
	void                 setReml(int);
    };
	
    class glmerResp : public modResp {
    protected:
	glm::glmFamily d_fam;
	VectorXd       d_eta, d_n;
	double         d_pwrss;
    public:
	glmerResp(Rcpp::List, NumericVector);

	VectorXd            devResid() const {return d_fam.devResid(d_mu, d_weights, d_y);}
	VectorXd               muEta() const {return d_fam.muEta(d_eta);}
        VectorXd           sqrtWrkWt() const;
	VectorXd            variance() const {return d_fam.variance(d_mu);}
        VectorXd           wrkResids() const {return (d_y - d_mu).cwiseQuotient(muEta());}
        VectorXd             wrkResp() const {return (d_eta - d_offset) + wrkResids();}
	const VectorXd&          eta() const {return d_eta;}
	const std::string&    family() const {return d_fam.fam();}
	const std::string&      link() const {return d_fam.lnk();}
	double               Laplace(double,double,double) const;
	double                 pwrss() const {return d_pwrss;}
	double                resDev() const {return devResid().sum();}
	double              updateMu(const VectorXd&);
	double             updateWts();
	void                setPwrss(double val) {d_pwrss = val;}
    };

    class nlmerResp : public modResp {
    protected:
	Rcpp::Environment     d_nlenv;
	Rcpp::Language        d_nlmod;
	Rcpp::CharacterVector d_pnames;
    public:
	nlmerResp(Rcpp::S4);

	const VectorXd&        mu() const {return d_mu;}
	const VectorXd&    offset() const {return d_offset;}
	const VectorXd&   sqrtXwt() const {return d_sqrtXwt;}
	const VectorXd&   sqrtrwt() const {return d_sqrtrwt;}
	const VectorXd&   weights() const {return d_weights;}
	const VectorXd&     wtres() const {return d_wtres;}
	const Map<VectorXd>&    y() const {return d_y;}
	double            Laplace(double, double, double) const;
	double           updateMu(const VectorXd&);
    };
}

extern "C" {
    SEXP glmerRespCreate(SEXP,SEXP);

    SEXP glmerRespLaplace(SEXP, SEXP, SEXP, SEXP);
    SEXP glmerRespdevResid(SEXP);
    SEXP glmerRespeta(SEXP);
    SEXP glmerRespfamily(SEXP);
    SEXP glmerResplink(SEXP);
    SEXP glmerRespresDev(SEXP);
    SEXP glmerRespsqrtWrkWt(SEXP);
    SEXP glmerRespsqrtXwt(SEXP);
    SEXP glmerRespupdateMu(SEXP,SEXP);
    SEXP glmerRespupdateWts(SEXP);
    SEXP glmerRespvariance(SEXP);
    SEXP glmerRespwrkResids(SEXP);
    SEXP glmerRespwrkResp(SEXP);

    SEXP lmerRespCreate(SEXP);

    SEXP lmerRespsetREML(SEXP, SEXP);
    SEXP lmerRespREML(SEXP);

    SEXP lmerRespLaplace(SEXP, SEXP, SEXP, SEXP);
    SEXP lmerRespupdateMu(SEXP, SEXP);

    SEXP modRespsetOffset(SEXP, SEXP);
    SEXP modRespsetWeights(SEXP, SEXP);

    SEXP modRespmu(SEXP);
    SEXP modRespoffset(SEXP);
    SEXP modRespweights(SEXP);
    SEXP modRespwrss(SEXP);
    SEXP modRespwtres(SEXP);
    SEXP modRespy(SEXP);
}

#endif
