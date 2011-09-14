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
    using Eigen::Map;
    using Eigen::VectorXd;

    using Rcpp::CharacterVector;
    using Rcpp::Environment;
    using Rcpp::Language;
    using Rcpp::NumericVector;

    using glm::glmFamily;

    class lmResp {
    protected:
	double                     d_wrss;
	const NumericVector        d_yR;
	const Map<VectorXd>        d_y;
	VectorXd                   d_weights, d_offset, d_mu, d_sqrtXwt, d_sqrtrwt, d_wtres;
    public:
	lmResp(NumericVector);

	const VectorXd&           mu() const {return d_mu;}
	const VectorXd&       offset() const {return d_offset;}
	const VectorXd&      sqrtXwt() const {return d_sqrtXwt;}
	const VectorXd&      sqrtrwt() const {return d_sqrtrwt;}
	const VectorXd&      weights() const {return d_weights;}
	const VectorXd&        wtres() const {return d_wtres;}
	const Map<VectorXd>&       y() const {return d_y;}
	double                  wrss() const {return d_wrss;}
	double              updateMu(const VectorXd&);
	double             updateWts()       {return updateWrss();}
	double            updateWrss();
	void               setOffset(const VectorXd&);
	void              setWeights(const VectorXd&);
    };

    class lmerResp : public lmResp {
    private:
	int d_reml;
    public:
	lmerResp(NumericVector);

	double               Laplace(double,double,double)const;
	int                     REML() const {return d_reml;}
	void                 setReml(int);
    };
	
    class glmResp : public lmResp {
    protected:
	glmFamily       d_fam;
	VectorXd        d_eta, d_n;
    public:
	glmResp(Rcpp::List, NumericVector);

	VectorXd            devResid() const;
	VectorXd               muEta() const;
        VectorXd           sqrtWrkWt() const;
	VectorXd            variance() const;
	VectorXd           wrkResids() const;
	VectorXd             wrkResp() const;

	const VectorXd&          eta() const {return d_eta;}
	const VectorXd&            n() const {return d_n;}

	const std::string&    family() const {return d_fam.fam();}
	const std::string&      link() const {return d_fam.lnk();}

	double               Laplace(double,double,double) const;
	double                resDev() const;
	double              updateMu(const VectorXd&);
	double             updateWts();

	void                    setN(const VectorXd&);
    };

    class nlmerResp : public lmResp {
    protected:
	Environment     d_nlenv;
	Language        d_nlmod;
	CharacterVector d_pnames;
    public:
	nlmerResp(NumericVector,Language,Environment,CharacterVector);

	double            Laplace(double, double, double) const;
	double           updateMu(const VectorXd&);
    };
}

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

    // nonlinear model (and nonlinear mixed model) class

    SEXP nlm_Create(SEXP,SEXP,SEXP,SEXP); // constructor

    SEXP nlm_Laplace(SEXP,SEXP,SEXP,SEXP); // methods
    SEXP nlm_updateMu(SEXP,SEXP);
}

#endif
