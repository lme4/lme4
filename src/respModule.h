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
	const Map<VectorXd>        d_y;
	Map<VectorXd>              d_weights, d_offset, d_mu, d_sqrtXwt, d_sqrtrwt, d_wtres;
    public:
	lmResp(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);

	const Map<VectorXd>& sqrtXwt() const {return d_sqrtXwt;}
	const Map<VectorXd>&      mu() const {return d_mu;}
	const Map<VectorXd>&  offset() const {return d_offset;}
	const Map<VectorXd>& sqrtrwt() const {return d_sqrtrwt;}
	const Map<VectorXd>& weights() const {return d_weights;}
	const Map<VectorXd>&   wtres() const {return d_wtres;}
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
	lmerResp(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);

	double               Laplace(double,double,double)const;
	int                     REML() const {return d_reml;}
	void                 setReml(int);
    };
	
    class glmResp : public lmResp {
    protected:
	glmFamily       d_fam;
	Map<VectorXd>   d_eta, d_n;
    public:
	glmResp(Rcpp::List,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);

	VectorXd            devResid() const;
	VectorXd               muEta() const;
        VectorXd           sqrtWrkWt() const;
	VectorXd            variance() const;
	VectorXd           wrkResids() const;
	VectorXd             wrkResp() const;

	const Map<VectorXd>&     eta() const {return d_eta;}
	const Map<VectorXd>&       n() const {return d_n;}

	const std::string&    family() const {return d_fam.fam();}
	const std::string&      link() const {return d_fam.lnk();}

	double               Laplace(double,double,double) const;
	double                resDev() const;
	double              updateMu(const VectorXd&);
	double             updateWts();

	void                    setN(const VectorXd&);
    };

    class nlsResp : public lmResp {
    protected:
	Environment     d_nlenv;
	Language        d_nlmod;
	CharacterVector d_pnames;
    public:
	nlsResp(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,Language,Environment,CharacterVector);

	double            Laplace(double, double, double) const;
	double           updateMu(const VectorXd&);
    };
}

#endif
