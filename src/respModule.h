// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// respModule.h: response modules using Eigen
//
// Copyright (C) 2011-2012 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#ifndef LME4_RESPMODULE_H
#define LME4_RESPMODULE_H

#include "glmFamily.h"

namespace lme4Eigen {
    typedef Eigen::Map<Eigen::VectorXd> MVec;

    using Rcpp::CharacterVector;
    using Rcpp::Environment;
    using Rcpp::Language;
    using Rcpp::NumericVector;

    using glm::glmFamily;

    class lmResp {
    protected:
	double            d_wrss;
	const MVec        d_y;
	MVec              d_weights, d_offset, d_mu, d_sqrtXwt, d_sqrtrwt, d_wtres;
    public:
	lmResp(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);

	const MVec&    sqrtXwt() const {return d_sqrtXwt;}
	const MVec&         mu() const {return d_mu;}
	const MVec&     offset() const {return d_offset;}
	const MVec&    sqrtrwt() const {return d_sqrtrwt;}
	const MVec&    weights() const {return d_weights;}
	const MVec&      wtres() const {return d_wtres;}
	const MVec&          y() const {return d_y;}
	double            wrss() const {return d_wrss;}
	double        updateMu(const Eigen::VectorXd&);
	double       updateWts()       {return updateWrss();}
	double      updateWrss();
	void         setOffset(const Eigen::VectorXd&);
	void        setWeights(const Eigen::VectorXd&);
    };

    class lmerResp : public lmResp {
    private:
	int d_reml;
    public:
	lmerResp(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);

	double         Laplace(double,double,double)const;
	int               REML() const {return d_reml;}
	void           setReml(int);
    };
	
    class glmResp : public lmResp {
    protected:
	glmFamily  d_fam;
	MVec       d_eta, d_n;
    public:
	glmResp(Rcpp::List,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);

	Eigen::VectorXd  devResid() const;
	Eigen::VectorXd     muEta() const;
        Eigen::VectorXd sqrtWrkWt() const;
	Eigen::VectorXd  variance() const;
	Eigen::VectorXd wrkResids() const;
	Eigen::VectorXd   wrkResp() const;

	const MVec&           eta() const {return d_eta;}
	const MVec&             n() const {return d_n;}

	const std::string& family() const {return d_fam.fam();}
	const std::string&   link() const {return d_fam.lnk();}

	double            Laplace(double,double,double) const;
	double             resDev() const;
	double           updateMu(const Eigen::VectorXd&);
	double          updateWts();

	void                 setN(const Eigen::VectorXd&);
    };

    class nlsResp : public lmResp {
    protected:
	MVec            d_gamma;
	Environment     d_nlenv;
	Language        d_nlmod;
	CharacterVector d_pnames;
    public:
	nlsResp(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);

	double            Laplace(double, double, double) const;
	double           updateMu(const Eigen::VectorXd&);
    };
}

#endif
