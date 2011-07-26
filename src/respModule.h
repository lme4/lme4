// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// respModule.h: response modules using Eigen
//
// Copyright (C)       2011 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#ifndef LME4_RESPMODULE_H
#define LME4_RESPMODULE_H

#include "eigen.h"
#include "glmFamily.h"

namespace lme4Eigen {
    class modResp {
    protected:
	double                     d_wrss;
	const Rcpp::NumericVector  d_yR;
	const MVectorXd            d_y;
	VectorXd                   d_weights, d_offset, d_mu, d_sqrtXwt, d_sqrtrwt, d_wtres;
    public:
//	modResp(Rcpp::S4)                                        throw (std::invalid_argument);
	modResp(Rcpp::NumericVector)                             throw (std::invalid_argument);
//	modResp(Rcpp::NumericVector, Rcpp::NumericVector)        throw (std::invalid_argument);
//	modResp(Rcpp::NumericVector, Rcpp::NumericVector,
//		Rcpp::NumericVector)                             throw (std::invalid_argument);

	const VectorXd&           mu() const {return d_mu;}
	const VectorXd&       offset() const {return d_offset;}
	const VectorXd&      sqrtXwt() const {return d_sqrtXwt;}
	const VectorXd&      sqrtrwt() const {return d_sqrtrwt;}
	const VectorXd&      weights() const {return d_weights;}
	const VectorXd&        wtres() const {return d_wtres;}
	const MVectorXd&           y() const {return d_y;}
	double                  wrss() const {return d_wrss;}
//	double              updateMu(const VectorXd&);
	double             updateWts()       {return updateWrss();}
	double            updateWrss();
//	void                    init();
	void               setOffset(const Rcpp::NumericVector&) throw (std::invalid_argument);
	void              setWeights(const Rcpp::NumericVector&) throw (std::invalid_argument);
    };

    class lmerResp : public modResp {
    private:
	int d_reml;
    public:
//	lmerResp(Rcpp::S4)                                       throw (std::invalid_argument);
	lmerResp(Rcpp::NumericVector)                       throw (std::invalid_argument);
//	lmerResp(int, Rcpp::NumericVector, Rcpp::NumericVector)  throw (std::invalid_argument);
//	lmerResp(int, Rcpp::NumericVector,
//		 Rcpp::NumericVector, Rcpp::NumericVector)       throw (std::invalid_argument);

	// These methods are repeated in subclasses of modResp because of current
	// restrictions on Rcpp module discovery
//	const VectorXd&           mu() const {return d_mu;}
//	const VectorXd&       offset() const {return d_offset;}
//	const VectorXd&      sqrtXwt() const {return d_sqrtXwt;}
//	const VectorXd&      sqrtrwt() const {return d_sqrtrwt;}
//	const VectorXd&      weights() const {return d_weights;}
//	const VectorXd&        wtres() const {return d_wtres;}
//	const MVectorXd&           y() const {return d_y;}
//	double                  wrss() const {return d_wrss;}
//	double             updateWts()       {return updateWrss();}
	double               Laplace(double,double,double)const;
	double              updateMu(const VectorXd&);
	int                     REML() const {return d_reml;}
//	void               setOffset(const VectorXd&)            throw (std::invalid_argument);
	void                 setReml(int)                        throw (std::invalid_argument);
//	void              setWeights(const VectorXd&)            throw (std::invalid_argument);
    };
	
    class glmerResp : public modResp {
    protected:
	glm::glmFamily d_fam;
	VectorXd       d_eta, d_n;
	double         d_pwrss;
    public:
//	glmerResp(Rcpp::S4)                               throw (std::invalid_argument);
	glmerResp(Rcpp::List, Rcpp::NumericVector)        throw (std::invalid_argument);
//	glmerResp(Rcpp::List, Rcpp::NumericVector,
//		  Rcpp::NumericVector)                    throw (std::invalid_argument);
//	glmerResp(Rcpp::List, Rcpp::NumericVector,
//		  Rcpp::NumericVector,
//		  Rcpp::NumericVector)                    throw (std::invalid_argument);
	// glmerResp(Rcpp::List, Rcpp::NumericVector,
	// 	  Rcpp::NumericVector,
	// 	  Rcpp::NumericVector,
	// 	  Rcpp::NumericVector n)                  throw (std::invalid_argument);
	// glmerResp(Rcpp::List, Rcpp::NumericVector,
	// 	  Rcpp::NumericVector,
	// 	  Rcpp::NumericVector,
	// 	  Rcpp::NumericVector,
	// 	  Rcpp::NumericVector)                    throw (std::invalid_argument);

	// Some extractor methods are repeated in subclasses of
	// modResp because of current restrictions on Rcpp module
	// discovery.
	VectorXd            devResid() const {return d_fam.devResid(d_mu, d_weights, d_y);}
	VectorXd               muEta() const {return d_fam.muEta(d_eta);}
        VectorXd           sqrtWrkWt() const;
	VectorXd            variance() const {return d_fam.variance(d_mu);}
        VectorXd           wrkResids() const {return ((d_y - d_mu).array()/muEta().array()).matrix();}
        VectorXd             wrkResp() const {return (d_eta - d_offset) + wrkResids();}
//	const MVectorXd&           y() const {return d_y;}
	const VectorXd&          eta() const {return d_eta;}
//	const VectorXd&           mu() const {return d_mu;}
//	const VectorXd&       offset() const {return d_offset;}
//	const VectorXd&      sqrtXwt() const {return d_sqrtXwt;}
//	const VectorXd&      sqrtrwt() const {return d_sqrtrwt;}
//	const VectorXd&      weights() const {return d_weights;}
//	const VectorXd&        wtres() const {return d_wtres;}
	const std::string&    family() const {return d_fam.fam();}
	const std::string&      link() const {return d_fam.lnk();}
	double               Laplace(double,double,double)const;
	double                 pwrss() const {return d_pwrss;}
	double                resDev() const {return devResid().sum();}
//	double                  wrss() const {return d_wrss;}
	double              updateMu(const VectorXd&);
	double             updateWts();
//	void               setOffset(const VectorXd&)     throw (std::invalid_argument);
	void                setPwrss(double val) {d_pwrss = val;}
//	void              setWeights(const VectorXd&)     throw (std::invalid_argument);
    };

    class nlmerResp : public modResp {
    protected:
	Rcpp::Environment     d_nlenv;
	Rcpp::Language        d_nlmod;
	Rcpp::CharacterVector d_pnames;
    public:
	nlmerResp(Rcpp::S4)                               throw (std::invalid_argument);

	const VectorXd&        mu() const {return d_mu;}
	const VectorXd&    offset() const {return d_offset;}
	const VectorXd&   sqrtXwt() const {return d_sqrtXwt;}
	const VectorXd&   sqrtrwt() const {return d_sqrtrwt;}
	const VectorXd&   weights() const {return d_weights;}
	const VectorXd&     wtres() const {return d_wtres;}
	const MVectorXd&        y() const {return d_y;}
	double            Laplace(double, double, double) const;
	double           updateMu(const VectorXd&) throw(std::invalid_argument);
    };
}

extern "C" {
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
