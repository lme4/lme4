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

namespace lme4 {
    typedef Eigen::Map<Eigen::VectorXd> MVec;

    using Rcpp::CharacterVector;
    using Rcpp::Environment;
    using Rcpp::Language;
    using Rcpp::NumericVector;

    using glm::glmFamily;

    class lmResp {
    protected:
	double d_wrss; /**< current weighted sum of squared residuals */
	MVec d_y,      /**< response vector */
	    d_weights, /**< prior weights - always present even if unity */
	    d_offset,  /**< offset in the model */
	    d_mu,      /**< mean response from current linear predictor */
	    d_sqrtXwt, /**< Square roots of the "X weights".  For
			* lmResp and lmerResp these are the same as
			* the sqrtrwt.  For glmResp and nlsResp they
			* incorporate the gradient of the eta to mu
			* mapping.*/
	    d_sqrtrwt,	  /**< Square roots of the residual weights */
	    d_wtres;	  /**< Current weighted residuals */
    public:
	lmResp(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);

	const MVec&    sqrtXwt() const {return d_sqrtXwt;}
				/**< return a const reference to d_sqrtXwt */
	const MVec&         mu() const {return d_mu;} 
				/**< return a const reference to d_mu */
	const MVec&     offset() const {return d_offset;}
				/**< return a const reference to d_offset */
	const MVec&    sqrtrwt() const {return d_sqrtrwt;}
				/**< return a const reference to d_sqrtrwt */
	const MVec&    weights() const {return d_weights;}
				/**< return a const reference to d_weights */
	const MVec&      wtres() const {return d_wtres;}
				/**< return a const reference to d_wtres */
	const MVec&          y() const {return d_y;}
				/**< return a const reference to d_y */
	double            wrss() const {return d_wrss;}
				/**< return the weighted sum of squared residuals */
	double        updateMu(const Eigen::VectorXd&);
	double       updateWts()       {return updateWrss();}
				/**< update the weights.  For a
				 * glmResp this done separately from
				 * updating the mean, because of the
				 * iterative reweighting. */
	double      updateWrss(); /**< update the weighted residuals and d_wrss */
	void         setOffset(const Eigen::VectorXd&);
				/**< set a new value of the offset */
	void           setResp(const Eigen::VectorXd&);
				/**< set a new value of the response, y */
	void        setWeights(const Eigen::VectorXd&);
				/**< set a new value of the prior weights */
    };

    class lmerResp : public lmResp {
    private:
	int d_reml;		/**< 0 for evaluating the deviance, p
				 * for evaluating the REML criterion. */
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

	Eigen::ArrayXd   devResid() const;
	Eigen::ArrayXd      muEta() const;
        Eigen::ArrayXd  sqrtWrkWt() const;
	Eigen::ArrayXd   variance() const;
	Eigen::ArrayXd  wrkResids() const;
	Eigen::ArrayXd    wrkResp() const;
	Eigen::ArrayXd  wtWrkResp() const;

	const MVec&           eta() const {return d_eta;}
	const MVec&             n() const {return d_n;}

	const std::string& family() const {return d_fam.fam();}
	const std::string&   link() const {return d_fam.lnk();}

	double                aic() const;
	double            Laplace(double,double,double) const;
	double             resDev() const;
	double              theta() const {return d_fam.theta();}
				//< negative binomial distribution only
	double           updateMu(const Eigen::VectorXd&);
	double          updateWts();

	void                 setN(const Eigen::VectorXd&);
	void             setTheta(const double& ntheta) {d_fam.setTheta(ntheta);}
				// negative binomial distribution only
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
