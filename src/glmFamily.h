// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// glmFamily.h: glm family class using Eigen
//
// Copyright (C)       2012 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.
#ifndef LME4_GLMFAMILY_H
#define LME4_GLMFAMILY_H

#include <RcppEigen.h>

namespace glm {
    using Eigen::ArrayXd;

#if 0
    class glmDist {
    public:
	virtual const ArrayXd variance(const ArrayXd&) const = 0;
	virtual const ArrayXd devResid(const ArrayXd&, const ArrayXd&, const ArrayXd&) const;
	virtual double             aic(const ArrayXd&, const ArrayXd&, const ArrayXd&,
				       const ArrayXd&, double) const;
	/**< in keeping with the botched up nomenclature in the R glm function, 
	 *   the value of aic is the deviance */
    };

    class glmLink {
    public:
	virtual const ArrayXd linkFun(const ArrayXd&) const = 0;
	virtual const ArrayXd linkInv(const ArrayXd&) const;
	virtual const ArrayXd   muEta(const ArrayXd&) const;
    };

    class logLink : public glmLink {
    public:
	const ArrayXd linkFun(const ArrayXd&) const;
	const ArrayXd linkInv(const ArrayXd&) const;
	const ArrayXd   muEta(const ArrayXd&) const;
    };

    class logitLink : public glmLink {
    public:
	const ArrayXd linkFun(const ArrayXd&) const;
	const ArrayXd linkInv(const ArrayXd&) const;
	const ArrayXd   muEta(const ArrayXd&) const;
    };


    class cloglogLink : public glmLink {
    public:
	const ArrayXd linkFun(const ArrayXd&) const;
	const ArrayXd linkInv(const ArrayXd&) const;
	const ArrayXd   muEta(const ArrayXd&) const;
    };

    class probitLink : public glmLink {
    public:
	const ArrayXd linkFun(const ArrayXd&) const;
	const ArrayXd linkInv(const ArrayXd&) const;
	const ArrayXd   muEta(const ArrayXd&) const;
    };

    class identityLink : public glmLink {
    public:
	const ArrayXd linkFun(const ArrayXd&) const;
	const ArrayXd linkInv(const ArrayXd&) const;
	const ArrayXd   muEta(const ArrayXd&) const;
    };

    class  inverseLink : public glmLink {
    public:
	const ArrayXd linkFun(const ArrayXd&) const;
	const ArrayXd linkInv(const ArrayXd&) const;
	const ArrayXd   muEta(const ArrayXd&) const;
    };
#endif    
    /** associative arrays of vector-valued functions */
    typedef std::map<std::string, ArrayXd(*)(const ArrayXd&)> fmap;
    typedef std::map<std::string, ArrayXd(*)(const ArrayXd&,const ArrayXd&,const ArrayXd&)> drmap;
    typedef std::map<std::string, double(*)(const ArrayXd&,const ArrayXd&,const ArrayXd&,
					    const ArrayXd&,double)> aicmap;

    class glmFamily {
    protected:
	std::string    d_family, d_link; /**< as in the R glm family */
				//@{ R functions from the family, as a fall-back
	Rcpp::Function d_devRes, d_linkfun, d_linkinv, d_muEta, d_variance, d_aic;
				//@}
	Rcpp::Environment d_rho;
    public:
	glmFamily(Rcpp::List);

	const std::string& fam() const {return d_family;}
	const std::string& lnk() const {return d_link;}

	// initialize the associative arrays of scalar functions
	void initMaps();
	
	//@{ Application of functions from the family using compiled code when available
	ArrayXd  linkFun(const ArrayXd&) const;
	ArrayXd  linkInv(const ArrayXd&) const;
	virtual ArrayXd devResid(const ArrayXd&, const ArrayXd&, const ArrayXd&) const;
	ArrayXd    muEta(const ArrayXd&) const;
	virtual ArrayXd variance(const ArrayXd&) const;
	virtual double        aic(const ArrayXd&, const ArrayXd&, const ArrayXd&,
				  const ArrayXd&, double) const;
	/**< in keeping with the botched up nomenclature in the R glm function, 
	 *   the value of aic is the deviance */
	//@}
    private:
	// Class members that are maps storing the scalar functions
	static fmap
	    lnks,	        /**< scalar link function */
	    linvs,		/**< scalar linkinv functions */
	    muEtas,		/**< scalar muEta functions */
	    varFuncs;		/**< scalar variance functions */
	static drmap devRes;	/**< scalar deviance residuals functions */
	static aicmap aics;	/**< scalar aic functions */
	
    };

    class negativeBinomial : public glmFamily {
    protected:
	double  d_theta;
    public:
	negativeBinomial(Rcpp::List);
	ArrayXd variance(const ArrayXd&) const;
	ArrayXd devResid(const ArrayXd&, const ArrayXd&, const ArrayXd&) const;
	double        aic(const ArrayXd&, const ArrayXd&, const ArrayXd&,
				  const ArrayXd&, double) const;
	double      theta() const {return d_theta;}
	void     setTheta(const double& ntheta) {d_theta = ntheta;}
    };
}
    
#endif /* LME4_GLMFAMILY_H */

