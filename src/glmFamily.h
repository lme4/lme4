// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_GLMFAMILY_H
#define LME4_GLMFAMILY_H

#include <RcppEigen.h>
#include <limits>
namespace glm {
    using Eigen::VectorXd;

    /** associative arrays of functions returning double from a double */
    typedef std::map<std::string, double(*)(const double&)> fmap;
    typedef std::map<std::string, double(*)(const double&,const double&,const double&)> drmap;
    typedef std::map<std::string, double(*)(const VectorXd&,const VectorXd&,const VectorXd&,
					    const VectorXd&,double)> aicmap;

    class glmFamily {
    public:

    typedef std::map<std::string, double(*)(const double&)> fmap; /**< associative array of functions returning double from a double */
    typedef std::map<std::string, double(*)(const double&,const double&,const double&)> drmap; /**< associative array of deviance residual functions */
    typedef std::map<std::string, double(*)(const VectorXd&,const VectorXd&,const VectorXd&,
					    const VectorXd&,double)> aicmap;
    protected:
	Rcpp::List     lst;		 /**< original list from R */
	std::string    d_family, d_link; /**< as in the R glm family */
				//@{ R functions from the family, as a fall-back
	Rcpp::Function d_devRes, d_linkfun, d_linkinv, d_muEta, d_variance, d_aic;
				//@}
    public:
	glmFamily(Rcpp::List);

	const std::string& fam() const {return d_family;}
	const std::string& lnk() const {return d_link;}

	// initialize the associative arrays of scalar functions
	void initMaps();
	
	// Application of functions from the family
	// The scalar transformations use compiled code when available 
	VectorXd  linkFun(const VectorXd&) const;
	VectorXd  linkInv(const VectorXd&) const;
	VectorXd devResid(const VectorXd&, const VectorXd&, const VectorXd&) const;
	VectorXd    muEta(const VectorXd&) const;
	VectorXd variance(const VectorXd&) const;
	double        aic(const VectorXd&, const VectorXd&, const VectorXd&,
			  const VectorXd&, double) const;
	/**< in keeping with the botched up nomenclature in the R glm function, 
	 *   the value of aic is the deviance */
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

    static double epsilon(std::numeric_limits<double>::epsilon());
	/**< Threshold for some comparisons */
}
    
#endif /* LME4_GLMFAMILY_H */

