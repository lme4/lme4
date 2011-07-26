// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_GLMFAMILY_H
#define LME4_GLMFAMILY_H

#include "eigen.h"
#include <Rmath.h>

namespace glm {
    /** associative arrays of functions returning double from a double */
    typedef std::map<std::string, double(*)(const double&)> fmap;
    typedef std::map<std::string, double(*)(const double&,const double&,const double&)> drmap;

    class glmFamily {
    protected:
	Rcpp::List     lst;		 /**< original list from R */
	std::string    d_family, d_link; /**< as in the R glm family */
				//@{ R functions from the family, as a fall-back
	Rcpp::Function d_devRes, d_linkfun, d_linkinv, d_muEta, d_variance;
				//@}
    public:
	glmFamily(Rcpp::List) throw (std::runtime_error);

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
    private:
	// Class members that are maps storing the scalar functions
	static fmap
	    lnks,	        /**< scalar link function */
	    linvs,		/**< scalar linkinv functions */
	    muEtas,		/**< scalar muEta functions */
	    varFuncs;		/**< scalar variance functions */
	static drmap devRes;	/**< scalardeviance residuals functions */
	
	static double epsilon;	/**< Threshold for some comparisons */

	/** 
	 * Utility to evaluate y * log(y/mu) with correct limit at y = 0
	 * 
	 * @param y numerator of log argument, if non-zero
	 * @param mu denominator of log argument
	 * 
	 * @return y * log(y/mu) with correct limit at at y = 0
	 */
	static inline double y_log_y(const double& y, const double& mu) {
	    return y ? y * std::log(y/mu) : 0;
	}


        //@{  Scalar functions for components
	static inline double         cubef(const double& x) {return x * x * x;}
	static inline double          expf(const double& x) {return std::exp(x);}
	static inline double        identf(const double& x) {return x;}
	static inline double     invderivf(const double& x) {return -1./(x * x);}
	static inline double      inversef(const double& x) {return 1./x;}
	static inline double          logf(const double& x) {return std::log(x);}
	static inline double          onef(const double& x) {return 1.;}
	static inline double          sqrf(const double& x) {return x * x;}
	static inline double         sqrtf(const double& x) {return std::sqrt(x);}
	static inline double         twoxf(const double& x) {return 2. * x;}
	static inline double         x1mxf(const double& x) {return std::max(epsilon, x*(1 - x));}
	static inline double  logitLinkInv(const double& x) {return ::Rf_plogis(x, 0., 1., 1, 0);}
	static inline double     logitLink(const double& x) {return ::Rf_qlogis(x, 0., 1., 1, 0);}
	static inline double    logitMuEta(const double& x) {return ::Rf_dlogis(x, 0., 1., 0);}
	static inline double probitLinkInv(const double& x) {return ::Rf_pnorm5(x, 0., 1., 1, 0);}
	static inline double    probitLink(const double& x) {return ::Rf_qnorm5(x, 0., 1., 1, 0);}
	static inline double   probitMuEta(const double& x) {return ::Rf_dnorm4(x, 0., 1., 0);}

	/** cumulative probability function of the complement of the Gumbel distribution
	  * 
          * (i.e. pgumbel(x,0.,1.,0) == 1-pgumbel2(-x,0.,1.,0))
          */
	static inline double	
	pgumbel2(const double& q, const double& loc, const double& scale, int lower_tail) {
	    double qq = (q - loc) / scale;
	    qq = -std::exp(qq);
	    return lower_tail ? -expm1(qq) : std::exp(qq);
	}
	static inline double cloglogLinkInv(const double& x) {
	    return pgumbel2(x, 0., 1., 1);
	}

	//density of the complement of the Gumbel distribution
	static inline double
	dgumbel2(const double& x, const double& loc, const double& scale, int give_log) {
	    double xx = (x - loc) / scale;
	    xx = xx - std::exp(xx) - std::log(scale);
	    return give_log ? xx : std::exp(xx);
	}
	static inline double   cloglogMuEta(const double& x) {
	    return dgumbel2(x, 0., 1., 0);
	}
	
        //@{ deviance residuals functions
	static inline double
	BinomialDevRes(const double& y, const double& mu, const double& wt) {
	    return 2 * wt * (y_log_y(y, mu) + y_log_y(1 - y, 1 - mu));
	}
	static inline double logYMu(const double& y, const double& mu) {
	    return y ? std::log(y/mu) : 0;
	}
	static inline double
	GammaDevRes   (const double& y, const double& mu, const double& wt) {
	    return -2 * wt * (logYMu(y, mu) - (y - mu)/mu);
	}
	static inline double
	GaussianDevRes(const double& y, const double& mu, const double& wt) {
	    double res = y - mu;
	    return wt * res * res;
	}
	static inline double
	PoissonDevRes (const double& y, const double& mu, const double& wt) {
	    return 2 * wt * (y_log_y(y, mu) - (y - mu));
	}
	//@}
    //@}
    };
}
    
#endif /* LME4_GLMFAMILY_H */

