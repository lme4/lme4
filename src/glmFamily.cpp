//
// glmFamily.cpp: implementation of glmFamily and related classes using Eigen
//
// Copyright (C) 2011-2012 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#include "glmFamily.h"
#include <Rmath.h>
#include <limits>
#include <cmath>

using namespace Rcpp;

namespace glm {
    /** Cumulative probability function of the complement of the Gumbel distribution
     * 
     * (i.e. pgumbel(q,0.,1.,0) == 1 - pgumbel2(-q,0.,1.,0))
     * 
     * @param q the quantile at which to evaluate the cumulative probability
     * @param loc location parameter
     * @param scale scale parameter
     * @param lower_tail when zero evaluate the complement of the cdf
     * 
     * @return Cumulative probability value or its complement, according to the value of lower_tail
     */
    static inline double	
    pgumbel2(const double& q, const double& loc, const double& scale, int lower_tail) {
	double qq = (q - loc) / scale;
	qq = -std::exp(qq);
	return lower_tail ? -expm1(qq) : std::exp(qq);
    }

    /** 
     * density of the complement of the Gumbel distribution
     * 
     * @param x numeric argument
     * @param loc location parameter
     * @param scale scale parameter
     * @param give_log should the logarithm of the density be returned
     * 
     * @return density or its logarithm, according to the value of give_log
     */    
    static inline double
    dgumbel2(const double& x, const double& loc, const double& scale, int give_log) {
	double xx = (x - loc) / scale;
	xx = xx - std::exp(xx) - std::log(scale);
	return give_log ? xx : std::exp(xx);
    }

    //@{ Templated scalar functors used in links, inverse links, etc.
    template <typename T>
    struct logN0 : public std::unary_function<T, T> {
	const T operator()(const T& x) const {return x ? std::log(x) : T();}
    };

    static inline ArrayXd Y_log_Y(const ArrayXd& y, const ArrayXd& mu) {
	return y * (y/mu).unaryExpr(logN0<double>());
    }

    template<typename T>
    struct Round : public std::unary_function<T, T> {
	const T operator()(const T& x) const {return nearbyint(x);}
    };

    template<typename T>
    struct x1mx : public std::unary_function<T, T> {
	const T operator() (const T& x) const {
	    return T(std::max(std::numeric_limits<T>::epsilon(), x * (1 - x)));
	}
    };

    template<typename T>
    struct Lgamma : public std::unary_function<T, T> {
	const T operator() (const T& x) const {
	    return lgamma(x);
	}
    };

    template<typename T>
    struct cauchitinv : public std::unary_function<T, T> {
	const T operator() (const T& x) const {
	    return T(::Rf_pcauchy(double(x), 0., 1., 1, 0));
	}
    };

    template<typename T>
    struct cauchit : public std::unary_function<T, T> {
	const T operator() (const T& x) const {
	    return T(::Rf_qcauchy(double(x), 0., 1., 1, 0));
	}
    };

    template<typename T>
    struct cauchitmueta : public std::unary_function<T, T> {
	const T operator() (const T& x) const {
	    return T(::Rf_dcauchy(double(x), 0., 1., 0));
	}
    };

    template<typename T>
    struct logitinv : public std::unary_function<T, T> {
	const T operator() (const T& x) const {
	    return T(::Rf_plogis(double(x), 0., 1., 1, 0));
	}
    };

    template<typename T>
    struct logit : public std::unary_function<T, T> {
	const T operator() (const T& x) const {
	    return T(::Rf_qlogis(double(x), 0., 1., 1, 0));
	}
    };

    template<typename T>
    struct logitmueta : public std::unary_function<T, T> {
	const T operator() (const T& x) const {
	    return T(::Rf_dlogis(double(x), 0., 1., 0));
	}
    };

    template<typename T>
    struct probitinv : public std::unary_function<T, T> {
	const T operator() (const T& x) const {
	    return T(::Rf_pnorm5(double(x), 0., 1., 1, 0));
	}
    };

    template<typename T>
    struct probit : public std::unary_function<T, T> {
	const T operator() (const T& x) const {
	    return T(::Rf_qnorm5(double(x), 0., 1., 1, 0));
	}
    };

    template<typename T>
    struct probitmueta : public std::unary_function<T, T> {
	const T operator() (const T& x) const {
	    return T(::Rf_dnorm4(double(x), 0., 1., 0));
	}
    };

    template<typename T>
    struct clogloginv : public std::unary_function<T, T> {
	const T operator() (const T& x) const {
	    return T(pgumbel2(double(x), 0., 1., 1));
	}
    };

    template<typename T>
    struct cloglogmueta : public std::unary_function<T, T> {
	const T operator() (const T& x) const {
	    return T(dgumbel2(double(x), 0., 1., 0));
	}
    };
    //@}

    //@{
    double                binomialDist::aic     (const ArrayXd& y, const ArrayXd& n, const ArrayXd& mu,
						 const ArrayXd& wt, double dev) const {
	ArrayXd    m((n > 1).any() ? n : wt);
	ArrayXd   yy((m * y).unaryExpr(Round<double>()));
	m = m.unaryExpr(Round<double>());
	double ans(0.);
	for (int i=0; i < mu.size(); ++i)
	    ans += (m[i] <= 0. ? 0. : wt[i]/m[i]) * ::Rf_dbinom(yy[i], m[i], mu[i], true);
	return (-2. * ans);
    }
    const ArrayXd         binomialDist::devResid(const ArrayXd& y, const ArrayXd& mu, const ArrayXd& wt) const {
	return 2. * wt * (Y_log_Y(y, mu) + Y_log_Y(1. - y, 1. - mu));
    }
    const ArrayXd         binomialDist::variance(const ArrayXd& mu) const {return mu.unaryExpr(x1mx<double>());}
    //@}

    //@{
    double                   gammaDist::aic     (const ArrayXd& y, const ArrayXd& n, const ArrayXd& mu,
						 const ArrayXd& wt, double dev) const {
	double   nn(wt.sum());
	double disp(dev/nn);
	double   ans(0), invdisp(1./disp);
	for (int i = 0; i < mu.size(); ++i)
	    ans += wt[i] * ::Rf_dgamma(y[i], invdisp, mu[i] * disp, true);
	return -2. * ans + 2.;
    }
    const ArrayXd            gammaDist::devResid(const ArrayXd& y, const ArrayXd& mu, const ArrayXd& wt) const {
	return -2. * wt * ((y/mu).unaryExpr(logN0<double>()) - (y - mu)/mu);
    }
    const ArrayXd            gammaDist::variance(const ArrayXd& mu) const {return mu.square();}
    //@}

    //@{
    double                GaussianDist::aic     (const ArrayXd& y, const ArrayXd& n, const ArrayXd& mu,
						 const ArrayXd& wt, double dev) const {
	double   nn(mu.size());
	return nn * (std::log(2. * M_PI * dev/nn) + 1.) + 2. - wt.log().sum();
    }
    const ArrayXd         GaussianDist::devResid(const ArrayXd& y, const ArrayXd& mu, const ArrayXd& wt) const {
	return wt * (y - mu).square();
    }
    const ArrayXd         GaussianDist::variance(const ArrayXd& mu) const {return ArrayXd::Ones(mu.size());}
    //@}

    //@{
    double         inverseGaussianDist::aic     (const ArrayXd& y, const ArrayXd& n, const ArrayXd& mu,
	const ArrayXd& wt, double dev) const {
	double wtsum(wt.sum());
	return wtsum * (std::log(dev/wtsum * 2. * M_PI) + 1.) + 3. * (y.log() * wt).sum() + 2.;
    }
    const ArrayXd  inverseGaussianDist::devResid(const ArrayXd& y, const ArrayXd& mu, const ArrayXd& wt) const {
	return wt * ((y - mu).square())/(y * mu.square());
    }
    const ArrayXd  inverseGaussianDist::variance(const ArrayXd& mu) const {return mu.cube();}
    //@}

    //@{
    double        negativeBinomialDist::aic     (const ArrayXd& y, const ArrayXd& n, const ArrayXd& mu,
						 const ArrayXd& wt, double dev) const {
	return 2. * (wt * (y + d_theta) * (mu + d_theta).log() -
		     y * mu.log() + (y + 1).unaryExpr(Lgamma<double>()) -
		     d_theta * std::log(d_theta) + lgamma(d_theta) -
		     (d_theta + y).unaryExpr(Lgamma<double>())).sum();
    }
    const ArrayXd negativeBinomialDist::devResid(const ArrayXd &y, const ArrayXd &mu, const ArrayXd &wt) const {
	return 2. * wt * (Y_log_Y(y, mu) - (y + d_theta) * ((y + d_theta)/(mu + d_theta)).log());
    }
    const ArrayXd negativeBinomialDist::variance(const ArrayXd &mu) const {
	return mu + mu.square()/d_theta;
    }
    //@}

    //@{
    double                 PoissonDist::aic     (const ArrayXd& y, const ArrayXd& n, const ArrayXd& mu,
						 const ArrayXd& wt, double dev) const {
	double ans(0.);
	for (int i = 0; i < mu.size(); ++i) ans += ::Rf_dpois(y[i], mu[i], true) * wt[i];
	return (-2. * ans);
    }
    const ArrayXd          PoissonDist::devResid(const ArrayXd& y, const ArrayXd& mu, const ArrayXd& wt) const {
	return 2. * wt * (y * (y/mu).unaryExpr(logN0<double>()) - (y - mu));
    }
    const ArrayXd          PoissonDist::variance(const ArrayXd& mu) const {return mu;}
    //@}

    //@{
    const ArrayXd  cauchitLink::linkFun(const ArrayXd&  mu) const {return  mu.unaryExpr(cauchit<double>());}
    const ArrayXd  cauchitLink::linkInv(const ArrayXd& eta) const {return eta.unaryExpr(cauchitinv<double>());}
    const ArrayXd  cauchitLink::muEta(  const ArrayXd& eta) const {return eta.unaryExpr(cauchitmueta<double>());}
    //@}

    //@{
    const ArrayXd      logLink::linkFun(const ArrayXd&  mu) const {return  mu.log();}
    const ArrayXd      logLink::linkInv(const ArrayXd& eta) const {return eta.exp();}
    const ArrayXd      logLink::muEta(  const ArrayXd& eta) const {return eta.exp();}
    //@}

    //@{
    const ArrayXd    logitLink::linkFun(const ArrayXd&  mu) const {return  mu.unaryExpr(logit<double>());}
    const ArrayXd    logitLink::linkInv(const ArrayXd& eta) const {return eta.unaryExpr(logitinv<double>());}
    const ArrayXd    logitLink::muEta(  const ArrayXd& eta) const {return eta.unaryExpr(logitmueta<double>());}
    //@}

    //@{
    const ArrayXd   probitLink::linkFun(const ArrayXd&  mu) const {return  mu.unaryExpr(probit<double>());}
    const ArrayXd   probitLink::linkInv(const ArrayXd& eta) const {return eta.unaryExpr(probitinv<double>());}
    const ArrayXd   probitLink::muEta(  const ArrayXd& eta) const {return eta.unaryExpr(probitmueta<double>());}
    //@}

    //@{
    const ArrayXd identityLink::linkFun(const ArrayXd&  mu) const {return  mu;}
    const ArrayXd identityLink::linkInv(const ArrayXd& eta) const {return eta;}
    const ArrayXd identityLink::muEta(  const ArrayXd& eta) const {return ArrayXd::Ones(eta.size());}
    //@}

    //@{
    const ArrayXd  inverseLink::linkFun(const ArrayXd&  mu) const {return  mu.inverse();}
    const ArrayXd  inverseLink::linkInv(const ArrayXd& eta) const {return eta.inverse();}
    const ArrayXd  inverseLink::muEta(  const ArrayXd& eta) const {return -(eta.inverse().square());}
    //@}

    //@{
//    const ArrayXd  cloglogLink::linkFun(const ArrayXd&  mu) const {return  mu.unaryExpr(cloglog<double>());}
    const ArrayXd  cloglogLink::linkInv(const ArrayXd& eta) const {return eta.unaryExpr(clogloginv<double>());}
    const ArrayXd  cloglogLink::muEta(  const ArrayXd& eta) const {return eta.unaryExpr(cloglogmueta<double>());}
    //@}

    glmDist::glmDist(Rcpp::List& ll)
	: d_devRes  (as<SEXP>(ll["dev.resids"])),
	  d_variance(as<SEXP>(ll["variance"])),
	  d_aic(     as<SEXP>(ll["aic"])),
	  d_rho(     d_aic.environment()) {
    }

    glmLink::glmLink(Rcpp::List& ll)
	: d_linkFun(as<SEXP>(ll["linkfun"])),
	  d_linkInv(as<SEXP>(ll["linkinv"])),
	  d_muEta(  as<SEXP>(ll["mu.eta"])),
	  d_rho(    d_linkFun.environment()) {
    }

    glmFamily::glmFamily(Rcpp::List ll)
	: d_family( as<std::string>(as<SEXP>(ll["family"]))),
	  d_linknam(as<std::string>(as<SEXP>(ll["link"]))),
	  d_dist(   new glmDist(ll)),
	  d_link(   new glmLink(ll)) {
	if (!ll.inherits("family"))
	    throw std::runtime_error("glmFamily requires a list of (S3) class \"family\"");

	if (d_linknam == "cauchit")  {delete d_link; d_link = new cauchitLink(ll);}
	if (d_linknam == "cloglog")  {delete d_link; d_link = new cloglogLink(ll);}
	if (d_linknam == "identity") {delete d_link; d_link = new identityLink(ll);}
	if (d_linknam == "inverse")  {delete d_link; d_link = new inverseLink(ll);}
	if (d_linknam == "log")      {delete d_link; d_link = new logLink(ll);}
	if (d_linknam == "logit")    {delete d_link; d_link = new logitLink(ll);}
	if (d_linknam == "probit")   {delete d_link; d_link = new probitLink(ll);}

	if (d_family  == "binomial")         {delete d_dist; d_dist = new binomialDist(ll);}
	if (d_family  == "Gamma")            {delete d_dist; d_dist = new gammaDist(ll);}
	if (d_family  == "gaussian")         {delete d_dist; d_dist = new GaussianDist(ll);}
	if (d_family  == "inverse.gaussian") {delete d_dist; d_dist = new inverseGaussianDist(ll);}
	if (d_family.substr(0, 18) ==
	    "Negative Binomial(")             {delete d_dist; d_dist = new negativeBinomialDist(ll);}
	if (d_family  == "poisson")          {delete d_dist; d_dist = new PoissonDist(ll);}
    }

    glmFamily::~glmFamily() {
	delete d_dist;
	delete d_link;
    }

    const ArrayXd glmFamily::devResid(const ArrayXd& y, const ArrayXd& mu, const ArrayXd& wt) const {
	return d_dist->devResid(y, mu, wt);
    }

    double glmFamily::aic(const ArrayXd& y, const ArrayXd& n, const ArrayXd& mu,
			  const ArrayXd& wt, double dev) const {
	return d_dist->aic(y, n, mu, wt, dev);
    }

    const ArrayXd glmLink::linkFun(const ArrayXd& mu) const {
	return as<ArrayXd>(::Rf_eval(::Rf_lang2(as<SEXP>(d_linkFun),
						as<SEXP>(Rcpp::NumericVector(mu.data(),
									     mu.data() + mu.size()))
					 ), d_rho));
    }

    const ArrayXd glmLink::linkInv(const ArrayXd& eta) const {
	return as<ArrayXd>(::Rf_eval(::Rf_lang2(as<SEXP>(d_linkInv),
						as<SEXP>(Rcpp::NumericVector(eta.data(),
									     eta.data() + eta.size()))
					 ), d_rho));
    }

    const ArrayXd glmLink::muEta(const ArrayXd &eta) const {
	return as<ArrayXd>(::Rf_eval(::Rf_lang2(as<SEXP>(d_muEta),
						as<SEXP>(Rcpp::NumericVector(eta.data(),
									     eta.data() + eta.size()))
					 ), d_rho));
    }
    
    const ArrayXd glmDist::variance(const ArrayXd &mu) const {
	return as<ArrayXd>(::Rf_eval(::Rf_lang2(as<SEXP>(d_variance),
						as<SEXP>(Rcpp::NumericVector(mu.data(),
									     mu.data() + mu.size()))
					 ), d_rho));
    }
    
    const ArrayXd glmDist::devResid(const ArrayXd &y, const ArrayXd &mu, const ArrayXd &wt) const {
	int n = mu.size();
	return as<ArrayXd>(::Rf_eval(::Rf_lang4(as<SEXP>(d_devRes),
						as<SEXP>(NumericVector(y.data(), y.data() + n)),
						as<SEXP>(NumericVector(mu.data(), mu.data() + n)),
						as<SEXP>(NumericVector(wt.data(), wt.data() + n))
					 ), d_rho));
    }

    double glmDist::aic(const ArrayXd& y, const ArrayXd& n, const ArrayXd& mu,
			const ArrayXd& wt, double dev) const {
	int nn = mu.size();
	double ans =
	    ::Rf_asReal(::Rf_eval(::Rf_lang6(as<SEXP>(d_aic),
					     as<SEXP>(NumericVector(y.data(), y.data() + nn)),
					     as<SEXP>(NumericVector(n.data(), n.data() + nn)),
					     as<SEXP>(NumericVector(mu.data(), mu.data() + nn)),
					     as<SEXP>(NumericVector(wt.data(), wt.data() + nn)),
					     PROTECT(::Rf_ScalarReal(dev))), d_rho));
	UNPROTECT(1);
	return ans;
    }
    
    negativeBinomialDist::negativeBinomialDist(Rcpp::List& ll)
	: glmDist(ll),
	  d_theta(::Rf_asReal(as<SEXP>(d_rho[".Theta"]))) {}


    double glmDist::theta() const {
	throw std::invalid_argument("theta accessor applies only to negative binomial");
    }

    void   glmDist::setTheta(const double& theta) {
	throw std::invalid_argument("setTheta applies only to negative binomial");
    }
}
