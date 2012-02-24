#include "glmFamily.h"
#include <Rmath.h>
#include <limits>
#include <cmath>

using namespace Rcpp;
using namespace std;

namespace glm {
    /** 
     * Evaluate y * log(y/mu) with correct limit at y = 0
     * 
     * @param y numerator of log argument, if non-zero
     * @param mu denominator of log argument
     * 
     * @return y * log(y/mu) with correct limit at at y = 0
     */
    static inline double y_log_y(const double& y, const double& mu) {
	return y ? y * std::log(y/mu) : 0;
    }

    /** 
     * Evaluate log(y/mu) returning 0 when y = 0
     * 
     * @param y numerator of log argument, if non-zero
     * @param mu denominator of log argument
     * 
     * @return log(y/mu) returning 0 when y = 0
     */
    static inline double logYMu(const double& y, const double& mu) {
	return y ? std::log(y/mu) : 0;
    }
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
    template<typename T>
    struct Round : public std::unary_function<T, T> {
	const T operator()(const T& x) const {return nearbyint(x);}
    };

    template <typename T>
    struct logN0 {
	const T operator()(const T& x) const {return x ? std::log(x) : T();}
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

#if 0
    const ArrayXd      logLink::linkFun(const ArrayXd&  mu) const {return  mu.log();}
    const ArrayXd      logLink::linkInv(const ArrayXd& eta) const {return eta.exp();}
    const ArrayXd      logLink::muEta(  const ArrayXd& eta) const {return eta.exp();}

    const ArrayXd    logitLink::linkFun(const ArrayXd&  mu) const {return  mu.unaryExpr(logit<double>());}
    const ArrayXd    logitLink::linkInv(const ArrayXd& eta) const {return eta.unaryExpr(logitinv<double>());}
    const ArrayXd    logitLink::muEta(  const ArrayXd& eta) const {return eta.unaryExpr(logitmueta<double>());}

    const ArrayXd   probitLink::linkFun(const ArrayXd&  mu) const {return  mu.unaryExpr(probit<double>());}
    const ArrayXd   probitLink::linkInv(const ArrayXd& eta) const {return eta.unaryExpr(probitinv<double>());}
    const ArrayXd   probitLink::muEta(  const ArrayXd& eta) const {return eta.unaryExpr(probitmueta<double>());}

    const ArrayXd identityLink::linkFun(const ArrayXd&  mu) const {return  mu;}
    const ArrayXd identityLink::linkInv(const ArrayXd& eta) const {return eta;}
    const ArrayXd identityLink::muEta(  const ArrayXd& eta) const {return ArrayXd::Ones(eta.size());}

    const ArrayXd  inverseLink::linkFun(const ArrayXd&  mu) const {return  mu.inverse();}
    const ArrayXd  inverseLink::linkInv(const ArrayXd& eta) const {return eta.inverse();}
    const ArrayXd  inverseLink::muEta(  const ArrayXd& eta) const {return -(eta.inverse().square());}

//    const ArrayXd  cloglogLink::linkFun(const ArrayXd&  mu) const {return  mu.unaryExpr(cloglog<double>());}
    const ArrayXd  cloglogLink::linkFun(const ArrayXd&  mu) const {return  mu;}
    const ArrayXd  cloglogLink::linkInv(const ArrayXd& eta) const {return eta.unaryExpr(clogloginv<double>());}
    const ArrayXd  cloglogLink::muEta(  const ArrayXd& eta) const {return eta.unaryExpr(cloglogmueta<double>());}

#endif
    //@{  Vector functions for links, inverse links, derivatives and variances
    static ArrayXd cubef(     const ArrayXd& x) {return x.cube();}
    static ArrayXd expf(      const ArrayXd& x) {return x.exp();}
    static ArrayXd identf(    const ArrayXd& x) {return x;}
    static ArrayXd invderivf( const ArrayXd& x) {return -(x.inverse().square());}
    static ArrayXd inversef(  const ArrayXd& x) {return x.inverse();}
    static ArrayXd logf(      const ArrayXd& x) {return x.log();}
    static ArrayXd onef(      const ArrayXd& x) {return ArrayXd::Ones(x.size());}
    static ArrayXd sqrf(      const ArrayXd& x) {return x.square();}
    static ArrayXd sqrtf(     const ArrayXd& x) {return x.sqrt();}
    static ArrayXd twoxf(     const ArrayXd& x) {return x * 2.;}
    static ArrayXd x1mxf(     const ArrayXd& x) {return x.unaryExpr(x1mx<double>());}
    static ArrayXd logitf(    const ArrayXd& x) {return x.unaryExpr(logit<double>());}
    static ArrayXd logitinvf( const ArrayXd& x) {return x.unaryExpr(logitinv<double>());}
    static ArrayXd logitmuf(  const ArrayXd& x) {return x.unaryExpr(logitmueta<double>());}
    static ArrayXd probitf(   const ArrayXd& x) {return x.unaryExpr(probit<double>());}
    static ArrayXd probitinvf(const ArrayXd& x) {return x.unaryExpr(probitinv<double>());}
    static ArrayXd probitmuf( const ArrayXd& x) {return x.unaryExpr(probitmueta<double>());}
    static ArrayXd clogloginf(const ArrayXd& x) {return x.unaryExpr(clogloginv<double>());}
    static ArrayXd cloglogmuf(const ArrayXd& x) {return x.unaryExpr(cloglogmueta<double>());}
    //@}
    
    //@{ deviance residuals functions
    static inline Eigen::ArrayXd Y_log_Y(const Eigen::ArrayXd& y, const Eigen::ArrayXd& mu) {
	return y * (y/mu).unaryExpr(logN0<double>());
    }

    static inline ArrayXd
    BinomialDevRes(const Eigen::ArrayXd& y, const Eigen::ArrayXd& mu, const Eigen::ArrayXd& wt) {
	return 2. * wt * (Y_log_Y(y, mu) + Y_log_Y(1. - y, 1. - mu));
    }
    static inline ArrayXd
    GammaDevRes   (const Eigen::ArrayXd& y, const Eigen::ArrayXd& mu, const Eigen::ArrayXd& wt) {
	return -2. * wt * ((y/mu).unaryExpr(logN0<double>()) - (y - mu)/mu);
    }
    static inline ArrayXd
    GaussianDevRes(const Eigen::ArrayXd& y, const Eigen::ArrayXd& mu, const Eigen::ArrayXd& wt) {
	return wt * (y - mu).square();
    }
    static inline ArrayXd
    PoissonDevRes (const Eigen::ArrayXd& y, const Eigen::ArrayXd& mu, const Eigen::ArrayXd& wt) {
	return 2. * wt * (y * (y/mu).unaryExpr(logN0<double>()) - (y - mu));
    }
    //@}

    //@{  AIC functions (which actually return the deviance, go figure)
    static inline double
    BinomialAIC   (const ArrayXd& y, const ArrayXd& n, const ArrayXd& mu,
		   const ArrayXd& wt, double dev) {
	Eigen::ArrayXd m((n > 1).any() ? n : wt);
	Eigen::ArrayXd yy((m * y).unaryExpr(Round<double>()));
	m = m.unaryExpr(Round<double>());
	double ans(0.);
	for (int i=0; i < mu.size(); ++i)
	    ans += (m[i] <= 0. ? 0. : (wt[i]/m[i]) * ::Rf_dbinom(yy[i], m[i], mu[i], true));
	return (-2. * ans);
    }

    static inline double
    PoissonAIC    (const ArrayXd& y, const ArrayXd& n, const ArrayXd& mu,
		   const ArrayXd& wt, double dev) {
	double ans(0.);
	for (int i=0; i < mu.size(); ++i) ans += ::Rf_dpois(y[i], mu[i], true) * wt[i];
	return (-2. * ans);
    }
    //@}
    
    // initialize the function maps (i.e. associative arrays of
    // functions).  Needed because these are static maps.
    drmap  glmFamily::devRes   = drmap();
    aicmap glmFamily::aics     = aicmap();

    fmap   glmFamily::linvs    = fmap();
    fmap   glmFamily::lnks     = fmap();
    fmap   glmFamily::muEtas   = fmap();
    fmap   glmFamily::varFuncs = fmap();
    
    /** 
     * Initialize the static maps.  The identity link is guaranteed to be initialized if
     * any of the maps are initialized.  FIXME: Should probably check the size of the map instead?
     * 
     */
    void glmFamily::initMaps() {
	lnks    ["log"]              = &logf;
	linvs   ["log"]              = &expf;
	muEtas  ["log"]              = &expf;
	    
	lnks    ["sqrt"]             = &sqrtf;
	linvs   ["sqrt"]             = &sqrf;
	muEtas  ["sqrt"]             = &twoxf;
	    
	lnks    ["identity"]         = &identf;
	linvs   ["identity"]         = &identf;
	muEtas  ["identity"]         = &onef;
	    
	lnks    ["inverse"]          = &inversef;
	linvs   ["inverse"]          = &inversef;
	muEtas  ["inverse"]          = &invderivf;
	    
	lnks    ["logit"]            = &logitf;
	linvs   ["logit"]            = &logitinvf;
	muEtas  ["logit"]            = &logitmuf;
	    
	lnks    ["probit"]           = &probitf;
	linvs   ["probit"]           = &probitinvf;
	muEtas  ["probit"]           = &probitmuf;
	    
	linvs   ["cloglog"]          = &clogloginf;
	muEtas  ["cloglog"]          = &cloglogmuf;
	    
	devRes  ["Gamma"]            = &GammaDevRes;
	varFuncs["Gamma"]            = &sqrf; // x^2

	aics    ["binomial"]         = &BinomialAIC;
	devRes  ["binomial"]         = &BinomialDevRes;
	varFuncs["binomial"]         = &x1mxf; // x * (1 - x)

	devRes  ["gaussian"]         = &GaussianDevRes;
	varFuncs["gaussian"]         = &onef; // 1

	varFuncs["inverse.gaussian"] = &cubef;  // x^3

	aics    ["poisson"]          = &PoissonAIC;
	devRes  ["poisson"]          = &PoissonDevRes;
	varFuncs["poisson"]          = &identf; // x
    }
    
    glmFamily::glmFamily(List ll)
	: d_family(  as<std::string>(as<SEXP>(ll["family"]))),
	  d_link(    as<std::string>(as<SEXP>(ll["link"]))),
	  d_devRes(  as<SEXP>(ll["dev.resids"])),
	  d_linkfun( as<SEXP>(ll["linkfun"])),
	  d_linkinv( as<SEXP>(ll["linkinv"])),
	  d_muEta(   as<SEXP>(ll["mu.eta"])),
	  d_variance(as<SEXP>(ll["variance"])),
	  d_aic(     as<SEXP>(ll["aic"])),
	  d_rho(     d_devRes.environment()) {
	if (!ll.inherits("family"))
	    throw std::runtime_error("glmFamily requires a list of (S3) class \"family\"");
	if (!lnks.count("identity")) initMaps();
    }

    ArrayXd glmFamily::linkFun(const ArrayXd& mu) const {
	if (lnks.count(d_link)) return lnks[d_link](mu);
	return as<ArrayXd>(d_linkfun(NumericVector(mu.data(), mu.data() + mu.size())));
    }

    ArrayXd glmFamily::linkInv(const ArrayXd& eta) const {
	if (linvs.count(d_link)) return linvs[d_link](eta);
	return as<ArrayXd>(d_linkinv(NumericVector(eta.data(), eta.data() + eta.size())));
    }

    ArrayXd glmFamily::muEta(const ArrayXd &eta) const {
	if (muEtas.count(d_link)) return muEtas[d_link](eta);
	return as<ArrayXd>(as<SEXP>(d_muEta(NumericVector(eta.data(), eta.data() + eta.size()))));
    }
    
    ArrayXd glmFamily::variance(const ArrayXd &mu) const {
	if (varFuncs.count(d_family)) return varFuncs[d_family](mu);
	return as<ArrayXd>(as<SEXP>(d_variance(NumericVector(mu.data(), mu.data() + mu.size()))));
    }
    
/// FIXME: change this so the looping is done in the devResid[d_family] function
    ArrayXd
    glmFamily::devResid(const ArrayXd &mu, const ArrayXd &weights, const ArrayXd &y) const {
	if (devRes.count(d_family)) return devRes[d_family](y, mu, weights);
	int n = mu.size();
	return as<ArrayXd>(as<SEXP>(d_devRes(NumericVector(y.data(), y.data() + n),
					      NumericVector(mu.data(), mu.data() + n),
					      NumericVector(weights.data(), weights.data() + n))));
    }

    double glmFamily::aic(const ArrayXd& y, const ArrayXd& n, const ArrayXd& mu,
			  const ArrayXd& wt, double dev) const {
	int nn = mu.size();
	if (aics.count(d_family)) return aics[d_family](y, n, mu, wt, dev);
	SEXP ans = d_aic(NumericVector(y.data(), y.data() + nn),
			 NumericVector(n.data(), n.data() + nn),
			 NumericVector(mu.data(), mu.data() + nn),
			 NumericVector(wt.data(), wt.data() + nn),
			 ::Rf_ScalarReal(dev));
	return Rcpp::as<double>(ans);
    }
    
    negativeBinomial::negativeBinomial(Rcpp::List ll)
	: glmFamily(ll),
	  d_theta(::Rf_asReal(as<SEXP>(d_rho[".Theta"]))) {}

    ArrayXd negativeBinomial::variance(const ArrayXd &mu) const {
	Rcpp::Rcout << "in negativeBinomial::variance with theta = " << d_theta << std::endl;
	return mu + mu.square()/d_theta;
    }

    ArrayXd negativeBinomial::devResid(const ArrayXd &mu, const ArrayXd &weights,
					const ArrayXd &y) const {
	return 2. * weights * (Y_log_Y(y, mu) - (y + d_theta) *
			       ((y + d_theta)/(mu + d_theta)).log());
    }

    double negativeBinomial::aic(const ArrayXd& y, const ArrayXd& n, const ArrayXd& mu,
				 const ArrayXd& wt, double dev) const {
	return 2. * (wt * (y + d_theta) * (mu + d_theta).log() -
		     y * mu.log() + (y + 1).unaryExpr(Lgamma<double>()) -
		     d_theta * std::log(d_theta) + lgamma(d_theta) -
		     (d_theta + y).unaryExpr(Lgamma<double>())).sum();
    }
}
