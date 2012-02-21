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

    //@{  Vector functions for links, inverse links, derivatives and variances
    static VectorXd cubef(     const VectorXd& x) {return x.array().cube();}
    static VectorXd expf(      const VectorXd& x) {return x.array().exp();}
    static VectorXd identf(    const VectorXd& x) {return x;}
    static VectorXd invderivf( const VectorXd& x) {return -(x.array().inverse().square());}
    static VectorXd inversef(  const VectorXd& x) {return x.array().inverse();}
    static VectorXd logf(      const VectorXd& x) {return x.array().log();}
    static VectorXd onef(      const VectorXd& x) {return VectorXd::Ones(x.size());}
    static VectorXd sqrf(      const VectorXd& x) {return x.array().square();}
    static VectorXd sqrtf(     const VectorXd& x) {return x.array().sqrt();}
    static VectorXd twoxf(     const VectorXd& x) {return x * 2.;}
    static VectorXd x1mxf(     const VectorXd& x) {return x.array().unaryExpr(x1mx<double>());}
    static VectorXd logitf(    const VectorXd& x) {return x.array().unaryExpr(logit<double>());}
    static VectorXd logitinvf( const VectorXd& x) {return x.array().unaryExpr(logitinv<double>());}
    static VectorXd logitmuf(  const VectorXd& x) {return x.array().unaryExpr(logitmueta<double>());}
    static VectorXd probitf(   const VectorXd& x) {return x.array().unaryExpr(probit<double>());}
    static VectorXd probitinvf(const VectorXd& x) {return x.array().unaryExpr(probitinv<double>());}
    static VectorXd probitmuf( const VectorXd& x) {return x.array().unaryExpr(probitmueta<double>());}
    static VectorXd clogloginf(const VectorXd& x) {return x.array().unaryExpr(clogloginv<double>());}
    static VectorXd cloglogmuf(const VectorXd& x) {return x.array().unaryExpr(cloglogmueta<double>());}
    //@}
    
    //@{ deviance residuals functions
    static inline Eigen::ArrayXd Y_log_Y(const Eigen::ArrayXd& y, const Eigen::ArrayXd& mu) {
	return y * (y/mu).unaryExpr(logN0<double>());
    }

    static inline VectorXd
    BinomialDevRes(const Eigen::ArrayXd& y, const Eigen::ArrayXd& mu, const Eigen::ArrayXd& wt) {
	return 2. * wt * (Y_log_Y(y, mu) + Y_log_Y(1. - y, 1. - mu));
    }
    static inline VectorXd
    GammaDevRes   (const Eigen::ArrayXd& y, const Eigen::ArrayXd& mu, const Eigen::ArrayXd& wt) {
	return -2. * wt * ((y/mu).unaryExpr(logN0<double>()) - (y - mu)/mu);
    }
    static inline VectorXd
    GaussianDevRes(const Eigen::ArrayXd& y, const Eigen::ArrayXd& mu, const Eigen::ArrayXd& wt) {
	return wt * (y - mu).square();
    }
    static inline VectorXd
    PoissonDevRes (const Eigen::ArrayXd& y, const Eigen::ArrayXd& mu, const Eigen::ArrayXd& wt) {
	return 2. * wt * (y * (y/mu).unaryExpr(logN0<double>()) - (y - mu));
    }
    //@}

    //@{  AIC functions (which actually return the deviance, go figure)
    static inline double
    BinomialAIC   (const VectorXd& y, const VectorXd& n, const VectorXd& mu,
		   const VectorXd& wt, double dev) {
	Eigen::ArrayXd m((n.array() > 1).any() ? n : wt);
	Eigen::ArrayXd yy((m * y.array()).unaryExpr(Round<double>()));
	m = m.unaryExpr(Round<double>());
	double ans(0.);
	for (int i=0; i < mu.size(); ++i)
	    ans += (m[i] <= 0. ? 0. : (wt[i]/m[i]) * ::Rf_dbinom(yy[i], m[i], mu[i], true));
	return (-2. * ans);
    }

    static inline double
    PoissonAIC    (const VectorXd& y, const VectorXd& n, const VectorXd& mu,
		   const VectorXd& wt, double dev) {
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

    VectorXd glmFamily::linkFun(const VectorXd& mu) const {
	if (lnks.count(d_link)) return lnks[d_link](mu);
	return as<VectorXd>(d_linkfun(NumericVector(mu.data(), mu.data() + mu.size())));
    }

    VectorXd glmFamily::linkInv(const VectorXd& eta) const {
	if (linvs.count(d_link)) return linvs[d_link](eta);
	return as<VectorXd>(d_linkinv(NumericVector(eta.data(), eta.data() + eta.size())));
    }

    VectorXd glmFamily::muEta(const VectorXd &eta) const {
	if (muEtas.count(d_link)) return muEtas[d_link](eta);
	return as<VectorXd>(as<SEXP>(d_muEta(NumericVector(eta.data(), eta.data() + eta.size()))));
    }
    
    VectorXd glmFamily::variance(const VectorXd &mu) const {
	if (varFuncs.count(d_family)) return varFuncs[d_family](mu);
	return as<VectorXd>(as<SEXP>(d_variance(NumericVector(mu.data(), mu.data() + mu.size()))));
    }
    
/// FIXME: change this so the looping is done in the devResid[d_family] function
    VectorXd
    glmFamily::devResid(const VectorXd &mu, const VectorXd &weights, const VectorXd &y) const {
	if (devRes.count(d_family)) return devRes[d_family](y, mu, weights);
	int n = mu.size();
	return as<VectorXd>(as<SEXP>(d_devRes(NumericVector(y.data(), y.data() + n),
					      NumericVector(mu.data(), mu.data() + n),
					      NumericVector(weights.data(), weights.data() + n))));
    }

    double glmFamily::aic(const VectorXd& y, const VectorXd& n, const VectorXd& mu,
			  const VectorXd& wt, double dev) const {
	int nn = mu.size();
	if (aics.count(d_family)) return aics[d_family](y, n, mu, wt, dev);
	SEXP ans = d_aic(NumericVector(y.data(), y.data() + nn),
			 NumericVector(n.data(), n.data() + nn),
			 NumericVector(mu.data(), mu.data() + nn),
			 NumericVector(wt.data(), wt.data() + nn),
			 ::Rf_ScalarReal(dev));
	return Rcpp::as<double>(ans);
    }
}
