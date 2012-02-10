#include "glmFamily.h"
#include <Rmath.h>
#include <cmath>

using namespace Rcpp;
using namespace std;

namespace glm {
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

    /** Cumulative probability function of the complement of the Gumbel distribution
     * 
     * (i.e. pgumbel(q,0.,1.,0) == 1-pgumbel2(-q,0.,1.,0))
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

    template<typename Scalar>
    struct Round {
	const Scalar operator()(const Scalar& x) const {return nearbyint(x);}
    };

    //@{  Scalar functions for components
    static inline double          cubef(const double& x) {return x * x * x;}
    static inline double           expf(const double& x) {return std::exp(x);}
    static inline double         identf(const double& x) {return x;}
    static inline double      invderivf(const double& x) {return -1./(x * x);}
    static inline double       inversef(const double& x) {return 1./x;}
    static inline double           logf(const double& x) {return std::log(x);}
    static inline double           onef(const double& x) {return 1.;}
    static inline double           sqrf(const double& x) {return x * x;}
    static inline double          sqrtf(const double& x) {return std::sqrt(x);}
    static inline double          twoxf(const double& x) {return 2. * x;}
    static inline double          x1mxf(const double& x) {return std::max(epsilon, x*(1 - x));}
    static inline double   logitLinkInv(const double& x) {return ::Rf_plogis(x, 0., 1., 1, 0);}
    static inline double      logitLink(const double& x) {return ::Rf_qlogis(x, 0., 1., 1, 0);}
    static inline double     logitMuEta(const double& x) {return ::Rf_dlogis(x, 0., 1., 0);}
    static inline double  probitLinkInv(const double& x) {return ::Rf_pnorm5(x, 0., 1., 1, 0);}
    static inline double     probitLink(const double& x) {return ::Rf_qnorm5(x, 0., 1., 1, 0);}
    static inline double    probitMuEta(const double& x) {return ::Rf_dnorm4(x, 0., 1., 0);}
    static inline double cloglogLinkInv(const double& x) {return pgumbel2(x, 0., 1., 1);}
    //@}
    
    static inline double   cloglogMuEta(const double& x) {return dgumbel2(x, 0., 1., 0);}
    
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
	    lnks["log"]                  = &logf;
	    muEtas["log"] = linvs["log"] = &expf;
	    
	    lnks["sqrt"]                 = &sqrtf;
	    linvs["sqrt"]                = &sqrf;
	    muEtas["sqrt"]               = &twoxf;
	    
	    lnks["identity"]             = &identf;
	    linvs["identity"]            = &identf;
	    muEtas["identity"]           = &onef;
	    
	    lnks["inverse"]              = &inversef;
	    linvs["inverse"]             = &inversef;
	    muEtas["inverse"]            = &invderivf;
	    
	    lnks["logit"]                = &logitLink;
	    linvs["logit"]               = &logitLinkInv;
	    muEtas["logit"]              = &logitMuEta;
	    
	    lnks["probit"]               = &probitLink;
	    linvs["probit"]              = &probitLinkInv;
	    muEtas["probit"]             = &probitMuEta;
	    
//	    lnks["cloglog"]              = &cloglogLink;
	    linvs["cloglog"]             = &cloglogLinkInv;
	    muEtas["cloglog"]            = &cloglogMuEta;
	    
	    devRes["Gamma"]              = &GammaDevRes;
	    varFuncs["Gamma"]            = &sqrf;   // x^2

	    aics["binomial"]             = &BinomialAIC;
	    devRes["binomial"]           = &BinomialDevRes;
	    varFuncs["binomial"]         = &x1mxf;  // x * (1 - x)

	    devRes["gaussian"]           = &GaussianDevRes;
	    varFuncs["gaussian"]         = &onef;   // 1

	    varFuncs["inverse.gaussian"] = &cubef;  // x^3

	    aics["poisson"]              = &PoissonAIC;
	    devRes["poisson"]            = &PoissonDevRes;
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

    VectorXd glmFamily::linkFun(const VectorXd &mu) const {
	VectorXd ans(mu.size());
	if (lnks.count(d_link)) {
	    std::transform(mu.data(), mu.data() + mu.size(), ans.data(), lnks[d_link]);
	} else {
	    NumericVector ans_R = d_linkfun(NumericVector(mu.data(), mu.data() + mu.size()));
	    std::copy(ans_R.begin(), ans_R.end(), ans.data());
	}
	return ans;
    }
    
    VectorXd glmFamily::linkInv(const VectorXd &eta) const {
	VectorXd ans(eta.size());
	if (linvs.count(d_link)) {
	    std::transform(eta.data(), eta.data() + eta.size(), ans.data(), linvs[d_link]);
	} else {
	    NumericVector ans_R = d_linkinv(NumericVector(eta.data(), eta.data() + eta.size()));
	    std::copy(ans_R.begin(), ans_R.end(), ans.data());
	}
	return ans;
    }

    VectorXd glmFamily::muEta(const VectorXd &eta) const {
	VectorXd ans(eta.size());
	if (muEtas.count(d_link)) {
	    std::transform(eta.data(), eta.data() + eta.size(), ans.data(), muEtas[d_link]);
	} else {
	    NumericVector ans_R = d_muEta(NumericVector(eta.data(), eta.data() + eta.size()));
	    std::copy(ans_R.begin(), ans_R.end(), ans.data());
	}
	return ans;
    }
    
    VectorXd glmFamily::variance(const VectorXd &mu) const {
	VectorXd ans(mu.size());
	if (varFuncs.count(d_link)) {
	    std::transform(mu.data(), mu.data() + mu.size(), ans.data(), varFuncs[d_link]);
	} else {
	    NumericVector ans_R = d_variance(NumericVector(mu.data(), mu.data() + mu.size()));
	    std::copy(ans_R.begin(), ans_R.end(), ans.data());
	}
	return ans;
    }
    
    VectorXd
    glmFamily::devResid(const VectorXd &mu, const VectorXd &weights, const VectorXd &y) const {
	int n = mu.size();
	VectorXd ans(n);
	if (devRes.count(d_family)) {
	    double (*f)(const double&, const double&, const double&) = devRes[d_family];
	    const double *mm = mu.data(), *ww = weights.data(), *yy = y.data();
	    double *aa = ans.data();
	    for (int i = 0; i < n; ++i) aa[i] = f(yy[i], mm[i], ww[i]);
	} else {
	    NumericVector ans_R = d_devRes(NumericVector(y.data(), y.data() + n),
					   NumericVector(mu.data(), mu.data() + n),
					   NumericVector(weights.data(), weights.data() + n));
	    std::copy(ans_R.begin(), ans_R.end(), ans.data());
	}
	return ans;
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
