#include "glmFamily.h"
#include <limits>

using namespace Rcpp;
using namespace std;

namespace glm {
    // Establish the values for the class constants
    double glmFamily::epsilon = numeric_limits<double>::epsilon();
    
    // initialize the function maps (i.e. associative arrays of functions)
    drmap glmFamily::devRes = drmap();

    fmap glmFamily::linvs = fmap();
    fmap glmFamily::lnks = fmap();
    fmap glmFamily::muEtas = fmap();
    fmap glmFamily::varFuncs = fmap();
    
    void glmFamily::initMaps() {
	// initialize the static maps.  The identity link is
	// guaranteed to be initialized if any maps are initialized
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

	    devRes["binomial"]           = &BinomialDevRes;
	    varFuncs["binomial"]         = &x1mxf;  // x * (1 - x)

	    devRes["gaussian"]           = &GaussianDevRes;
	    varFuncs["gaussian"]         = &onef;   // 1

	    varFuncs["inverse.gaussian"] = &cubef;  // x^3

	    devRes["poisson"]            = &PoissonDevRes;
	    varFuncs["poisson"]          = &identf; // x
    }
    
    glmFamily::glmFamily(List ll) throw (std::runtime_error)
	: lst(ll),
	  //d_family(as<std::string>(CharacterVector(ll["family"]))),
	  //d_link(  as<std::string>(CharacterVector(ll["link"]))),
// I haven't been able to work out an expression to initialize the
// Functions from list components.  This is a placeholder until I can
// do so.
	  d_devRes("c"), d_linkfun("c"), d_linkinv("c"),
	  d_muEta("c"), d_variance("c") {
	  // d_devRes(wrap(ll["dev.resids"])),
	  // d_linkfun(wrap(ll["linkfun"])),
	  // d_linkinv(wrap(ll["linkinv"])),
	  // d_muEta(wrap(ll["mu.eta"])),
	  // d_variance(wrap(ll["variance"])) {
	if (!lst.inherits("family"))
	    throw std::runtime_error("glmFamily requires a list of (S3) class \"family\"");
 	CharacterVector ff = lst["family"], lnk = lst["link"];
 	d_family = as<std::string>(ff);
 	d_link = as<std::string>(lnk);
 	d_linkinv = ll["linkinv"];
 	d_linkfun = ll["linkfun"];
 	d_muEta = ll["mu.eta"];
 	d_variance = ll["variance"];
 	d_devRes = ll["dev.resids"];

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
}
