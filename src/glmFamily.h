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

    class glmDist {
    protected:
                                //@{ R functions from the family, as a fall-back
        Rcpp::Function d_devRes, d_variance, d_aic;
                                //@}
        Rcpp::Environment d_rho;
        std::string d_familyName;    /**< $family, for error messages only */
        bool   d_hasDispersionField; /**< does the R family list have a $dispersion
                                       *   component (?family, R >= 4.3.0)? */
        double d_dispersionField;    /**< its value, if present: NA_real_ (free) or
                                       *   a fixed numeric value; meaningless otherwise */
    public:
        glmDist(Rcpp::List&);
        virtual ~glmDist() {}

        virtual const ArrayXd variance(const ArrayXd&) const;
        virtual const ArrayXd devResid(const ArrayXd&, const ArrayXd&, const ArrayXd&) const;
        virtual double             aic(const ArrayXd&, const ArrayXd&, const ArrayXd&,
                                       const ArrayXd&, double) const;
        /**< in keeping with the botched up nomenclature in the R glm function,
         *   the value of aic is the deviance */
        virtual double           theta() const;
        virtual void          setTheta(const double&);
        /**< whether this family has a free/estimated dispersion parameter
         *   (Gamma, gaussian, inverse.gaussian) vs. a fixed one (poisson,
         *   binomial, negative binomial) -- mirrors hasNoScale() in
         *   R/utilities.R. The base implementation is reached only for
         *   families glmFamily's constructor didn't recognize by name
         *   (i.e. user-supplied/external family objects). Those are
         *   resolved generically via the family's own $dispersion
         *   component (?family, R >= 4.3.0: NA_real_ if free, a fixed
         *   numeric value otherwise) rather than guessed; if that
         *   component is itself absent (pre-4.3.0-style or incomplete
         *   custom family objects), there is no way to tell, so this
         *   throws an error. */
        virtual bool  hasFreeDispersion() const {
            if (!d_hasDispersionField)
                throw std::runtime_error(
                    "can't tell if dispersion is fixed or free for family '" +
                    d_familyName + "': it has no 'dispersion' component "
                    "(see the 'dispersion' item under 'Value' in ?stats::family)");
            return ISNA(d_dispersionField);
        }
    };

    class binomialDist : public glmDist {
    public:
        binomialDist(Rcpp::List& ll) : glmDist(ll) {}
        const ArrayXd variance(const ArrayXd&) const;
        const ArrayXd devResid(const ArrayXd&, const ArrayXd&, const ArrayXd&) const;
        double             aic(const ArrayXd&, const ArrayXd&, const ArrayXd&,
                               const ArrayXd&, double) const;
        bool  hasFreeDispersion() const {return false;}
    };

    class gammaDist : public glmDist {
    public:
        gammaDist(Rcpp::List& ll) : glmDist(ll) {}
        const ArrayXd variance(const ArrayXd&) const;
        const ArrayXd devResid(const ArrayXd&, const ArrayXd&, const ArrayXd&) const;
        double             aic(const ArrayXd&, const ArrayXd&, const ArrayXd&,
                               const ArrayXd&, double) const;
        bool  hasFreeDispersion() const {return true;}
    };

    class GaussianDist : public glmDist {
    public:
        GaussianDist(Rcpp::List& ll) : glmDist(ll) {}
        const ArrayXd variance(const ArrayXd&) const;
        const ArrayXd devResid(const ArrayXd&, const ArrayXd&, const ArrayXd&) const;
        double             aic(const ArrayXd&, const ArrayXd&, const ArrayXd&,
                               const ArrayXd&, double) const;
        bool  hasFreeDispersion() const {return true;}
    };

    class inverseGaussianDist : public glmDist {
    public:
        inverseGaussianDist(Rcpp::List& ll) : glmDist(ll) {}
        const ArrayXd variance(const ArrayXd&) const;
        const ArrayXd devResid(const ArrayXd&, const ArrayXd&, const ArrayXd&) const;
        double             aic(const ArrayXd&, const ArrayXd&, const ArrayXd&,
                               const ArrayXd&, double) const;
        bool  hasFreeDispersion() const {return true;}
    };

    class negativeBinomialDist : public glmDist {
    protected:
        double  d_theta;
    public:
        negativeBinomialDist (Rcpp::List& ll);
        const ArrayXd variance(const ArrayXd&) const;
        const ArrayXd devResid(const ArrayXd&, const ArrayXd&, const ArrayXd&) const;
        double             aic(const ArrayXd&, const ArrayXd&, const ArrayXd&,
                               const ArrayXd&, double) const;
        double           theta() const {return d_theta;}
        void          setTheta(const double& ntheta) {d_theta = ntheta;}
        bool  hasFreeDispersion() const {return false;}
    };

    class PoissonDist : public glmDist {
    public:
        PoissonDist(Rcpp::List& ll) : glmDist(ll) {}
        const ArrayXd variance(const ArrayXd&) const;
        const ArrayXd devResid(const ArrayXd&, const ArrayXd&, const ArrayXd&) const;
        double             aic(const ArrayXd&, const ArrayXd&, const ArrayXd&,
                               const ArrayXd&, double) const;
        bool  hasFreeDispersion() const {return false;}
    };

    class glmLink {
    protected:
                                //@{ R functions from the family, as a fall-back
        Rcpp::Function    d_linkFun, d_linkInv, d_muEta;
                                //@}
        Rcpp::Environment d_rho;
    public:
        glmLink(Rcpp::List&);
        virtual ~glmLink() {}

        virtual const ArrayXd linkFun(const ArrayXd&) const;
        virtual const ArrayXd linkInv(const ArrayXd&) const;
        virtual const ArrayXd   muEta(const ArrayXd&) const;
    };

    class cauchitLink : public glmLink {
    public:
        cauchitLink(Rcpp::List& ll) : glmLink(ll) {}

        const ArrayXd linkFun(const ArrayXd&) const;
        const ArrayXd linkInv(const ArrayXd&) const;
        const ArrayXd   muEta(const ArrayXd&) const;
    };

    class cloglogLink : public glmLink {
    public:
        cloglogLink(Rcpp::List& ll) : glmLink(ll) {}

//      const ArrayXd linkFun(const ArrayXd&) const;
        const ArrayXd linkInv(const ArrayXd&) const;
        const ArrayXd   muEta(const ArrayXd&) const;
    };

    class identityLink : public glmLink {
    public:
        identityLink(Rcpp::List& ll) : glmLink(ll) {}

        const ArrayXd linkFun(const ArrayXd&) const;
        const ArrayXd linkInv(const ArrayXd&) const;
        const ArrayXd   muEta(const ArrayXd&) const;
    };

    class  inverseLink : public glmLink {
    public:
        inverseLink(Rcpp::List& ll) : glmLink(ll) {}

        const ArrayXd linkFun(const ArrayXd&) const;
        const ArrayXd linkInv(const ArrayXd&) const;
        const ArrayXd   muEta(const ArrayXd&) const;
    };

    class logLink : public glmLink {
    public:
        logLink(Rcpp::List& ll) : glmLink(ll) {}

        const ArrayXd linkFun(const ArrayXd&) const;
        const ArrayXd linkInv(const ArrayXd&) const;
        const ArrayXd   muEta(const ArrayXd&) const;
    };

    class logitLink : public glmLink {
    public:
        logitLink(Rcpp::List& ll) : glmLink(ll) {}

        const ArrayXd linkFun(const ArrayXd&) const;
        const ArrayXd linkInv(const ArrayXd&) const;
        const ArrayXd   muEta(const ArrayXd&) const;
    };

    class probitLink : public glmLink {
    public:
        probitLink(Rcpp::List& ll) : glmLink(ll) {}

        const ArrayXd linkFun(const ArrayXd&) const;
        const ArrayXd linkInv(const ArrayXd&) const;
        const ArrayXd   muEta(const ArrayXd&) const;
    };

    class glmFamily {
    protected:
        std::string  d_family, d_linknam; /**< as in the R glmFamily object */
        glmDist     *d_dist;
        glmLink     *d_link;
    public:
        glmFamily(Rcpp::List ll);
        ~glmFamily();           /**< explicit destructor to call delete on d_dist and d_link */
        const std::string& fam() const {return d_family;}
        const std::string& lnk() const {return d_linknam;}

        //@{ Application of functions from the family using compiled code when available
        const ArrayXd devResid(const ArrayXd&, const ArrayXd&, const ArrayXd&) const;
        const ArrayXd  linkFun(const ArrayXd&  mu) const {return d_link->linkFun(mu);}
        const ArrayXd  linkInv(const ArrayXd& eta) const {return d_link->linkInv(eta);}
        const ArrayXd    muEta(const ArrayXd& eta) const {return d_link->muEta(eta);}
        const ArrayXd variance(const ArrayXd&  mu) const {return d_dist->variance(mu);}
        double             aic(const ArrayXd&, const ArrayXd&, const ArrayXd&,
                               const ArrayXd&, double) const;
        double           theta() const {return d_dist->theta();}
        void          setTheta(const double& theta) {d_dist->setTheta(theta);}
        bool  hasFreeDispersion() const {return d_dist->hasFreeDispersion();}
        //@}
    };

}
    
#endif /* LME4_GLMFAMILY_H */

