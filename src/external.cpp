// external.cpp: externally .Call'able functions in lme4
//
// Copyright (C)       2011-2014 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#include <iomanip>
#include "predModule.h"
#include "respModule.h"
#include "optimizer.h"

extern "C" {
    typedef   Eigen::VectorXi                     iVec;
    typedef   Eigen::Map<iVec>                   MiVec;
    typedef   Eigen::MatrixXd                      Mat;
    typedef   Eigen::Map<Mat>                     MMat;
    typedef   Eigen::VectorXd                      Vec;
    typedef   Eigen::Map<Vec>                     MVec;
    typedef   Eigen::ArrayXd                       Ar1;
    typedef   Eigen::Map<Ar1>                     MAr1;
    typedef   Eigen::ArrayXXd                      Ar2;
    typedef   Eigen::Map<Ar2>                     MAr2;
    typedef   Eigen::MappedSparseMatrix<double> MSpMat;
    typedef   Eigen::SparseMatrix<double>        SpMat;
    typedef   Eigen::SimplicialLLT<SpMat>         SLLT;
    typedef   Eigen::LLT<Mat>                      LLT;

    using      Rcpp::CharacterVector;
    using      Rcpp::Environment;
    using      Rcpp::IntegerVector;
    using      Rcpp::Language;
    using      Rcpp::List;
    using      Rcpp::Named;
    using      Rcpp::NumericVector;
    using      Rcpp::NumericMatrix;
    using      Rcpp::Rcout;
    using      Rcpp::XPtr;
    using      Rcpp::as;
    using      Rcpp::wrap;

    using       glm::glmFamily;

    using      lme4::glmResp;
    using      lme4::lmResp;
    using      lme4::lmerResp;
    using      lme4::merPredD;
    using      lme4::nlsResp;

    using optimizer::Golden;
    using optimizer::Nelder_Mead;
    using optimizer::nm_status;

    using      std::runtime_error;
    using      std::invalid_argument;

    // utilities

    // check the sizes of mapped vectors and, if consistent, copy src to dest
    static inline void sized_copy(MVec& dest, const MVec& src) {
        int n = dest.size();
        if (src.size() != n) throw invalid_argument("size mismatch");
        dest = src;
    }

    SEXP allPerm_int(SEXP v_) {
        BEGIN_RCPP;
        iVec     v(as<iVec>(v_));   // forces a copy
        int     sz(v.size());
        std::vector<iVec> vec;

        std::sort(v.data(), v.data() + sz);
        do {
            vec.push_back(iVec(v));
        } while (std::next_permutation(v.data(), v.data() + sz));

        int  nperm(vec.size());
        List allPerm(nperm);
        for (int j = 0; j < nperm; ++j) allPerm[j] = wrap(vec[j]);
        return allPerm;
        END_RCPP;
    }

    SEXP Eigen_SSE() {
        BEGIN_RCPP;
        return wrap(Eigen::SimdInstructionSetsInUse());
        END_RCPP;
    }

    // generalized linear model (and generalized linear mixed model) response

    SEXP glm_Create(SEXP fam, SEXP y, SEXP weights, SEXP offset, SEXP mu,
                    SEXP sqrtXwt, SEXP sqrtrwt, SEXP wtres, SEXP eta, SEXP n) {
        BEGIN_RCPP;
        glmResp *ans = new glmResp(List(fam), y, weights, offset, mu,
                                   sqrtXwt, sqrtrwt, wtres, eta, n);
        return wrap(XPtr<glmResp>(ans, true));
        END_RCPP;
    }

    SEXP glm_aic(SEXP ptr_) {
        BEGIN_RCPP;
        return ::Rf_ScalarReal(XPtr<glmResp>(ptr_)->aic());
        END_RCPP;
    }

    SEXP glm_setN(SEXP ptr_, SEXP n) {
        BEGIN_RCPP;
        XPtr<glmResp>(ptr_)->setN(as<MVec>(n));
        END_RCPP;
    }

    SEXP glm_devResid(SEXP ptr_) {
        BEGIN_RCPP;
        return wrap(XPtr<glmResp>(ptr_)->devResid());
        END_RCPP;
    }

    SEXP glm_family(SEXP ptr_) {
        BEGIN_RCPP;
        return wrap(XPtr<glmResp>(ptr_)->family());
        END_RCPP;
    }

    SEXP glm_link(SEXP ptr_) {
        BEGIN_RCPP;
        return wrap(XPtr<glmResp>(ptr_)->link());
        END_RCPP;
    }

    SEXP glm_muEta(SEXP ptr_) {
        BEGIN_RCPP;
        return wrap(XPtr<glmResp>(ptr_)->muEta());
        END_RCPP;
    }

    SEXP glm_resDev(SEXP ptr_) {
        BEGIN_RCPP;
        return ::Rf_ScalarReal(XPtr<glmResp>(ptr_)->resDev());
        END_RCPP;
    }

    SEXP glm_setTheta(SEXP ptr, SEXP newtheta) {
        BEGIN_RCPP;
        XPtr<glmResp>(ptr)->setTheta(::Rf_asReal(newtheta));
        END_RCPP;
    }

    SEXP glm_sqrtWrkWt(SEXP ptr_) {
        BEGIN_RCPP;
        return wrap(XPtr<glmResp>(ptr_)->sqrtWrkWt());
        END_RCPP;
    }

    SEXP glm_theta(SEXP ptr) {
        BEGIN_RCPP;
        return ::Rf_ScalarReal(XPtr<glmResp>(ptr)->theta());
        END_RCPP;
    }

    SEXP glm_updateWts(SEXP ptr_) {
        BEGIN_RCPP;
        return ::Rf_ScalarReal(XPtr<glmResp>(ptr_)->updateWts());
        END_RCPP;
    }

    SEXP glm_variance(SEXP ptr_) {
        BEGIN_RCPP;
        return wrap(XPtr<glmResp>(ptr_)->variance());
        END_RCPP;
    }

    SEXP glm_wrkResids(SEXP ptr_) {
        BEGIN_RCPP;
        return wrap(XPtr<glmResp>(ptr_)->wrkResids());
        END_RCPP;
    }

    SEXP glm_wrkResp(SEXP ptr_) {
        BEGIN_RCPP;
        return wrap(XPtr<glmResp>(ptr_)->wrkResp());
        END_RCPP;
    }

    SEXP glm_wtWrkResp(SEXP ptr_) {
        BEGIN_RCPP;
        return wrap(XPtr<glmResp>(ptr_)->wtWrkResp());
        END_RCPP;
    }

    SEXP glm_Laplace(SEXP ptr_, SEXP ldL2, SEXP ldRX2, SEXP sqrL) {
        BEGIN_RCPP;
        return ::Rf_ScalarReal(XPtr<glmResp>(ptr_)->Laplace(::Rf_asReal(ldL2),
                                                            ::Rf_asReal(ldRX2),
                                                            ::Rf_asReal(sqrL)));
        END_RCPP;
    }

    SEXP glm_updateMu(SEXP ptr_, SEXP gamma) {
        BEGIN_RCPP;
        return ::Rf_ScalarReal(XPtr<glmResp>(ptr_)->updateMu(as<MVec>(gamma)));
        END_RCPP;
    }

    // glm family objects

    SEXP glmFamily_Create(SEXP fam_) {
        BEGIN_RCPP;
        glmFamily *ans = new glmFamily(List(fam_));
        return wrap(XPtr<glmFamily>(ans, true));
        END_RCPP;
    }

    SEXP glmFamily_link(SEXP ptr, SEXP mu) {
        BEGIN_RCPP;
        return wrap(XPtr<glmFamily>(ptr)->linkFun(as<MVec>(mu)));
        END_RCPP;
    }

    SEXP glmFamily_linkInv(SEXP ptr, SEXP eta) {
        BEGIN_RCPP;
        return wrap(XPtr<glmFamily>(ptr)->linkInv(as<MVec>(eta)));
        END_RCPP;
    }

    SEXP glmFamily_devResid(SEXP ptr, SEXP y, SEXP mu, SEXP wt) {
        BEGIN_RCPP;
        return wrap(XPtr<glmFamily>(ptr)->devResid(as<MVec>(y),
                                                   as<MVec>(mu),
                                                   as<MVec>(wt)));
        END_RCPP;
    }

    SEXP glmFamily_aic(SEXP ptr, SEXP y, SEXP n, SEXP mu, SEXP wt, SEXP dev) {
        BEGIN_RCPP;
        return ::Rf_ScalarReal(XPtr<glmFamily>(ptr)->aic(as<MVec>(y),
                                                         as<MVec>(n),
                                                         as<MVec>(mu),
                                                         as<MVec>(wt),
                                                         ::Rf_asReal(dev)));
        END_RCPP;
    }

    SEXP glmFamily_muEta(SEXP ptr, SEXP eta) {
        BEGIN_RCPP;
        return wrap(XPtr<glmFamily>(ptr)->muEta(as<MVec>(eta)));
        END_RCPP;
    }

    SEXP glmFamily_setTheta(SEXP ptr, SEXP ntheta) {
        BEGIN_RCPP;
        XPtr<glmFamily>(ptr)->setTheta(::Rf_asReal(ntheta));
        END_RCPP;
    }

    SEXP glmFamily_theta(SEXP ptr) {
        BEGIN_RCPP;
        return ::Rf_ScalarReal(XPtr<glmFamily>(ptr)->theta());
        END_RCPP;
    }

    SEXP glmFamily_variance(SEXP ptr, SEXP mu) {
        BEGIN_RCPP;
        return wrap(XPtr<glmFamily>(ptr)->variance(as<MVec>(mu)));
        END_RCPP;
    }

    static inline double pwrss(lmResp *rp, merPredD *pp, double fac) {
        return rp->wrss() + (fac ? pp->sqrL(fac) : pp->u0().squaredNorm());
    }

    static double internal_glmerWrkIter(merPredD *pp, glmResp *rp, bool uOnly) {
        int debug=0; // !=0 to enable
        if (debug) Rcpp::Rcout << "(igWI, pre-updateXwts) Xwts: min: " <<
                       pp->Xwts().minCoeff() <<
                       " sqrtWrkWt: min: " <<
                       rp->sqrtWrkWt().minCoeff() << std::endl;
        pp->updateXwts(rp->sqrtWrkWt());
        if (debug) Rcpp::Rcout << "(igWI) Xwts: min: " <<
                       pp->Xwts().minCoeff() <<
                       " max: " << pp->Xwts().maxCoeff() << std::endl;
        pp->updateDecomp();
        pp->updateRes(rp->wtWrkResp());
        if (uOnly) pp->solveU();
        else pp->solve();
        if (debug) {
            Rcpp::Rcout << "(igWI)" <<
                " delu_min: " << pp->delu().minCoeff() <<
                "; delu_max: " << pp->delu().maxCoeff() <<
                "; delb_min: " << pp->delb().minCoeff() <<
                "; delb_max: " << pp->delb().maxCoeff() <<
                std::endl; // if (verb)
        }
        rp->updateMu(pp->linPred(1.));
        if (debug) Rcpp::Rcout << "(igWI) mu: min: " << rp->mu().minCoeff() <<
                       " max: " << rp->mu().maxCoeff() << std::endl;
        return rp->resDev() + pp->sqrL(1.);
    }

    // FIXME: improve verbose output (remove code, even commented,
    //     intended for finding pointer/referencing/updating bugs;
    //     leave code that allows end-user to see what's going on
    //     in PWRSS iterations)
    //
    //   Separate verb settings for min/max delu, delb 
    //      (and which_min, which_max) vs entire
    //       delu/delb vectors? note length(delu) will generally be 
    //       >> length(delb) ...
    //
    // FIXME: sufficient to print just before/after update?
    //
    // FIXME: allow user-set maxit
    static void pwrssUpdate(glmResp *rp, merPredD *pp, bool uOnly,
                            double tol, int verbose) {
        double oldpdev=std::numeric_limits<double>::max();
        double pdev;
        int maxit = 30, maxstephalfit = 10;
        bool   cvgd = false, verb = verbose > 2, moreverb = verbose > 10;

        pdev = oldpdev; // define so debugging statements work on first step
        for (int i = 0; i < maxit; i++) {
            if (verb) {
                Rcpp::Rcout << "*** pwrssUpdate step " << i << std::endl;
                // Rcpp::Rcout << "\nmin delu at iteration " << i << ": " << pp->delu().minCoeff() << std::endl;
                // Rcpp::Rcout << "\nmax delu at iteration " << i << ": " << pp->delu().maxCoeff() << std::endl;
                // Rcpp::Rcout << "\nresDev before dels, iter:  " << i << ",  " << rp->resDev() << std::endl;
                // FIXME: would like to print this in row, not column, format
                //
                // Rcpp::Rcout << "before update:" << "pdev = " << pdev << std::endl; // if (verb)
            }
            Vec   olddelu(pp->delu()), olddelb(pp->delb());
            pdev = internal_glmerWrkIter(pp, rp, uOnly);
            if (verb) {
                Rcpp::Rcout << "pdev=" << pdev <<
                    "; delu_min: " << pp->delu().minCoeff() <<
                    "; delu_max: " << pp->delu().maxCoeff() <<
                    "; delb_min: " << pp->delb().minCoeff() <<
                    "; delb_max: " << pp->delb().maxCoeff() <<
                    std::endl; // if (verb)
            }
            if (std::abs((oldpdev - pdev) / pdev) < tol) {cvgd = true; break;}

            // if (pdev != pdev) Rcpp::Rcout << "nan detected" << std::endl;
            // if (isnan(pdev)) Rcpp::Rcout << "nan detected" << std::endl;

            // trying to detect nan; may be hard to do it completely portably,
            // and hard to detect in advance (i.e. what conditions lead to
            // nan from internal_glmerWrkIter ... ?)
            // http://stackoverflow.com/questions/570669/checking-if-a-double-or-float-is-nan-in-c
            // check use of isnan() in base R code, or other Rcpp code??
#define isNAN(a)  (a!=a)
            if (isNAN(pdev) || (pdev > oldpdev)) {
                // PWRSS step led to _larger_ deviation, or nan; try step halving
                if (verb) Rcpp::Rcout <<
                              "\npwrssUpdate: Entering step halving loop" 
                                      << std::endl;
                for (int k = 0; k < maxstephalfit &&
                         (isNAN(pdev) || (pdev > oldpdev)); k++) {
                    pp->setDelu((olddelu + pp->delu())/2.);
                    if (!uOnly) pp->setDelb((olddelb + pp->delb())/2.);
                    rp->updateMu(pp->linPred(1.));
                    pdev = rp->resDev() + pp->sqrL(1.);
                    if (moreverb) {
                        Rcpp::Rcout << "step-halving iteration " <<
                            k << ":  pdev=" << pdev <<
                            "; delu_min: " << pp->delu().minCoeff() <<
                            "; delu_max: " << pp->delu().maxCoeff() <<
                            "; delb_min: " << pp->delb().minCoeff() <<
                            "; delb_max: " << pp->delb().maxCoeff() <<
                            std::endl;
                    } // if (moreverb)
                }
                if (isNAN(pdev) || ((pdev - oldpdev) > tol) )
                    // FIXME: fill in max halfsetp iters in error statement
                    throw runtime_error("(maxstephalfit) PIRLS step-halvings failed to reduce deviance in pwrssUpdate");
            } // step-halving
            oldpdev = pdev;
        } // pwrss loop
        if (!cvgd)
            // FIXME: fill in max iters in error statement
            throw runtime_error("pwrssUpdate did not converge in (maxit) iterations");
    }

    SEXP glmerLaplace(SEXP pp_, SEXP rp_, SEXP nAGQ_, SEXP tol_, SEXP verbose_) {
        BEGIN_RCPP;
        XPtr<glmResp>  rp(rp_);
        XPtr<merPredD> pp(pp_);

        if ( ::Rf_asInteger(verbose_) >100) {
            Rcpp::Rcout << "\nglmerLaplace resDev:  " << rp->resDev() << std::endl;
            Rcpp::Rcout << "\ndelb 1:  " << pp->delb() << std::endl;
        }
        pwrssUpdate(rp, pp, ::Rf_asInteger(nAGQ_), ::Rf_asReal(tol_), ::Rf_asInteger(verbose_));
        return ::Rf_ScalarReal(rp->Laplace(pp->ldL2(), pp->ldRX2(), pp->sqrL(1.)));
        END_RCPP;
    }

    static Ar1 devcCol(const MiVec& fac, const Ar1& u, const Ar1& devRes) {
        Ar1  ans(u.square());
        for (int i = 0; i < devRes.size(); ++i) ans[fac[i] - 1] += devRes[i];
        return ans;
    }

    static double sqrt2pi = std::sqrt(2. * PI);

    SEXP glmerAGQ(SEXP pp_, SEXP rp_, SEXP tol_, SEXP GQmat_, SEXP fac_, SEXP verbose_) {
        BEGIN_RCPP;

        XPtr<glmResp>     rp(rp_);
        XPtr<merPredD>    pp(pp_);
        const MiVec      fac(as<MiVec>(fac_));
        double           tol(::Rf_asReal(tol_));
        double          verb(::Rf_asReal(verbose_));
        if (fac.size() != rp->mu().size())
            throw std::invalid_argument("size of fac must match dimension of response vector");

        pwrssUpdate(rp, pp, true, tol, verb); // should be a no-op
        const Ar1      devc0(devcCol(fac, pp->u(1.), rp->devResid()));

        const unsigned int q(pp->u0().size());
        if (pp->L().factor()->nzmax !=  q)
            throw std::invalid_argument("AGQ only defined for a single scalar random-effects term");
        const Ar1         sd(MAr1((double*)pp->L().factor()->x, q).inverse());

        const MMat     GQmat(as<MMat>(GQmat_));
        Ar1             mult(q);

        mult.setZero();
        for (int i = 0; i < GQmat.rows(); ++i) {
            double     zknot(GQmat(i, 0));
            if (zknot == 0)
                mult += Ar1::Constant(q, GQmat(i, 1));
            else {
                pp->setU0(zknot * sd); // to be added to current delu
                rp->updateMu(pp->linPred(1.));
                mult += (-0.5 * (devcCol(fac, pp->u(1.), rp->devResid()) - devc0) -
                         GQmat(i, 2)).exp() * GQmat(i, 1)/sqrt2pi;
            }
        }
        pp->setU0(Vec::Zero(q)); // restore settings from pwrssUpdate;
        rp->updateMu(pp->linPred(1.));
        return ::Rf_ScalarReal(devc0.sum() + pp->ldL2() - 2 * std::log(mult.prod()));
        END_RCPP;
    }

    void nstepFac(nlsResp *rp, merPredD *pp, int verb) {
        double prss0(pwrss(rp, pp, 0.));

        for (double fac = 1.; fac > 0.001; fac /= 2.) {
            double prss1 = rp->updateMu(pp->linPred(fac)) + pp->sqrL(fac);
            if (verb > 3)
                ::Rprintf("prss0=%10g, diff=%10g, fac=%6.4f\n",
                          prss0, prss0 - prss1, fac);
            if (prss1 < prss0) {
                pp->installPars(fac);
                return;
            }
        }
        throw runtime_error("step factor reduced below 0.001 without reducing pwrss");
    }

#define NMAXITER 300
    static void prssUpdate(nlsResp *rp, merPredD *pp, int verb, bool uOnly, double tol) {
        bool cvgd(false);
        for (int it=0; it < NMAXITER; ++it) {
            rp->updateMu(pp->linPred(0.));
            pp->updateXwts(rp->sqrtXwt());
            pp->updateDecomp();
            pp->updateRes(rp->wtres());
            double ccrit((uOnly ? pp->solveU() : pp->solve())/pwrss(rp, pp, 0.));
            if (verb > 3)
                ::Rprintf("ccrit=%10g, tol=%10g\n", ccrit, tol);
            if (ccrit < tol) {
                cvgd = true;
                break;
            }
            nstepFac(rp, pp, verb);
        }
        if (!cvgd) throw runtime_error("prss failed to converge in 300 iterations");
    }

    SEXP nlmerLaplace(SEXP pp_, SEXP rp_, SEXP theta_, SEXP u0_, SEXP beta0_,
                      SEXP verbose_, SEXP uOnly_, SEXP tol_) {
        BEGIN_RCPP;

        XPtr<nlsResp>     rp(rp_);
        XPtr<merPredD>    pp(pp_);
        pp->setTheta(as<MVec>(theta_));
        pp->setU0(as<MVec>(u0_));
        pp->setBeta0(as<MVec>(beta0_));
        prssUpdate(rp, pp, ::Rf_asInteger(verbose_), ::Rf_asLogical(uOnly_),
                    ::Rf_asReal(tol_));
        return ::Rf_ScalarReal(rp->Laplace(pp->ldL2(), pp->ldRX2(), pp->sqrL(1.)));

        END_RCPP;
    }

    SEXP golden_Create(SEXP lower_, SEXP upper_) {
        BEGIN_RCPP;
        Golden *ans = new Golden(::Rf_asReal(lower_), ::Rf_asReal(upper_));
        return wrap(XPtr<Golden>(ans, true));
        END_RCPP;
    }

    SEXP golden_newf(SEXP ptr_, SEXP f_) {
        BEGIN_RCPP;
        XPtr<Golden>(ptr_)->newf(::Rf_asReal(f_));
        END_RCPP;
    }

    SEXP golden_xeval(SEXP ptr_) {
        BEGIN_RCPP;
        return wrap(XPtr<Golden>(ptr_)->xeval());
        END_RCPP;
    }

    SEXP golden_value(SEXP ptr_) {
        BEGIN_RCPP;
        return wrap(XPtr<Golden>(ptr_)->value());
        END_RCPP;
    }

    SEXP golden_xpos(SEXP ptr_) {
        BEGIN_RCPP;
        return wrap(XPtr<Golden>(ptr_)->xpos());
        END_RCPP;
    }

    SEXP isNullExtPtr(SEXP Ptr) {
        void *ptr = R_ExternalPtrAddr(Ptr);
        return ::Rf_ScalarLogical(ptr == (void*)NULL);
    }

    // linear model response (also the base class for other response classes)

    SEXP lm_updateWrss(Environment rho) {
        const MAr1  y(as<MAr1>(rho.get("y")));
        const MAr1  mu(as<MAr1>(rho.get("mu")));
        const MAr1  sqrtrwt(as<MAr1>(rho.get("sqrtrwt")));
        MVec        wtres(as<MVec>(rho.get("wtres")));
        wtres = sqrtrwt.cwiseProduct(y - mu);
        SEXP        ans = ::Rf_ScalarReal(wtres.squaredNorm());
        rho.assign("wrss", ans);
        return ans;
    }

    SEXP lm_setOffset(SEXP rho_, SEXP oo_) {
        BEGIN_RCPP;
        Environment  rho(rho_);
        const MAr1   oo(as<MAr1>(oo_));
        MAr1         offset(as<MAr1>(rho.get("offset")));
        if (oo.size() != offset.size())
            throw invalid_argument("size mismatch in setOffset");
        offset = oo;
        return lm_updateWrss(rho);
        END_RCPP;
    }

    SEXP lm_setResp(SEXP rho_, SEXP rr_) {
        BEGIN_RCPP;
        Environment  rho(rho_);
        const MAr1   rr(as<MAr1>(rr_));
        MAr1         y(as<MAr1>(rho.get("y")));
        if (rr.size() != y.size())
            throw invalid_argument("size mismatch in setResp");
        y = rr;
        return lm_updateWrss(rho);
        END_RCPP;
    }

    SEXP lm_setWeights(SEXP rho_, SEXP wts_) {
        BEGIN_RCPP;
        Environment  rho(rho_);
        const MAr1   wts(as<MAr1>(wts_));
        MAr1         weights(as<MAr1>(rho.get("weights")));
        if (wts.size() != weights.size())
            throw invalid_argument("size mismatch in setWeights");
        weights = wts;
        as<MAr1>(rho.get("sqrtrwt")) = wts.sqrt();
        rho.assign("ldW",wrap(wts.log().sum()));
        return lm_updateWrss(rho);
        END_RCPP;
    }

    SEXP lm_updateMu(SEXP rho_, SEXP gam_) {
        BEGIN_RCPP;
        Environment  rho(rho_);
        const MAr1   gam(as<MAr1>(gam_));
        MAr1         mu(as<MAr1>(rho.get("mu")));
        if (gam.size() != mu.size())
            throw invalid_argument("size mismatch in updateMu");
        mu = gam;
        return lm_updateWrss(rho);
        END_RCPP;
    }

    // linear mixed-effects model response

    SEXP lmer_Laplace(SEXP rho_, SEXP ldL2_, SEXP ldRX2_, SEXP sqrL_, SEXP sigma_sq_) {
        BEGIN_RCPP;
        Environment   rho(rho_);
        int n     = ::Rf_length(rho.get("y")),
            REML  = ::Rf_asInteger(rho.get("REML"));
        double
            ldL2  = ::Rf_asReal(ldL2_),
            ldRX2 = ::Rf_asReal(ldRX2_),
            nmp   = n - REML,
            sigsq = ::Rf_isNull(sigma_sq_) ? 1.0 : ::Rf_asReal(sigma_sq_),
            sqrL  = ::Rf_asReal(sqrL_),
            wrss  = ::Rf_asReal(rho.get("wrss"));

        double result = nmp * (2.0 * M_LN_SQRT_2PI + std::log(sigsq)); // (2pi sigma_sq)^-df/2
        result += (wrss + sqrL) / sigsq; // exp(-1/2sigma_sq x |pwrss|)
        result += ldL2 + (REML > 0 ? ldRX2 : 0.0); // det|LL'|^-1/2 and similar REML penalty
        result -= ::Rf_asReal(rho.get("ldW")); // subtract prior weights factor
        return ::Rf_ScalarReal(result);
        END_RCPP;
    }

    static double lmer_dev(Environment pred, Environment resp, const MVec& theta) {
        int debug=0;
        double val = 1.;
/*
        ppt->setTheta(theta);

        ppt->updateXwts(rpt->sqrtXwt());
        ppt->updateDecomp();
        rpt->updateMu(ppt->linPred(0.));
        ppt->updateRes(rpt->wtres());
        ppt->solve();
        rpt->updateMu(ppt->linPred(1.));
        val=rpt->Laplace(ppt->ldL2(), ppt->ldRX2(), ppt->sqrL(1.));
        if (debug) {
            Rcpp::Rcout.precision(10);
            Rcpp::Rcout << "lmer_dev: theta=" <<
                ppt->theta() << ", val=" << val << std::endl;
        }
*/
        return val;
    }

    SEXP lmer_Deviance(SEXP pptr_, SEXP rptr_, SEXP theta_) {
        BEGIN_RCPP;
        XPtr<lmerResp>   rpt(rptr_);
        XPtr<merPredD>   ppt(pptr_);
        return ::Rf_ScalarReal(lmer_dev(ppt, rpt, as<MVec>(theta_)));
        END_RCPP;
    }

    SEXP lmer_opt1(SEXP pptr_, SEXP rptr_, SEXP lower_, SEXP upper_) {
        BEGIN_RCPP;
        XPtr<lmerResp>     rpt(rptr_);
        XPtr<merPredD>     ppt(pptr_);
        Eigen::VectorXd     th(1);
        optimizer::Golden gold(::Rf_asReal(lower_), ::Rf_asReal(upper_));
        for (int i = 0; i < 30; ++i) {
            th[0] = gold.xeval();
            gold.newf(lmer_dev(ppt, rpt, th));
        }
        return List::create(Named("theta") = ::Rf_ScalarReal(gold.xpos()),
                            Named("objective") = ::Rf_ScalarReal(gold.value()));
        END_RCPP;
    }

    // dense predictor module for mixed-effects models
				// methods that can be called by other methods
    SEXP predD_mkL(SEXP rho_) {
	// create the SLLT and LLT objects and assign external pointers in rho_
        BEGIN_RCPP;
        Environment  rho(rho_);
        const MSpMat Zt(as<MSpMat>(rho.get("Zt")));
        SLLT        *L = new SLLT();
        L->setShift(1.0);
        L->compute(Zt * Zt.adjoint());
        rho.assign("Lptr",wrap(XPtr<SLLT>(L)));
        LLT         *RX = new LLT();
        const MMat   X(as<MMat>(rho.get("X")));
        RX->compute(X.adjoint() * X);
        rho.assign("RXpt", wrap(XPtr<LLT>(RX)));
        END_RCPP;
    }

    SEXP predD_updtL(SEXP rho_) {
	// update Lambdat, LambdatZt and L
        BEGIN_RCPP;
        Environment   rho(rho_);
        const MSpMat  Zt(as<MSpMat>(rho.get("Zt")));
        MSpMat        Lambdat(as<MSpMat>(rho.get("Lambdat")));
        const MiVec   Lind(as<MiVec>(rho.get("Lind")));
        const MVec    theta(as<MVec>(rho.get("theta")));
        if (Lind.size() != Lambdat.nonZeros()) throw invalid_argument("updtL: size mismatch");
        double *LamX = Lambdat.valuePtr();
        for (int i = 0; i < Lind.size(); ++i) {
            LamX[i] = theta(Lind(i) - 1);
        }
        SpMat         LamtUt(Lambdat * Zt);
	rho.assign("LamtUt",wrap(LamtUt));
        XPtr<SLLT>(rho.get("Lptr"))->factorize(LamtUt * LamtUt.adjoint());
        END_RCPP;
    }

    SEXP predD_updtRX(SEXP rho_) {
        BEGIN_RCPP;
        Environment   rho(rho_);
        const MMat    UtV(as<MMat>(rho.get("UtV")));
        const MSpMat  Lambdat(as<MSpMat>(rho.get("Lambdat")));
        SLLT         *L(XPtr<SLLT>(rho.get("Lptr")));
        MMat          RZX(as<MMat>(rho.get("RZX")));
        Mat           pr(L->permutationP() * (Lambdat * UtV));
        RZX = L->matrixL().solve(pr);
        XPtr<LLT>(rho.get("RXpt"))->compute(as<MMat>(rho.get("VtV")) - RZX.adjoint() * RZX);
        END_RCPP;
    }
                                // setters
    SEXP predD_setBeta0(SEXP rho_, SEXP b0_) {
        BEGIN_RCPP;
        MVec        beta0(as<MVec>(Environment(rho_).get("beta0")));
        sized_copy(beta0, as<MVec>(b0_));
        END_RCPP;
    }

    SEXP predD_setDelb(SEXP rho_, SEXP db_) {
        BEGIN_RCPP;
        MVec        delb(as<MVec>(Environment(rho_).get("delb")));
        sized_copy(delb, as<MVec>(db_));
        END_RCPP;
    }

    SEXP predD_setDelu(SEXP rho_, SEXP du_) {
        BEGIN_RCPP;
        MVec        delu(as<MVec>(Environment(rho_).get("delu")));
        sized_copy(delu, as<MVec>(du_));
        END_RCPP;
    }

    SEXP predD_setTheta(SEXP rho_, SEXP th_) {
        BEGIN_RCPP;
        MVec        theta(as<MVec>(Environment(rho_).get("theta")));
        sized_copy(theta, as<MVec>(th_));
	return predD_updtL(rho_);
        END_RCPP;
    }
                                // getters
    SEXP predD_L(SEXP rho_) {
        BEGIN_RCPP;
        return wrap(XPtr<SLLT>(Environment(rho_).get("Lptr"))->matrixL().derived());
        END_RCPP;
    }

    SEXP predD_P(SEXP rho_) {
        BEGIN_RCPP;
        return wrap(XPtr<SLLT>(Environment(rho_).get("Lptr"))->permutationP().indices());
        END_RCPP;
    }

    SEXP predD_RX(SEXP rho_) {
        BEGIN_RCPP;
        return wrap(Mat(XPtr<LLT>(Environment(rho_).get("RXpt"))->matrixU()));
        END_RCPP;
    }

    SEXP predD_RXdiag(SEXP rho_) {
        BEGIN_RCPP;
        return wrap(Mat(XPtr<LLT>(Environment(rho_).get("RXpt"))->matrixU()).diagonal());
        END_RCPP;
    }

    SEXP predD_RXi(SEXP rho_) {
        BEGIN_RCPP;
        LLT           *RXpt(XPtr<LLT>(Environment(rho_).get("RXpt")));
        int            p(RXpt->rows());
        return wrap(Mat(RXpt->matrixU().solve(Mat::Identity(p,p))));
        END_RCPP;
    }

    SEXP predD_ldL2(SEXP rho_) {
        BEGIN_RCPP;
        return ::Rf_ScalarReal(log(XPtr<SLLT>(Environment(rho_).get("Lptr"))->determinant()));
        END_RCPP;
    }

    SEXP predD_ldRX2(SEXP rho_) {
        BEGIN_RCPP;
        return ::Rf_ScalarReal(
            Mat(XPtr<LLT>(Environment(rho_).get("RXpt"))->matrixL()).diagonal().log().sum());
        END_RCPP;
    }

    SEXP predD_unsc(SEXP rho_) {
        BEGIN_RCPP;
        LLT           *RXpt(XPtr<LLT>(Environment(rho_).get("RXpt")));
        int            p(RXpt->rows());
        return wrap(Mat(RXpt->solve(Mat::Identity(p,p))));
        END_RCPP;
    }

                                // methods
    SEXP predD_CcNumer(SEXP rho_) {
        BEGIN_RCPP;
// FIXME: needs a body	
//        return ::Rf_ScalarReal(XPtr<merPredD>(ptr)->CcNumer());
        END_RCPP;
    }

    SEXP predD_condVar(SEXP rho_, SEXP respEnv_) {
        BEGIN_RCPP;
// FIXME: needs a body
//        return wrap(XPtr<merPredD>(ptr)->condVar(Environment(rho)));
        END_RCPP;
    }

    SEXP predD_solve(SEXP rho_) {
        BEGIN_RCPP;
	Environment    rho(rho_);
	const MSpMat   Lamt(as<MSpMat>(rho.get("Lambdat")));
	const SLLT    *Lp(XPtr<SLLT>(rho.get("Lptr")));
	const Vec      PLamtUtr(Lp->permutationP()*(Lamt * as<MVec>(rho.get("Utr"))));
	Vec            cu(Lp->matrixL().solve(PLamtUtr));
	const MMat     RZX(as<MMat>(rho.get("RZX")));
	MVec           delb(as<MVec>(rho.get("delb")));
	delb = XPtr<LLT>(rho.get("RXpt"))->solve(as<MVec>(rho.get("Vtr")) - RZX.adjoint() * cu);
	MVec           delu(as<MVec>(rho.get("delu")));
	delu = Lp->permutationPinv() * Lp->matrixU().solve(cu - RZX * delb);
        END_RCPP;
    }

    SEXP predD_solveU(SEXP rho_) {
        BEGIN_RCPP;
	Environment    rho(rho_);
	const SLLT    *Lp(XPtr<SLLT>(rho.get("Lptr")));
	MVec           delu(as<MVec>(rho.get("delu")));
	delu = Lp->solve(as<MVec>(rho.get("Utr")));
        END_RCPP;
    }

    static inline Vec internal_uvec(Environment rho, double fac) {
        return as<Vec>(rho.get("u0")) += fac * as<MVec>(rho.get("delu"));
    }

    static inline Vec internal_b(Environment rho, double fac) {
	return as<MSpMat>(rho.get("Lambdat")).adjoint() * internal_uvec(rho, fac);
    }

    static inline Vec internal_beta(Environment rho, double fac) {
	return as<Vec>(rho.get("beta0")) += fac * as<MVec>(rho.get("delb"));
    }

    static inline Vec internal_linpred(Environment rho, double fac) {
	return as<MMat>(rho.get("X")) * internal_beta(rho,fac) +
	    as<MSpMat>(rho.get("Zt")).adjoint() * internal_b(rho,fac);
    }

    SEXP predD_b(SEXP rho_, SEXP fac_) {
        BEGIN_RCPP;
	return wrap(internal_b(Environment(rho_), ::Rf_asReal(fac_)));
        END_RCPP;
    }

    SEXP predD_beta(SEXP rho_, SEXP fac_) {
        BEGIN_RCPP;
	return wrap(internal_beta(Environment(rho_), ::Rf_asReal(fac_)));
        END_RCPP;
    }

    SEXP predD_installPars(SEXP rho_, SEXP fac_) {
        BEGIN_RCPP;
        Environment rho(rho_);
        double      fac(::Rf_asReal(fac_));
        MVec        beta0(as<MVec>(rho.get("beta0")));
        MVec        delb(as<MVec>(rho.get("delb")));
        MVec        u0(as<MVec>(rho.get("u0")));
        MVec        delu(as<MVec>(rho.get("delu")));
        beta0 += fac * delb;
        u0 += fac * delu;
        delb.setZero();
        delu.setZero();
        END_RCPP;
    }

    SEXP predD_linPred(SEXP rho_, SEXP fac_) {
        BEGIN_RCPP;
	return wrap(internal_linpred(Environment(rho_), ::Rf_asReal(fac_)));
        END_RCPP;
    }

    SEXP predD_sqrL(SEXP rho_, SEXP fac_) {
        BEGIN_RCPP;
        return ::Rf_ScalarReal(internal_uvec(Environment(rho_),::Rf_asReal(fac_)).squaredNorm());
        END_RCPP;
    }

    SEXP predD_u(SEXP rho_, SEXP fac_) {
        BEGIN_RCPP;
        return wrap(internal_uvec(Environment(rho_),::Rf_asReal(fac_)));
        END_RCPP;
    }

    SEXP predD_updtRes(SEXP rho_, SEXP wres_) {
        BEGIN_RCPP;
	Environment rho(rho_);
	const MVec  wres(as<MVec>(wres_));
// in-place update works for Vtr but not Utr, not sure why
	rho.assign("Utr", wrap(as<MSpMat>(rho.get("Zt")) * wres));
	MVec        Vtr(as<MVec>(rho.get("Vtr")));
	Vtr = as<MMat>(rho.get("X")).adjoint() * wres;
        END_RCPP;
    }

    SEXP predD_updtXwts(SEXP rho_, SEXP wts_) {
        BEGIN_RCPP;
// needs a body
        END_RCPP;
    }

    SEXP NelderMead_Create(SEXP lb_, SEXP ub_, SEXP xstep0_, SEXP x_, SEXP xtol_) {
        BEGIN_RCPP;
        MVec  lb(as<MVec>(lb_)), ub(as<MVec>(ub_)), xstep0(as<MVec>(xstep0_)), x(as<MVec>(x_));
        Nelder_Mead *ans =
            new Nelder_Mead(lb, ub, xstep0, x, optimizer::nl_stop(as<MVec>(xtol_)));
        return wrap(XPtr<Nelder_Mead>(ans, true));
        END_RCPP;
    }

    SEXP NelderMead_newf(SEXP ptr_, SEXP f_) {
        BEGIN_RCPP;
        switch (XPtr<Nelder_Mead>(ptr_)->newf(::Rf_asReal(f_))) {
        case optimizer::nm_evals:         return ::Rf_ScalarInteger(-4);
        case optimizer::nm_forced:        return ::Rf_ScalarInteger(-3);
        case optimizer::nm_nofeasible:    return ::Rf_ScalarInteger(-2);
        case optimizer::nm_x0notfeasible: return ::Rf_ScalarInteger(-1);
        case optimizer::nm_active:        return ::Rf_ScalarInteger(0);
        case optimizer::nm_minf_max:      return ::Rf_ScalarInteger(1);
        case optimizer::nm_fcvg:          return ::Rf_ScalarInteger(2);
        case optimizer::nm_xcvg:          return ::Rf_ScalarInteger(3);
        }
        END_RCPP;
    }

    SEXP NelderMead_setForce_stop(SEXP ptr_, SEXP stp_) {
        BEGIN_RCPP;
        XPtr<Nelder_Mead>(ptr_)->setForce_stop(::Rf_asLogical(stp_));
        END_RCPP;
    }

    SEXP NelderMead_setFtol_abs(SEXP ptr_, SEXP fta_) {
        BEGIN_RCPP;
        XPtr<Nelder_Mead>(ptr_)->setFtol_rel(::Rf_asReal(fta_));
        END_RCPP;
    }

    SEXP NelderMead_setFtol_rel(SEXP ptr_, SEXP ftr_) {
        BEGIN_RCPP;
        XPtr<Nelder_Mead>(ptr_)->setFtol_rel(::Rf_asReal(ftr_));
        END_RCPP;
    }

    SEXP NelderMead_setIprint(SEXP ptr_, SEXP ip_) {
        BEGIN_RCPP;
        XPtr<Nelder_Mead>(ptr_)->set_Iprint(::Rf_asInteger(ip_));
        END_RCPP;
    }

    SEXP NelderMead_setMaxeval(SEXP ptr_, SEXP mm_) {
        BEGIN_RCPP;
        XPtr<Nelder_Mead>(ptr_)->set_Maxeval(::Rf_asInteger(mm_));
        END_RCPP;
    }

    SEXP NelderMead_setMinf_max(SEXP ptr_, SEXP mm_) {
        BEGIN_RCPP;
        XPtr<Nelder_Mead>(ptr_)->setMinf_max(::Rf_asReal(mm_));
        END_RCPP;
    }

    SEXP NelderMead_xeval(SEXP ptr_) {
        BEGIN_RCPP;
        return wrap(XPtr<Nelder_Mead>(ptr_)->xeval());
        END_RCPP;
    }

    SEXP NelderMead_value(SEXP ptr_) {
        BEGIN_RCPP;
        return wrap(XPtr<Nelder_Mead>(ptr_)->value());
        END_RCPP;
    }

    SEXP NelderMead_xpos(SEXP ptr_) {
        BEGIN_RCPP;
        return wrap(XPtr<Nelder_Mead>(ptr_)->xpos());
        END_RCPP;
    }

                                // return the number of function evaluations performed
    SEXP NelderMead_evals(SEXP ptr_) { 
        BEGIN_RCPP;
        return wrap(int(XPtr<Nelder_Mead>(ptr_)->evals()));
        END_RCPP;
    }

    // nonlinear model response (also the base class for other response classes)

    SEXP nls_Create(SEXP y, SEXP weights, SEXP offset, SEXP mu, SEXP sqrtXwt,
                    SEXP sqrtrwt, SEXP wtres, SEXP gamma, SEXP mod, SEXP env, SEXP pnms) {
        BEGIN_RCPP;
        nlsResp *ans =
            new nlsResp(y, weights, offset, mu, sqrtXwt, sqrtrwt, wtres, gamma, mod, env, pnms);
        return wrap(XPtr<nlsResp>(ans, true));
        END_RCPP;
    }

    SEXP nls_Laplace(SEXP ptr_, SEXP ldL2, SEXP ldRX2, SEXP sqrL) {
        BEGIN_RCPP;
        return ::Rf_ScalarReal(XPtr<nlsResp>(ptr_)->
                               Laplace(::Rf_asReal(ldL2),
                                       ::Rf_asReal(ldRX2),
                                       ::Rf_asReal(sqrL)));
        END_RCPP;
    }

    SEXP nls_updateMu(SEXP ptr_, SEXP gamma) {
        BEGIN_RCPP;
        return ::Rf_ScalarReal(XPtr<nlsResp>(ptr_)->updateMu(as<MVec>(gamma)));
        END_RCPP;
    }

    SEXP showlocation(SEXP obj) {
        int ll = Rf_length(obj);
        if (Rf_isReal(obj)) {
            double *vv = REAL(obj);
            Rcpp::Rcout << "Numeric vector of length " << ll
                        << " at location: " << vv << std::endl;
            if (ll > 0) {
                Rcpp::Rcout << "Values: " << vv[0];
                for(int i = 1; i < std::min(ll, 5); ++i)
                    Rcpp::Rcout << "," << vv[i];
                if (ll > 8) Rcpp::Rcout << ",...,";
                for (int i = std::max(5, ll - 3); i < ll; ++i)
                    Rcpp::Rcout << "," << vv[i];
                Rcpp::Rcout << std::endl;
            }
        }
        if (Rf_isInteger(obj)) {
            int *vv = INTEGER(obj);
            Rcpp::Rcout << "Numeric vector of length " << ll
                        << " at location: " << vv << std::endl;
            if (ll > 0) {
                Rcpp::Rcout << "Values: " << vv[0];
                for(int i = 1; i < std::min(ll, 5); ++i)
                    Rcpp::Rcout << "," << vv[i];
                if (ll > 8) Rcpp::Rcout << ",...,";
                for (int i = std::max(5,ll - 3); i < ll; ++i)
                    Rcpp::Rcout << "," << vv[i];
                Rcpp::Rcout << std::endl;
            }
        }
        return R_NilValue;
    }
}

#include <R_ext/Rdynload.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {

    CALLDEF(Eigen_SSE, 0),

    CALLDEF(allPerm_int, 1),

    CALLDEF(glm_Create, 10),    // generate external pointer

    CALLDEF(glm_setN, 2),       // setters
    CALLDEF(glm_setTheta,       2),

    CALLDEF(glm_aic,            1), // getters
    CALLDEF(glm_devResid,       1),
    CALLDEF(glm_family,         1),
    CALLDEF(glm_link,           1),
    CALLDEF(glm_muEta,          1),
    CALLDEF(glm_resDev,         1),
    CALLDEF(glm_sqrtWrkWt,      1),
    CALLDEF(glm_theta,          1),
    CALLDEF(glm_variance,       1),
    CALLDEF(glm_wtWrkResp,      1),
    CALLDEF(glm_wrkResids,      1),
    CALLDEF(glm_wrkResp,        1),

    CALLDEF(glm_Laplace,        4), // methods
    CALLDEF(glm_updateMu,       2),
    CALLDEF(glm_updateWts,      1),

    CALLDEF(glmFamily_Create,   1), // generate external pointer

    CALLDEF(glmFamily_aic,      6), // methods
    CALLDEF(glmFamily_link,     2),
    CALLDEF(glmFamily_linkInv,  2),
    CALLDEF(glmFamily_devResid, 4),
    CALLDEF(glmFamily_muEta,    2),
    CALLDEF(glmFamily_setTheta, 2),
    CALLDEF(glmFamily_theta,    1),
    CALLDEF(glmFamily_variance, 2),

    CALLDEF(glmerAGQ,           6),
    CALLDEF(glmerLaplace,       5),

    CALLDEF(golden_Create,      2),
    CALLDEF(golden_newf,        2),
    CALLDEF(golden_value,       1),
    CALLDEF(golden_xeval,       1),
    CALLDEF(golden_xpos,        1),

    CALLDEF(isNullExtPtr, 1),

    CALLDEF(lm_setOffset,       2), // setters
    CALLDEF(lm_setResp,         2),
    CALLDEF(lm_setWeights,      2),

    CALLDEF(lm_updateMu,        2), // methods
    CALLDEF(lm_updateWrss,      1),

    CALLDEF(lmer_Deviance,      3), // methods
    CALLDEF(lmer_Laplace,       5),
    CALLDEF(lmer_opt1,          4),

    CALLDEF(predD_setBeta0,     2), // setters
    CALLDEF(predD_setDelb,      2),
    CALLDEF(predD_setDelu,      2),
    CALLDEF(predD_setTheta,     2),

    CALLDEF(predD_CcNumer,      1), // getters
    CALLDEF(predD_L,            1),
    CALLDEF(predD_P,            1),
    CALLDEF(predD_RX,           1),
    CALLDEF(predD_RXdiag,       1),
    CALLDEF(predD_RXi,          1),
    CALLDEF(predD_ldL2,         1),
    CALLDEF(predD_ldRX2,        1),
    CALLDEF(predD_unsc,         1),

    CALLDEF(predD_b,            2), // methods
    CALLDEF(predD_beta,         2),
    CALLDEF(predD_condVar,      2),
    CALLDEF(predD_linPred,      2),
    CALLDEF(predD_installPars,  2),
    CALLDEF(predD_mkL,          1),
    CALLDEF(predD_solve,        1),
    CALLDEF(predD_solveU,       1),
    CALLDEF(predD_sqrL,         2),
    CALLDEF(predD_u,            2),
    CALLDEF(predD_updtL,        1),
    CALLDEF(predD_updtRX,       1),
    CALLDEF(predD_updtRes,      2),
    CALLDEF(predD_updtXwts,     2),

    CALLDEF(NelderMead_Create,  5),
    CALLDEF(NelderMead_newf,    2),
    CALLDEF(NelderMead_setForce_stop, 2),
    CALLDEF(NelderMead_setFtol_abs, 2),
    CALLDEF(NelderMead_setFtol_rel, 2),
    CALLDEF(NelderMead_setIprint, 2),
    CALLDEF(NelderMead_setMaxeval, 2),
    CALLDEF(NelderMead_setMinf_max, 2),
    CALLDEF(NelderMead_value,   1),
    CALLDEF(NelderMead_xeval,   1),
    CALLDEF(NelderMead_xpos,    1),

    CALLDEF(nlmerLaplace,       8),

    CALLDEF(nls_Create,        11), // generate external pointer

    CALLDEF(nls_Laplace,        4), // methods
    CALLDEF(nls_updateMu,       2),

    CALLDEF(showlocation,       1),
    {NULL, NULL, 0}
};

/** Initializer for lme4, called upon loading the package.
 *
 *  Register routines that can be called directly from R.
 *  Initialize CHOLMOD and require the LL' form of the factorization.
 *  Install the symbols to be used by functions in the package.
 */
extern "C"
void R_init_lme4(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, (Rboolean)FALSE);
}
