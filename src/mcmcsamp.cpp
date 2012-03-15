//
// mcmcsamp.cpp: implementation of mcmcsamp and related classes using Eigen
//
// Copyright (C) 2011-2012 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#include "mcmcsamp.h"

namespace lme4 {
    using Rcpp::as;

    static inline double pwrss(lme4::merPredD *pred, lme4::lmResp *resp) {
	return pred->sqrL(1.) + resp->wrss();
    }

    static inline double sigmaML(lme4::merPredD *pred, lme4::lmResp *resp) {
	return std::sqrt(pwrss(pred, resp)/double(resp->y().size()));
    }

    mcmcsamp::mcmcsamp(lme4::merPredD *pred, lme4::lmResp *resp,
		       SEXP dev, SEXP fixef, SEXP sigma, SEXP ranef)
	: d_dev(  as<MVec>(dev)),
	  d_fixef(as<MMat>(fixef)),
	  d_sigma(as<MVec>(sigma)),
	  d_ranef(as<MMat>(ranef))
    {
	Rcpp::RNGScope scope;	// handles the getRNGstate/putRNGstate
	bool           sig(  d_sigma.size() > 0);
	bool           rr(   d_ranef.rows() > 0);
	int            n(    resp->y().size());
	int            nsamp(d_dev.size());
	int            nth(  pred->theta().size());
	int            p(    pred->beta0().size());
	int            q(    pred->u0().size());
	double         npq(  n + q);
	double         lsigma(sig ? sigmaML(pred, resp) : 1.);
	double         wrss;

	if (d_fixef.cols() != nsamp || d_fixef.rows() != p || (sig && d_sigma.size() != nsamp) ||
	    (ranef && (d_ranef.cols() != nsamp || d_ranef.rows() != p)))
	    throw std::invalid_argument("dimension mismatch");
	if (nth > 1) ::Rf_error("only handling the simple (nth == 1) cases now");

	for (int k = 0; k < nsamp; ++k) {
	    pred->updateDecomp();
	    pred->solve();
	    pred->MCMC_beta_u(lsigma);
	    d_fixef.col(k) = pred->beta(1.);
	    if (rr) d_ranef.col(k) = pred->b(1.);
	    wrss = resp->updateMu(pred->linPred(1.));
	    if (sig) d_sigma[k] = lsigma = std::sqrt((pred->sqrL(1.) + wrss)/::Rf_rchisq(npq));
	}
    }

#if 0
/**
 * Generate a Markov-chain Monte Carlo sample from an mer object
 *
 * @param x pointer to an merMCMC object
 * @param fm pointer to an mer object
 *
 * @return x with samples filled in
 */
    SEXP mer_MCMCsamp(SEXP , SEXP fm)
    {
	SEXP devsamp = GET_SLOT(x, lme4_devianceSym);
	int *dims = DIMS_SLOT(x), nsamp = LENGTH(devsamp);
	int n = dims[n_POS], np = dims[np_POS],
	    p = dims[p_POS], q = dims[q_POS];
	double
	    *STsamp = REAL(GET_SLOT(x, lme4_STSym)),
	    *d = DEV_SLOT(fm), *dev = REAL(devsamp),
	    *sig = SLOT_REAL_NULL(x, lme4_sigmaSym),
	    *fixsamp = FIXEF_SLOT(x), *resamp = RANEF_SLOT(x);

	GetRNGstate();
	/* The first column of storage slots contains the fitted values */
	for (int i = 1; i < nsamp; i++) {
/* FIXME: This is probably wrong for a model with weights. */
	    if (sig) 		/* update and store sigma */
		sig[i] = sqrt(d[pwrss_POS]/rchisq((double)(n + q)));
	    /* update L, RX, beta, u, eta, mu, res and d */
	    MCMC_beta_u(fm, sig ? sig[i] : 1, fixsamp + i * p,
			resamp + (resamp ? i : 0) * q);
	    dev[i] = d[ML_POS];
	    /* update theta_T, theta_S and A */
	    MCMC_T(fm, sig ? sig[i] : 1);
	    MCMC_S(fm, sig ? sig[i] : 1);
	    ST_getPars(fm, STsamp + i * np); /* record theta */
	}
	PutRNGstate();
	/* Restore pars from the first columns of the samples */
	Memcpy(FIXEF_SLOT(fm), fixsamp, p);
	ST_setPars(fm, STsamp);
	update_ranef(fm);

	return x;
    }

/**
 * Update the fixed effects and the orthogonal random effects in an MCMC sample
 * from an mer object.
 *
 * @param x an mer object
 * @param sigma current standard deviation of the per-observation
 *        noise terms.
 * @param fvals pointer to memory in which to store the updated beta
 * @param rvals pointer to memory in which to store the updated b (may
 *              be (double*)NULL)
 */
    static void MCMC_beta_u(SEXP x, double sigma, double *fvals, double *rvals)
    {
	int *dims = DIMS_SLOT(x);
	int i1 = 1, p = dims[p_POS], q = dims[q_POS];
	double *V = V_SLOT(x), *fixef = FIXEF_SLOT(x), *muEta = MUETA_SLOT(x),
	    *u = U_SLOT(x), mone[] = {-1,0}, one[] = {1,0};
	CHM_FR L = L_SLOT(x);
	double *del1 = Calloc(q, double), *del2 = Alloca(p, double);
	CHM_DN sol, rhs = N_AS_CHM_DN(del1, q, 1);
	R_CheckStack();

	if (V || muEta) {
	    error(_("Update not yet written"));
	} else {			/* Linear mixed model */
	    update_L(x);
	    update_RX(x);
	    lmm_update_fixef_u(x);
	    /* Update beta */
	    for (int j = 0; j < p; j++) del2[j] = sigma * norm_rand();
	    F77_CALL(dtrsv)("U", "N", "N", &p, RX_SLOT(x), &p, del2, &i1);
	    for (int j = 0; j < p; j++) fixef[j] += del2[j];
	    /* Update u */
	    for (int j = 0; j < q; j++) del1[j] = sigma * norm_rand();
	    F77_CALL(dgemv)("N", &q, &p, mone, RZX_SLOT(x), &q,
			    del2, &i1, one, del1, &i1);
	    sol = M_cholmod_solve(CHOLMOD_Lt, L, rhs, &c);
	    for (int j = 0; j < q; j++) u[j] += ((double*)(sol->x))[j];
	    M_cholmod_free_dense(&sol, &c);
	    update_mu(x);	     /* and parts of the deviance slot */
	}
	Memcpy(fvals, fixef, p);
	if (rvals) {
	    update_ranef(x);
	    Memcpy(rvals, RANEF_SLOT(x), q);
	}
	Free(del1);
    }

/**
 * Update the theta_T parameters from the ST arrays in place.
 *
 * @param x an mer object
 * @param sigma current standard deviation of the per-observation
 *        noise terms.
 */
/* FIXME: Probably should fold this function into MCMC_S */
    static void MCMC_T(SEXP x, double sigma)
    {
	int *Gp = Gp_SLOT(x), nt = (DIMS_SLOT(x))[nt_POS];
	double **st = Alloca(nt, double*);
	int *nc = Alloca(nt, int), *nlev = Alloca(nt, int);
	R_CheckStack();

	if (ST_nc_nlev(GET_SLOT(x, lme4_STSym), Gp, st, nc, nlev) < 2) return;
	error("Code for non-trivial theta_T not yet written");
    }

/**
 * Update the theta_S parameters from the ST arrays in place.
 *
 * @param x an mer object
 * @param sigma current standard deviation of the per-observation
 *        noise terms.
 */
    static void MCMC_S(SEXP x, double sigma)
    {
	CHM_SP A = A_SLOT(x), Zt = Zt_SLOT(x);
	int *Gp = Gp_SLOT(x), *ai = (int*)(A->i),
	    *ap = (int*)(A->p), *dims = DIMS_SLOT(x), *perm = PERM_VEC(x);
	int annz = ap[A->ncol], info, i1 = 1, n = dims[n_POS],
	    nt = dims[nt_POS], ns, p = dims[p_POS], pos,
	    q = dims[q_POS], znnz = ((int*)(Zt->p))[Zt->ncol];
	double *R, *ax = (double*)(A->x), *b = RANEF_SLOT(x),
	    *eta = ETA_SLOT(x), *offset = OFFSET_SLOT(x),
	    *rr, *ss, one = 1, *u = U_SLOT(x), *y = Y_SLOT(x);
	int *nc = Alloca(nt, int), *nlev = Alloca(nt, int),
	    *spt = Alloca(nt + 1, int);
	double **st = Alloca(nt, double*);
	R_CheckStack();

	ST_nc_nlev(GET_SLOT(x, lme4_STSym), Gp, st, nc, nlev);
	ns = 0;			/* ns is length(theta_S) */
	spt[0] = 0;			/* pointers into ss for terms */
	for (int i = 0; i < nt; i++) {
	    ns += nc[i];
	    spt[i + 1] = spt[i] + nc[i];
	}

	if (annz == znnz) { /* Copy Z' to A unless A has new nonzeros */
	    Memcpy(ax, (double*)(Zt->x), znnz);
	} else error("Code not yet written for MCMC_S with NLMMs");
	/* Create T'Zt in A */
	Tt_Zt(A, Gp, nc, nlev, st, nt);
	/* Create P'u in ranef slot */
	for (int i = 0; i < q; i++) b[perm[i]] = u[i];
	/* Create X\beta + offset in eta slot */
	for (int i = 0; i < n; i++) eta[i] = offset ? offset[i] : 0;
	F77_CALL(dgemv)("N", &n, &p, &one, X_SLOT(x), &n,
			FIXEF_SLOT(x), &i1, &one, eta, &i1);
	/* Allocate R, rr and ss */
	R = Alloca(ns * ns, double); /* crossproduct matrix then factor */
	rr = Alloca(ns, double);	 /* row of model matrix for theta_S */
	ss = Alloca(ns, double);	 /* right hand side, then theta_S */
	R_CheckStack();
	AZERO(R, ns * ns);
	AZERO(ss, ns);
	/* Accumulate crossproduct from pseudo-data part of model matrix */
	for (int i = 0; i < q; i++) {
	    int sj = theta_S_ind(i, nt, Gp, nlev, spt);
	    AZERO(rr, ns);
	    rr[sj] = b[i];
	    F77_CALL(dsyr)("U", &ns, &one, rr, &i1, R, &ns);
	}
	/* Accumulate crossproduct and residual product of the model matrix. */
	/* This is done one row at a time.  Rows of the model matrix
	 * correspond to columns of T'Zt */
	for (int j = 0; j < n; j++) { /* jth column of T'Zt */
	    AZERO(rr, ns);
	    for (int p = ap[j]; p < ap[j + 1]; p++) {
		int i = ai[p];	/* row in T'Zt */
		int sj = theta_S_ind(i, nt, Gp, nlev, spt);

		rr[sj] += ax[p] * b[i];
		ss[sj] += rr[sj] * (y[j] - eta[j]);
	    }
	    F77_CALL(dsyr)("U", &ns, &one, rr, &i1, R, &ns);
	}
	F77_CALL(dposv)("U", &ns, &i1, R, &ns, ss, &ns, &info);
	if (info)
	    error(_("Model matrix for theta_S is not positive definite, %d."), info);
	for (int j = 0; j < ns; j++) rr[j] = sigma * norm_rand();
	/* Sample from the conditional Gaussian distribution */
	F77_CALL(dtrsv)("U", "N", "N", &ns, R, &ns, rr, &i1);
	for (int j = 0; j < ns; j++) ss[j] += rr[j];
	/* Copy positive part of solution onto diagonals of ST */
	pos = 0;
	for (int i = 0; i < nt; i++) {
	    for (int j = 0; j < nc[i]; j++) {
		st[i][j * (nc[i] + 1)] = (ss[pos] > 0) ? ss[pos] : 0;
		pos++;
	    }
	}
	update_A(x);
    }

#endif
}
