// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// RSC.cpp: Implementation of the regular sparse matrix representation
//
// Copyright (C) 2012-2013 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#include "RSC.h"
#include "CHM.h"

using namespace Rcpp;

namespace lme4 {
    RSC::RSC(S4 obj)
        : d_obj(obj),
          d_i(as<SEXP>(d_obj.slot("i"))),
          d_x(as<SEXP>(d_obj.slot("x"))),
          d_theta(d_obj.slot("theta")),
          d_lower(d_obj.slot("lower")),
          d_A4(S4(as<SEXP>(d_obj.slot("A")))),
          d_ubeta(d_obj.slot("ubeta")),
          d_uboff(d_obj.slot("uboff")),
          d_k(d_i.nrow()),
          d_kpp(d_x.nrow()),
          d_n(d_x.ncol()),
          d_p(d_kpp - d_k),
          d_qpp(d_ubeta.size()),
          d_q(d_qpp - d_p) {
//	Rcout << "k = " << d_k << ", kpp = " << d_kpp << ", n = " << d_n
//	      << ", p = " << d_p << ", qpp = " << d_qpp << ", p = " << d_p
//	      << std::endl;
        if (d_k != std::count_if(d_lower.begin(), d_lower.end(), ::R_finite))
            stop("number of finite elements in lower does not match nrow(i)");
        const IntegerVector dd(d_A4.slot("Dim"));
        if (dd[0] != d_qpp) stop("size of A must be q + p");
    }

    NumericVector& RSC::apply_lambda(NumericVector& dest, double wt) const {
        int dpos(-1);
        double *rr(0);          // initialize to NULL pointer
        for (int k = 0; k < d_kpp; ++k) {
            if (d_lower[k] == 0.) { // diagonal element of a factor
                dest[++dpos] *= d_theta[k];
                rr = &dest[k];
            }
            else dest[dpos] += d_theta[k] * *++rr;
        }
	if (wt != 1.0) for (int k = 0; k < dest.size(); ++k) dest[k] *= wt;
        return dest;
    }
    
    // This version reevaluates d_A from d_i, d_x and weights
    List RSC::update_A(const NumericVector& theta,
		       const NumericVector& resid,
		       const NumericVector& weights) {
	if (theta.size() != d_theta.size()) stop("Incorrect size of theta");
	for (int i = 0; i < theta.size(); ++i)
	    if (theta[i] < d_lower[i]) stop("theta violates lower bounds");
	std::copy(theta.begin(), theta.end(), d_theta.begin());
        if (resid.size() != d_n) stop("Dimension of resid should be n");
        if (weights.size() != d_n) stop("Dimension of weights should be n");
        const IntegerVector &rowval(d_A4.slot("i")), &colptr(d_A4.slot("p"));
        NumericVector nzval(d_A4.slot("x"));
        NumericVector w(d_kpp);
        NumericMatrix ZtXtr(d_qpp, 1);
        // initializations
        std::fill(nzval.begin(), nzval.end(), double(0)); // zero the contents of A
        for (int i = 0; i < d_q; ++i) {  // initialize Z part of A to the identity
            int ll(colptr[i + 1] - 1); // index of last element in column i
            if (rowval[ll] != i) stop("A is not stored as the upper triangle");
            nzval[ll] = 1.;
        }
        for (int j = 0; j < d_n; ++j)   { // iterate over columns of ZtXt
            std::copy(&d_x(0,j), &d_x(0,j) + d_kpp, w.begin());
            apply_lambda(w, weights[j]);
            for (int i = 0; i < d_kpp; ++i) ZtXtr(rv(i,j),0) += resid[j] * w[i];
            // scan up j'th column of ZtXt (easier to evaluate A's upper triangle)
            for (int i = d_kpp - 1; i >= 0; --i) {
                int ii(rv(i, j));      // row of ZtXt (column of m) 
                int cpi(colptr[ii]), ll(colptr[ii + 1] - 1); // ll is location of diagonal
                nzval[ll] += w[i] * w[i];
                for (int l = i - 1; l >= 0 && ll >= cpi; --l) { // off-diagonals in m
                    int ii1(rv(l,j)); // row in ZtXt
                    while (rowval[ll] > ii1) --ll; // up to the desired row
                    if (rowval[ll] != ii1) stop("Pattern mismatch");
                    nzval[ll] += w[i] * w[l];
                }
            }
        }
        CHM::dsCMatrix A(d_A4);
        A.update_factors();
        NumericVector ld(A.Ldiag());
        double ldL2(0.), ldRX2(0.);
        for (int i = 0; i < d_q; ++i) ldL2 += std::log(ld[i] * ld[i]);
        for (int i = d_q; i < d_qpp; ++i) ldRX2 += std::log(ld[i] * ld[i]);
        NumericVector sol(A.solve(ZtXtr, CHOLMOD_A));
        return List::create(Named("ldL2", ldL2),
                            Named("ldRX2", ldRX2),
                            Named("del_ubeta", sol),
                            Named("linpred", fitted(sol)));
    }

    // This version does not incorporate weights and updates the
    // Cholesky factor by scaling a copy of d_A.
    List RSC::scale_update_A(const NumericVector& theta,
			     const NumericMatrix& ZtXtr) {
	if (theta.size() != d_theta.size()) stop("Incorrect size of theta");
	for (int i = 0; i < theta.size(); ++i) {
	    if (d_lower[i] >= 0.) stop("Vector-valued random effects not allowed");
	    if (theta[i] < 0.) stop("theta violates lower bounds");
	}
	std::copy(theta.begin(), theta.end(), d_theta.begin());
	if (ZtXtr.nrow() != d_qpp) stop("Incorrect size for ZtXtr");
        CHM::dsCMatrix A(d_A4);
        A.scale_update_factors(d_uboff, theta);
        NumericVector ld(A.Ldiag());
        double ldL2(0.), ldRX2(0.);
        for (int i = 0; i < d_q; ++i) ldL2 += std::log(ld[i] * ld[i]);
        for (int i = d_q; i < d_qpp; ++i) ldRX2 += std::log(ld[i] * ld[i]);
				// need to add scaling of ZtXtr here.
        NumericVector sol(A.solve(ZtXtr, CHOLMOD_A));
        return List::create(Named("ldL2", ldL2),
                            Named("ldRX2", ldRX2),
                            Named("del_ubeta", sol),
                            Named("linpred", fitted(sol)));
    }

    NumericVector RSC::Ldiag() {
        CHM::dsCMatrix A(d_A4);
        return A.Ldiag();
    }

    NumericVector RSC::fitted(const NumericVector& ub) const {
        NumericVector w(d_kpp);
        NumericVector ans(d_n);
        for (int j = 0; j < d_n; ++j) {
            double *pt(&ans[j]);
            *pt = 0.;		// not sure this is necessary
            std::copy(&d_x(0,j), &d_x(0,j) + d_kpp, w.begin());
            apply_lambda(w, 1.);
            for (int i = 0; i < d_kpp; ++i) *pt += w[i] * ub[rv(i,j)];
        }
        return ans;
    }
}
