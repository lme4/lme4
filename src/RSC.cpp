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
    RSC::RSC(const IntegerMatrix& rv, const NumericMatrix& xv,
	     const NumericVector& theta, const NumericVector& lower,
	     S4& A4, NumericVector& ubeta)
	: d_rv(rv),
	  d_xv(xv),
	  d_theta(theta),
	  d_lower(lower),
	  d_A4(A4),
	  d_ubeta(ubeta),
	  d_k(d_rv.nrow()),
	  d_kpp(d_xv.nrow()),
	  d_n(d_xv.ncol()),
	  d_p(d_kpp - d_k),
	  d_q((*std::max_element(d_rv.begin(), d_rv.end())) + 1 - d_p) {
	Rcout << "k = " << d_k << ", kpp = " << d_kpp << ", n = " << d_n
	      << ", p = " << d_p << ", q = " << d_q << std::endl;
	if (d_rv.ncol() != d_n || d_p < 0) stop("dimension mismatch of rv and xv");
	if (d_k != std::count_if(d_lower.begin(), d_lower.end(), ::R_finite))
	    stop("dimension mismatch of rv and lower");
	const IntegerVector  dd(d_A4.slot("Dim"));
	if (dd[0] != d_q + d_p) stop("size of A must be q + p");
    }

    NumericVector& RSC::apply_lambda(NumericVector &dest) const {
	int dpos(-1);
	double *rr(0);		// initialize to NULL pointer
	for (int k = 0; k < d_kpp; ++k) {
	    if (d_lower[k] == 0.) { // diagonal element of a factor
		dest[++dpos] *= d_theta[k];
		rr = &dest[k];
	    }
	    else dest[dpos] += d_theta[k] * *++rr;
	}
	return dest;
    }

    void RSC::update_A(const NumericVector &resid) {
	if (resid.size() != d_n) stop("Dimension of resid should be n");
	const IntegerVector &rowval(d_A4.slot("i")), &colptr(d_A4.slot("p"));
	NumericVector nzval(d_A4.slot("x"));
	NumericVector w(d_kpp);
	// initializations
	std::fill(nzval.begin(), nzval.end(), double(0)); // zero the contents of A
	std::fill(d_ubeta.begin(), d_ubeta.end(), double(0)); // and ubeta
	for (int i = 0; i < d_q; ++i) {  // initialize Z part of A to the identity
	    int ll(colptr[i + 1] - 1); // index of last element in column i
	    int ii(rowval[ll]);	   // should be the diagonal element
	    if (ii != i) stop("A is not stored as the upper triangle");
	    nzval[ii] = 1.;
	}
	for (int j = 0; j < d_n; ++j)	{ // iterate over columns of ZtXt
	    double rj(resid[j]);
	    std::copy(&(d_xv(0,j)), &(d_xv(0,j)) + d_kpp, w.begin()); // copy j'th column
	    apply_lambda(w);
	    for (int i = 0; i < d_kpp; ++i) d_ubeta[d_rv(i,j)] += rj * w[i];
	    // scan up the j'th column of ZtXt, which makes
	    // it easier to evaluate the upper triangle of A
	    for (int i = d_kpp - 1; i >= 0; --i) {
		int ii(d_rv(i, j));	// row of ZtXt (column of m) 
		int cpi(colptr[ii]), ll(colptr[ii + 1] - 1); // ll should be location of diagonal
		if (rowval[ll] != ii) stop("u is not upper triangular");
		nzval[ll] += w[i] * w[i];
		for (int l = i - 1; l >= 0 && ll >= cpi; --l) { // off-diagonals in m
		    int ii1(d_rv(l,j)); // row in ZtXt
		    while (rowval[ll] > ii1) --ll; // up to the desired row
		    if (rowval[ll] != ii1) stop("Pattern mismatch");
		    nzval[ll] += w[i] * w[l];
		}
	    }
	}
    }
}
