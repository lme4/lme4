// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// RSC.h: regular sparse matrix representation
//
// Copyright (C)       2012-2013 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#ifndef LME4_RSC_H
#define LME4_RSC_H

#include <Rcpp.h>

namespace lme4 {    
    class RSC { /**< const parts of mixed-effects predictor in regular sparse column format */
    protected:
	const Rcpp::IntegerMatrix rv; /**< rowvals matrix for Zt */
	const Rcpp::NumericMatrix xv; /**< xvals matrix for ZtXt */
	const Rcpp::NumericVector lower; /**< lower bounds for theta components */
	const int k;  /**< number of random effects per observation */
	const int kpp;    /**< number of rows in xv = k + p */
	const int n;      /**< number of observations */
	const int p;      /**< number of fixed-effects coefficients */
	const int q;      /**< total number of random effects */
    public:
	RSC(const SEXP, const SEXP, const SEXP);
	Rcpp::NumericVector &apply_lambda(const Rcpp::NumericVector&,
					  Rcpp:: NumericVector&) const;
	void update_A(const Rcpp::NumericVector&,
		      const Rcpp::NumericVector&,
		      Rcpp::S4 &AA,
		      Rcpp::NumericVector&) const;
    };
}

#endif // LME4_RSC_H
