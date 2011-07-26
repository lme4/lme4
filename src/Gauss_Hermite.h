// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_GHQ_H
#define LME4_GHQ_H

#include <Rcpp.h>

namespace lme4Eigen {
    class GHQ {
	Rcpp::NumericVector d_xvals, d_wts;
    public:
	GHQ(int);

	const Rcpp::NumericVector& xvals() const {return d_xvals;}
	const Rcpp::NumericVector&   wts() const {return d_wts;}	
    };
}

#endif
