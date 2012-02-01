// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// reTrms.h: random-effects terms objects
//
// Copyright (C)       2012 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#ifndef LME4_RETRMS_H
#define LME4_RETRMS_H

#include <RcppEigen.h>

namespace lme4Eigen {
    class reTrms {
    protected:
	Rcpp::List            d_cnms;
	Rcpp::List           d_flist;
	Rcpp::IntegerVector d_assign;
    public:
	reTrms(SEXP, SEXP);

	const Rcpp::List&          flist() const {return d_flist;}
	const Rcpp::List&           cnms() const {return d_cnms;}

	Eigen::VectorXi            nctot() const;
	Eigen::VectorXi            ncols() const;
	Eigen::VectorXi            nlevs() const;
	Eigen::VectorXi          offsets() const;
	Eigen::VectorXi            terms(const int&) const;
    };
}

#endif // LME4_RETRMS_H
