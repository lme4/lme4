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
//    void parchk(int kind, int m, double alpha, double beta);
//    void cdgqf(int nt, int kind, double alpha, double beta, std::vector<double>& t, 
//	       std::vector<double>& wts);
    void cgqf(int nt, int kind, double alpha, double beta, double a, double b, 
	      std::vector<double>& t, std::vector<double>& wts);
//    double class_matrix(int kind, int m, double alpha, double beta, std::vector<double>& aj, 
//			std::vector<double>& bj);
//    void imtqlx(int n, std::vector<double>& d, std::vector<double>& e, std::vector<double>& z);
//    void sgqf(int nt, std::vector<double>& aj, std::vector<double>& bj, double zemu, std::vector<double>& t, 
//	      std::vector<double>& wts);
}

#endif
