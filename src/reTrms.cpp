//
// reTrms.cpp: random-effects terms objects
//
// Copyright (C)      2012 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#include "predModule.h"

namespace lme4Eigen {
    reTrms::reTrms(SEXP cnms, SEXP flist)
	: d_cnms(cnms),
	  d_flist(flist),
	  d_assign(d_flist.attr("assign")) {
    }

    Eigen::VectorXi reTrms::ncols() const {
	int nt = d_cnms.size();
	Eigen::VectorXi ans(nt);
	int *ap = ans.data();
	for (int i = 0; i < nt; i++) {
	    Rcpp::CharacterVector nmi = d_cnms(i);
	    ap[i] = nmi.size();
	}
	return ans;
    }

    Eigen::VectorXi reTrms::nctot() const {
	Eigen::VectorXi ans(d_flist.size()), nc = ncols();
	const int *asg = d_assign.begin(), *ncp = nc.data(), nnc(nc.size());
	int *ansp = ans.data();
	ans.setZero();
	for (int i = 0; i < nnc; i++) ansp[asg[i] - 1] += ncp[i];
	return ans;
    }

    Eigen::VectorXi reTrms::nlevs() const {
	int nfac = d_flist.size();
	Eigen::VectorXi ans(nfac);
	int *ap = ans.data();
	for (int i = 0; i < nfac; i++) {
	    Rcpp::IntegerVector   fi = d_flist(i);
	    Rcpp::CharacterVector ll = fi.attr("levels");
	    ap[i] = ll.size();
	}
	return ans;
    }

    Eigen::VectorXi reTrms::offsets() const {
	const int *asgn = d_assign.begin(), ntrm(d_cnms.size());
	Eigen::VectorXi ans(ntrm), nc = ncols(), nl = nlevs();
	int *ansp = ans.data(), offset = 0;
	const int *ncp = nc.data(),  *nlp = nl.data();
	for (int i = 0; i < ntrm; ++i) {
	    ansp[i] = offset;
	    offset += ncp[i] * nlp[asgn[i] - 1];
	}
	return ans;
    }

}
