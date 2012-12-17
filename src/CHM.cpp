// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// CHM.cpp: Implementation of the sparse matrix and factor classes
//
// Copyright (C) 2012-2013 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#include "CHM.h"
#include <Rcpp.h>		// for Rcout
using namespace Rcpp;

extern "C" {
    void R_cholmod_error(int status, const char *file, int line,
			 const char *message) {
	if(status < 0) {
	    Rcpp::Rcout << message << " at file " << file
		  << ", line " << line << std::endl;
	    Rcpp::stop("Cholmod error");
	}
	else Rf_warning("Cholmod warning '%s' at file '%s', line %d",
			message, file, line);
    }

    int cholmod_start(CHM_CM Common) {
	int val;
	static int(*fun)(CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (int(*)(CHM_CM))
		R_GetCCallable("Matrix", "cholmod_start");
	val = fun(Common);
	Common->print_function = NULL;
	Common->error_handler = R_cholmod_error;
	return val;
    }

    int cholmod_free_sparse(CHM_SP *A, CHM_CM Common) {
	static int(*fun)(CHM_SP*,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (int(*)(CHM_SP*,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_free_sparse");
	return fun(A, Common);
    }

}

namespace CHM {
    dsCMatrix::dsCMatrix(S4& A)
	: d_Dim(A.slot("Dim")),
	  d_upper(CharacterVector(A.slot("Dim"))[0] == "U"),
	  d_colptr(A.slot("p")),
	  d_rowval(A.slot("i")),
	  d_factors(A.slot("factors")),
	  d_nzval(A.slot("x")),
	  d_sp(new cholmod_sparse) {
	d_sp->nrow = d_Dim[0];
	d_sp->ncol = d_Dim[1];
	d_sp->nzmax = nnz();
	d_sp->p = (void*)&d_colptr[0];
	d_sp->i = (void*)&d_rowval[0];
	d_sp->nz = NULL;
	d_sp->x = (void*)&d_nzval[0];
	d_sp->z = NULL;
	d_sp->stype = d_upper ? 1 : -1;
	d_sp->itype = CHOLMOD_INT;
	d_sp->xtype = CHOLMOD_REAL;
	d_sp->dtype = CHOLMOD_DOUBLE;
	d_sp->sorted = 1;
	d_sp->packed = 1;
    }
    
    dsCMatrix::~dsCMatrix() {delete d_sp;}
    
// still need to implement dsCMatrix::update_factors;
    
    CHMfactor::CHMfactor(S4& L)
	: d_colcount(L.slot("colcount")),
	  d_perm(L.slot("perm")),
	  d_type(L.slot("type")),
	  d_fr(new cholmod_factor) {}
    
    CHMfactor::~CHMfactor() {delete d_fr;}
    
    CHMsimpl::CHMsimpl(S4& L)
	: CHMfactor(L),
	  d_p(L.slot("p")),
	  d_i(L.slot("i")),
	  d_nz(L.slot("nz")),
	  d_nxt(L.slot("nxt")),
	  d_prv(L.slot("prv")) {}
    
    dCHMsimpl::dCHMsimpl(S4 &L) : CHMsimpl(L), d_x(L.slot("x")) {}
    
    CHMsuper::CHMsuper(S4 &L)
	: CHMfactor(L),
	  d_super(L.slot("super")),
	  d_pi(L.slot("pi")),
	  d_px(L.slot("px")),
	  d_s(L.slot("s")) {}
    
    dCHMsuper::dCHMsuper(S4 &L) : CHMsuper(L), d_x(L.slot("x")) {}
}
