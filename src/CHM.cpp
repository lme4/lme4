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
	    Rcpp::Rcerr << message << " at file " << file
			<< ", line " << line << std::endl;
	    Rcpp::stop("Cholmod error");
	}
	else ::Rf_warning("Cholmod warning '%s' at file '%s', line %d",
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

    CHM_SP cholmod_aat(const_CHM_SP A, int *fset, size_t fsize,
		       int mode, CHM_CM Common) {
	static CHM_SP(*fun)(const_CHM_SP,int*,size_t,
			    int,CHM_CM) = NULL;
	if(fun == NULL)
	    fun = (CHM_SP(*)(const_CHM_SP,int*,size_t,
			     int,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_aat");
	return fun(A, fset, fsize, mode, Common);
    }

    int cholmod_band_inplace(CHM_SP A, int k1, int k2, int mode,
			     CHM_CM Common) {
	static int(*fun)(CHM_SP,int,int,int,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (int(*)(CHM_SP,int,int,int,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_band_inplace");
	return fun(A, k1, k2, mode, Common);
    }

    CHM_SP cholmod_add(const_CHM_SP A, const_CHM_SP B,
		       double alpha[2], double beta[2], int values,
		       int sorted, CHM_CM Common) {
	static CHM_SP(*fun)(const_CHM_SP,const_CHM_SP,
			    double*,double*,int,int,
			    CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (CHM_SP(*)(const_CHM_SP,const_CHM_SP,
			     double*,double*,int,int,
			     CHM_CM))
		R_GetCCallable("Matrix", "cholmod_add");
	return fun(A, B, alpha, beta, values, sorted, Common);
    }

    CHM_DN cholmod_allocate_dense(size_t nrow, size_t ncol, size_t d,
				  int xtype, CHM_CM Common) {
	static CHM_DN(*fun)(size_t,size_t,size_t,
			    int,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (CHM_DN(*)(size_t,size_t,size_t,
			     int,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_allocate_dense");
	return fun(nrow, ncol, d, xtype, Common);
    }

    CHM_SP cholmod_allocate_sparse(size_t nrow, size_t ncol, size_t nzmax,
				   int sorted, int packed, int stype,
				   int xtype, CHM_CM Common) {
	static CHM_SP(*fun)(size_t,size_t,size_t,int,int,
			    int,int,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (CHM_SP(*)
		   (size_t,size_t,size_t,int,int,int,int,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_allocate_sparse");
	return fun(nrow,ncol,nzmax,sorted,packed,stype,xtype,Common);
    }

    CHM_TR cholmod_allocate_triplet(size_t nrow, size_t ncol, size_t nzmax,
				    int stype, int xtype, CHM_CM Common) {
	static CHM_TR(*fun)(size_t,size_t,size_t, int,int,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (CHM_TR(*)(size_t,size_t,size_t,int,int,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_allocate_triplet");
	return fun(nrow,ncol,nzmax,stype,xtype,Common);
    }

    CHM_SP cholmod_triplet_to_sparse(const CHM_TR T, int nzmax,
				     CHM_CM Common) {
	static CHM_SP(*fun)(const CHM_TR,int,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (CHM_SP(*)(const CHM_TR,int,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_triplet_to_sparse");
	return fun(T, nzmax, Common);
    }

    CHM_TR cholmod_sparse_to_triplet(const_CHM_SP A, CHM_CM Common) {
	static CHM_TR(*fun)(const_CHM_SP,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (CHM_TR(*)(const_CHM_SP,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_sparse_to_triplet");
	return fun(A, Common);
    }

    CHM_DN cholmod_sparse_to_dense(const_CHM_SP A, CHM_CM Common) {
	static CHM_DN(*fun)(const_CHM_SP,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (CHM_DN(*)(const_CHM_SP,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_sparse_to_dense");
	return fun(A, Common);
    }

    CHM_FR cholmod_analyze(const_CHM_SP A, CHM_CM Common) {
	static CHM_FR(*fun)(const_CHM_SP,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (CHM_FR(*)(const_CHM_SP,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_analyze");
	return fun(A, Common);
    }

    CHM_FR cholmod_analyze_p(const_CHM_SP A, int *Perm, int *fset,
			     size_t fsize, CHM_CM Common) {
	static CHM_FR(*fun)(const_CHM_SP,int*,int*,size_t,
			    CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (CHM_FR(*)(const_CHM_SP,int*,int*,
			     size_t,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_analyze_p");
	return fun(A, Perm, fset, fsize, Common);
    }

    CHM_SP cholmod_copy(const_CHM_SP A, int stype,
			int mode, CHM_CM Common) {
	static CHM_SP(*fun)(const_CHM_SP,int,int,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (CHM_SP(*)(const_CHM_SP,int,int,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_copy");
	return fun(A, stype, mode, Common);
    }

    CHM_DN cholmod_copy_dense(const_CHM_DN  A, CHM_CM Common) {
	static CHM_DN(*fun)(const_CHM_DN,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (CHM_DN(*)(const_CHM_DN,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_copy_dense");
	return fun(A, Common);
    }

    CHM_FR cholmod_copy_factor(const_CHM_FR L, CHM_CM Common) {
	static CHM_FR(*fun)(const_CHM_FR,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (CHM_FR(*)(const_CHM_FR,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_copy_factor");
	return fun(L, Common);
    }

    int cholmod_change_factor(int to_xtype, int to_ll, int to_super, int to_packed,
			      int to_monotonic, CHM_FR L, CHM_CM Common) {
	static int(*fun)(int,int,int,int,int,CHM_FR,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (int(*)(int,int,int,int,int,CHM_FR,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_change_factor");
	return fun(to_xtype, to_ll, to_super, to_packed, to_monotonic, L, Common);
    }

    CHM_SP cholmod_copy_sparse(const_CHM_SP A, CHM_CM Common) {
	static CHM_SP(*fun)(const_CHM_SP,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (CHM_SP(*)(const_CHM_SP,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_copy_sparse");
	return fun(A, Common);
    }

    CHM_SP cholmod_factor_to_sparse(const_CHM_FR L, CHM_CM Common) {
	static CHM_SP(*fun)(const_CHM_FR,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (CHM_SP(*)(const_CHM_FR,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_factor_to_sparse");
	return fun(L, Common);
    }

    CHM_SP cholmod_submatrix(const_CHM_SP A, int *rset, int rsize, int *cset,
			     int csize, int values, int sorted, CHM_CM Common) {
	static CHM_SP(*fun)(const_CHM_SP,int*,int,int*,int,
			    int,int,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (CHM_SP(*)(const_CHM_SP,int*,int,int*,
			     int,int,int,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_submatrix");
	return fun(A, rset, rsize, cset, csize, values, sorted, Common);
    }

    CHM_SP cholmod_dense_to_sparse(const_CHM_DN  X, int values, CHM_CM Common) {
	static CHM_SP(*fun)(const_CHM_DN,int,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (CHM_SP(*)(const_CHM_DN,int,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_dense_to_sparse");
	return fun(X, values, Common);
    }

    int cholmod_factorize(const_CHM_SP A, CHM_FR L, CHM_CM Common) {
	static int(*fun)(const_CHM_SP,CHM_FR,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (int(*)(const_CHM_SP,CHM_FR,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_factorize");
	return fun(A, L, Common);
    }

    int cholmod_factorize_p(const_CHM_SP A, double *beta, int *fset,
			    size_t fsize, CHM_FR L, CHM_CM Common) {
	static int(*fun)(const_CHM_SP,double*,int*,size_t,
			 CHM_FR,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (int(*)(const_CHM_SP,double*,int*,size_t,
			  CHM_FR,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_factorize_p");
	return fun(A, beta, fset, fsize, L, Common);
    }

    int cholmod_finish(CHM_CM Common) {
	static int(*fun)(CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (int(*)(CHM_CM))
		R_GetCCallable("Matrix", "cholmod_finish");
	return fun(Common);
    }

    int cholmod_sort(CHM_SP A, CHM_CM Common) {
	static int(*fun)(CHM_SP,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (int(*)(CHM_SP,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_sort");
	return fun(A, Common);
    }

    int cholmod_free_dense(CHM_DN  *A, CHM_CM Common) {
	static int(*fun)(CHM_DN*,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (int(*)(CHM_DN*,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_free_dense");
	return fun(A, Common);
    }

    int cholmod_free_factor(CHM_FR *L, CHM_CM Common) {
	static int(*fun)(CHM_FR*,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (int(*)(CHM_FR*,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_free_factor");
	return fun(L, Common);
    }

    int cholmod_free_sparse(CHM_SP *A, CHM_CM Common) {
	static int(*fun)(CHM_SP*,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (int(*)(CHM_SP*,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_free_sparse");
	return fun(A, Common);
    }

    int cholmod_free_triplet(cholmod_triplet **T, CHM_CM Common) {
	static int(*fun)(cholmod_triplet**,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (int(*)(cholmod_triplet**,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_free_triplet");
	return fun(T, Common);
    }

    long cholmod_nnz(const_CHM_SP A, CHM_CM Common) {
	static long(*fun)(const_CHM_SP,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (long(*)(const_CHM_SP,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_nnz");
	return fun(A, Common);
    }

    int cholmod_sdmult(const_CHM_SP A, int transpose,
		       const double *alpha, const double *beta,
		       const_CHM_DN X, CHM_DN  Y,
		       CHM_CM Common) {
	static int(*fun)(const_CHM_SP,int,const double*,
			 const double*,const_CHM_DN,
			 CHM_DN,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (int(*)(const_CHM_SP,int,const double*,
			  const double*, const_CHM_DN,
			  CHM_DN,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_sdmult");
	return fun(A, transpose, alpha, beta, X, Y, Common);
    }

    CHM_SP cholmod_ssmult(const_CHM_SP A, const_CHM_SP B,
			  int stype, int values, int sorted, CHM_CM Common) {
	static CHM_SP(*fun)(const_CHM_SP,const_CHM_SP,
			    int,int,int,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (CHM_SP(*)(const_CHM_SP,const_CHM_SP,
			     int,int,int,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_ssmult");
	return fun(A, B, stype, values, sorted, Common);
    }

    CHM_DN cholmod_solve(int sys, const_CHM_FR L, const_CHM_DN B, CHM_CM Common) {
	static CHM_DN(*fun)(int,const_CHM_FR,const_CHM_DN,
			    CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (CHM_DN(*)(int,const_CHM_FR,const_CHM_DN,
			     CHM_CM))
		R_GetCCallable("Matrix", "cholmod_solve");
	return fun(sys, L, B, Common);
    }

    CHM_SP cholmod_speye(size_t nrow, size_t ncol, int xtype, CHM_CM Common) {
	static CHM_SP(*fun)(size_t,size_t,int,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (CHM_SP(*)(size_t,size_t,int,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_speye");
	return fun(nrow, ncol, xtype, Common);
    }

    CHM_SP cholmod_spsolve(int sys, const_CHM_FR L, const_CHM_SP B, CHM_CM Common) {
	static CHM_SP(*fun)(int,const_CHM_FR,
			    const_CHM_SP, CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (CHM_SP(*)(int,const_CHM_FR,
			     const_CHM_SP, CHM_CM))
		R_GetCCallable("Matrix", "cholmod_spsolve");
	return fun(sys, L, B, Common);
    }

    int cholmod_defaults (CHM_CM Common) {
	static int(*fun)(CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (int(*)(CHM_CM))
		R_GetCCallable("Matrix", "cholmod_defaults");
	return fun(Common);
    }

    int cholmod_updown(int update, const_CHM_SP C, const_CHM_FR L, CHM_CM Common) {
	static int(*fun)(int,const_CHM_SP,const_CHM_FR,
			 CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (int(*)(int,const_CHM_SP,const_CHM_FR,
			  CHM_CM))
		R_GetCCallable("Matrix", "cholmod_updown");
	return fun(update, C, L, Common);
    }

    CHM_SP cholmod_transpose(const_CHM_SP A, int values, CHM_CM Common) {
	static CHM_SP(*fun)(const_CHM_SP,int,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (CHM_SP(*)(const_CHM_SP,int,CHM_CM))
		R_GetCCallable("Matrix", "cholmod_transpose");
	return fun(A, values, Common);
    }

    CHM_SP cholmod_vertcat(const_CHM_SP A, const_CHM_SP B, int values, CHM_CM Common) {
	static CHM_SP(*fun)(const_CHM_SP,const_CHM_SP,int,CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (CHM_SP(*)(const_CHM_SP,const_CHM_SP, int, CHM_CM))
		R_GetCCallable("Matrix", "cholmod_vertcat");
	return fun(A, B, values, Common);
    }

    int cholmod_scale(const_CHM_DN S, int scale, CHM_SP A, CHM_CM Common) {
	static int(*fun)(const_CHM_DN,int,CHM_SP, CHM_CM) = NULL;
	if (fun == NULL)
	    fun = (int(*)(const_CHM_DN,int,CHM_SP, CHM_CM))
		R_GetCCallable("Matrix", "cholmod_scale");
	return fun(S, scale, A, Common);
    }
}

namespace CHM {
    dsCMatrix::dsCMatrix(S4& A)
	: d_Dim(A.slot("Dim")),
	  d_upper(CharacterVector(A.slot("uplo"))[0] == "U"),
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
    
    void dsCMatrix::update_factors() {
	cholmod_common* cm = new cholmod_common;
	if (!cholmod_start(cm)) stop("call to cholmod_start failed");
	cm->final_asis = 0; 
	cm->final_monotonic = 1; cm->final_pack = 1;
	for (int i = 0; i < d_factors.size(); ++i) {
	    S4 faci = d_factors[i];
	    if (!faci.is("CHMfactor")) continue;
	    if (faci.is("dCHMsimpl")) {
		CHM::dCHMsimpl ff(faci);
		cm->final_ll = ff.is_ll();
		if (!cholmod_factorize(spp(), ff.frp(), cm))
		    stop("failure in cholmod_factorize");
	    } else if (faci.is("dCHMsuper")) {
		CHM::dCHMsuper ff(faci);
		if (!cholmod_factorize(spp(), ff.frp(), cm))
		    stop("failure in cholmod_factorize");
	    } else Rf_warning("A@factors contains a CHMfactor that is neither \"dCHMsimpl\" nor \"dCHMsuper\"");
	}
	delete cm;
    }
   
    NumericMatrix dsCMatrix::solve(const NumericMatrix& b, int which) const {
	if (!n_factors()) stop("dsCMatrix::solve requires n_factors() > 0");
	S4 f0(d_factors[0]);
	if (!f0.is("CHMfactor")) stop("dsCMatrix::solve for CHMfactor only");
// Should use virtual methods here
	if (f0.is("dCHMsimpl")) {
	    const CHM::dCHMsimpl  ff(f0);
	    return ff.solve(b, which);
	}
	// if (f0.is("dCHMsuper")) {
	//     const CHM::dCHMsuper ff(f0);
	//     return ff.solve(b, which);
	// }
	stop("Unrecognized numeric Cholesky factorization");
	return NumericMatrix(1,1); // -Wall
    }

    NumericVector dsCMatrix::Ldiag() { // can't make this const unless S4.is is const
	for (int i = 0; i < n_factors(); ++i) {
	    S4 faci = d_factors[i];
	    if (!faci.is("CHMfactor")) continue;
	    if (faci.is("dCHMsimpl")) return CHM::dCHMsimpl(faci).Ldiag();
	    if (faci.is("dCHMsuper")) return CHM::dCHMsuper(faci).Ldiag();
	}
	stop("At least one CHMfactor must be defined for Ldiag");
	return NumericVector(0); // -Wall
    }
	
    CHMfactor::CHMfactor(S4& L)
	: d_colcount(L.slot("colcount")),
	  d_perm(L.slot("perm")),
	  d_type(L.slot("type")),
	  d_Dim(L.slot("Dim")),
	  d_fr(new cholmod_factor) {
	d_fr->n = d_Dim[0];
	d_fr->Perm = (void*)&d_perm[0];
	d_fr->ColCount = (void*)&d_colcount[0];
	d_fr->itype = CHOLMOD_INT;
	d_fr->xtype = CHOLMOD_PATTERN;
	d_fr->dtype = CHOLMOD_DOUBLE;
	d_fr->ordering = d_type[0];	/* unravel the type */
	d_fr->is_ll = (d_type[1] ? 1 : 0);
	d_fr->is_super = (d_type[2] ? 1 : 0);
	d_fr->is_monotonic = (d_type[3] ? 1 : 0);
    }
    
    CHMfactor::~CHMfactor() {delete d_fr;}
    
    CHMsimpl::CHMsimpl(S4& L)
	: CHMfactor(L),
	  d_p(L.slot("p")),
	  d_i(L.slot("i")),
	  d_nz(L.slot("nz")),
	  d_nxt(L.slot("nxt")),
	  d_prv(L.slot("prv")) {
	if (is_super()) stop("attempt to create a CHMsimpl object from a supernodal decomposition");
	d_fr->p = (void*)&d_p[0];
	d_fr->i = (void*)&d_i[0];
	d_fr->nz = (void*)&d_nz[0];
	d_fr->next = (void*)&d_nxt[0];
	d_fr->prev = (void*)&d_prv[0];
    }
    
    dCHMsimpl::dCHMsimpl(S4 &L)
	: CHMsimpl(L),
	  d_x(L.slot("x")) {
	d_fr->xtype = CHOLMOD_REAL;
	d_fr->x = (void*)&d_x[0];
    }
    
    NumericMatrix dCHMsimpl::solve(const Rcpp::NumericMatrix& b, int which) const {
	cholmod_common* cm = new cholmod_common;
	if (!cholmod_start(cm)) stop("call to cholmod_start failed");
	cholmod_dense* dn = new cholmod_dense;
	dn->nrow = b.nrow(); dn->ncol = b.ncol(); dn->nzmax = b.nrow() * b.ncol();
	dn->d = b.nrow(); dn->x = (void*)&b[0]; dn->xtype = CHOLMOD_REAL;
	dn->dtype = CHOLMOD_DOUBLE;
	cholmod_dense* xp = cholmod_solve(which, d_fr, dn, cm);
	double* dp = (double*) xp->x;
	NumericMatrix x = clone<NumericMatrix>(b);
	std::copy(dp, dp + b.nrow() * b.ncol(), x.begin());
	if (!cholmod_free_dense(&xp, cm)) stop("failure in cholmod_free_dense");
	delete cm;
	return x;
    }

    NumericVector dCHMsimpl::Ldiag() const {
	NumericVector ans(ncol());
	for (int i = 0; i < ncol(); ++i) ans[i] = d_x[d_p[i]];
	if (!d_fr->is_ll) return sqrt(ans);
	return ans;
    }

    CHMsuper::CHMsuper(S4 &L)
	: CHMfactor(L),
	  d_super(L.slot("super")),
	  d_pi(L.slot("pi")),
	  d_px(L.slot("px")),
	  d_s(L.slot("s")) {
	if (!is_super()) stop("attempt to create a CHMsuper object from a simplicial decomposition");
	d_fr->xsize = d_perm.size();
	d_fr->maxcsize = d_type[4];
	d_fr->maxesize = d_type[5];
	d_fr->nsuper = d_super.size() - 1;
	d_fr->super = (void*)&d_super[0];
	d_fr->pi = (void*)&d_pi[0];
	d_fr->px = (void*)&d_px[0];
	d_fr->ssize = d_s.size();
	d_fr->s = (void*)&d_s[0];
    }
    
    dCHMsuper::dCHMsuper(S4 &L)
	: CHMsuper(L),
	  d_x(L.slot("x")) {
	d_fr->xtype = CHOLMOD_REAL;
	d_fr->x = (void*)&d_x[0];
    }

    NumericMatrix dCHMsuper::solve(const Rcpp::NumericMatrix& b, int which) const {
	cholmod_common* cm = new cholmod_common;
	if (!cholmod_start(cm)) stop("call to cholmod_start failed");
	cholmod_dense* dn = new cholmod_dense;
	dn->nrow = b.nrow(); dn->ncol = b.ncol(); dn->nzmax = b.nrow() * b.ncol();
	dn->d = b.nrow(); dn->x = (void*)&b[0]; dn->xtype = CHOLMOD_REAL;
	dn->dtype = CHOLMOD_DOUBLE;
	cholmod_dense* xp = cholmod_solve(which, d_fr, dn, cm);
	double* dp = (double*) xp->x;
	NumericMatrix x = clone<NumericMatrix>(b);
	std::copy(dp, dp + b.nrow() * b.ncol(), x.begin());
	if (!cholmod_free_dense(&xp, cm)) stop("failure in cholmod_free_dense");
	delete cm;
	return x;
    }

    NumericVector dCHMsuper::Ldiag() const {
	stop("Ldiag not yet written for dCHMsuper");
	NumericVector ans(ncol());
	// for (int i = 0; i < ncol(); ++i) ans[i] = d_x[d_p[i]];
	// if (!d_fr->is_ll)
	//     std::transform(ans.begin(), ans.end(), ans.begin(), std::sqrt);
	return ans;
    }

}
