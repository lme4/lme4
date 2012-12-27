// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// CHM.h: sparse matrix classes paralleling those in the Matrix package
//
// Copyright (C) 2012-2013 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#ifndef LME4_CHM_H
#define LME4_CHM_H

#include <Rcpp.h>
#include "cholmod.h"

extern "C" {
    int cholmod_start(CHM_CM);
    CHM_SP cholmod_allocate_sparse(size_t nrow, size_t ncol,
                                   size_t nzmax, int sorted,
                                   int packed, int stype, int xtype,
                                   CHM_CM);
    int cholmod_free_factor(CHM_FR*, CHM_CM);
    int cholmod_free_dense(CHM_DN*, CHM_CM);
    int cholmod_free_sparse(CHM_SP*, CHM_CM);
    int cholmod_free_triplet(CHM_TR*, CHM_CM);
    long cholmod_nnz(const_CHM_SP, CHM_CM);
    CHM_SP cholmod_speye(size_t nrow, size_t ncol, int xtype, CHM_CM);
    CHM_SP cholmod_transpose(const_CHM_SP, int values, CHM_CM);
    int cholmod_sort(CHM_SP A, CHM_CM);
    CHM_SP cholmod_vertcat(const_CHM_SP, const_CHM_SP, int values, CHM_CM);
    CHM_SP cholmod_copy(const_CHM_SP, int stype, int mode, CHM_CM);
    CHM_SP cholmod_add(const_CHM_SP, const_CHM_SP, double alpha [2], double beta [2],
                       int values, int sorted, CHM_CM);
    int cholmod_finish(CHM_CM);
    CHM_DN cholmod_solve(int, const_CHM_FR, const_CHM_DN, CHM_CM);
    CHM_SP cholmod_spsolve(int, const_CHM_FR, const_CHM_SP, CHM_CM);
    int cholmod_sdmult(const_CHM_SP, int, const double*, const double*,
                       const_CHM_DN, CHM_DN Y, CHM_CM);
    CHM_SP cholmod_ssmult(const_CHM_SP, const_CHM_SP, int, int, int,
                          CHM_CM);
    int cholmod_factorize(const_CHM_SP, CHM_FR L, CHM_CM);
    int cholmod_factorize_p(const_CHM_SP, double *beta, int *fset,
                            size_t fsize, CHM_FR L, CHM_CM);
    CHM_SP cholmod_copy_sparse(const_CHM_SP, CHM_CM);
    CHM_DN cholmod_copy_dense(const_CHM_DN, CHM_CM);
    CHM_SP cholmod_aat(const_CHM_SP, int *fset, size_t fsize, int mode,
                       CHM_CM);
    int cholmod_band_inplace(CHM_SP A, int k1, int k2, int mode, CHM_CM);
    CHM_SP cholmod_add(const_CHM_SP, const_CHM_SP, double alpha[2], double beta[2],
                       int values, int sorted, CHM_CM);
    CHM_DN cholmod_allocate_dense(size_t nrow, size_t ncol, size_t d,
                                  int xtype, CHM_CM);
    CHM_FR cholmod_analyze(const_CHM_SP, CHM_CM);
    CHM_FR cholmod_analyze_p(const_CHM_SP, int *Perm, int *fset,
                             size_t fsize, CHM_CM);
    int cholmod_change_factor(int to_xtype, int to_ll, int to_super,
                              int to_packed, int to_monotonic,
                              CHM_FR L, CHM_CM);
    CHM_FR cholmod_copy_factor(const_CHM_FR, CHM_CM);
    CHM_SP cholmod_factor_to_sparse(const_CHM_FR, CHM_CM);
    CHM_SP cholmod_dense_to_sparse(const_CHM_DN, int values, CHM_CM);
    int cholmod_defaults(CHM_CM);
    CHM_SP cholmod_triplet_to_sparse(const CHM_TR, int nzmax, CHM_CM);
    CHM_SP cholmod_submatrix(const_CHM_SP, int *rset, int rsize, int *cset,
                             int csize, int values, int sorted,
                             CHM_CM);
    CHM_TR cholmod_sparse_to_triplet(const_CHM_SP, CHM_CM);
    CHM_DN cholmod_sparse_to_dense(const_CHM_SP, CHM_CM);
    CHM_TR cholmod_allocate_triplet (size_t nrow, size_t ncol, size_t nzmax,
                                     int stype, int xtype, CHM_CM);
    
    int cholmod_scale(const_CHM_DN, int scale, CHM_SP, CHM_CM);
}


namespace CHM {
// add classes CsparseMatrix, gCMatrix, sCMatrix, tCMatrix, dgCMatrix, dsCMatrix, dtCMatrix
    class dsCMatrix {
    protected:
        const Rcpp::IntegerVector d_Dim;
        const bool                d_upper;
        Rcpp::IntegerVector       d_colptr;
        Rcpp::IntegerVector       d_rowval;
        const Rcpp::List          d_factors;
        Rcpp::NumericVector       d_nzval;
        cholmod_sparse*           d_sp;
    public:
        dsCMatrix(Rcpp::S4&);
        ~dsCMatrix();
        // extractor methods
        const Rcpp::IntegerVector&    Dim() const {return d_Dim;}
        int                          nrow() const {return d_Dim[0];}
        int                          ncol() const {return d_Dim[1];}
        int                           nnz() const {return d_rowval.size();}
        bool                        upper() const {return d_upper;}
        const Rcpp::IntegerVector& colptr() const {return d_colptr;}
        Rcpp::IntegerVector&       colptr()       {return d_colptr;}
        const Rcpp::IntegerVector& rowval() const {return d_rowval;}
        Rcpp::IntegerVector&       rowval()       {return d_rowval;}
        const Rcpp::NumericVector&  nzval() const {return d_nzval;}
        Rcpp::NumericVector&        nzval()       {return d_nzval;}
        const Rcpp::List&         factors() const {return d_factors;}
        const cholmod_sparse*         spp() const {return d_sp;}
        cholmod_sparse*               spp()       {return d_sp;}
        int                     n_factors() const {return d_factors.size();}
        Rcpp::NumericVector         Ldiag();

        void               update_factors();

        Rcpp::NumericMatrix solve(const Rcpp::NumericMatrix&, int) const;
    };

    class CHMfactor {
    protected:
        const Rcpp::IntegerVector d_colcount;
        const Rcpp::IntegerVector d_perm;
        const Rcpp::IntegerVector d_type;
        const Rcpp::IntegerVector d_Dim;
        cholmod_factor*           d_fr;
    public:
        CHMfactor(Rcpp::S4&);
        ~CHMfactor();
        const Rcpp::IntegerVector& colcount() const {return d_colcount;}
        const Rcpp::IntegerVector&     perm() const {return d_perm;}
        const Rcpp::IntegerVector&     type() const {return d_type;}
        int                            nrow() const {return d_Dim[0];}
        int                            ncol() const {return d_Dim[1];}
        bool                          is_ll() const {return bool(d_fr->is_ll);}
        bool                       is_super() const {return bool(d_fr->is_super);}
        bool                   is_monotonic() const {return bool(d_fr->is_monotonic);}
        cholmod_factor*                 frp()       {return d_fr;}
        const cholmod_factor*           frp() const {return d_fr;}
    };
    
    class CHMsimpl : public CHMfactor {
    protected:
        const Rcpp::IntegerVector d_p;
        const Rcpp::IntegerVector d_i;
        const Rcpp::IntegerVector d_nz;
        const Rcpp::IntegerVector d_nxt;
        const Rcpp::IntegerVector d_prv;
    public:
        CHMsimpl(Rcpp::S4&);
        const Rcpp::IntegerVector&   p() const {return d_p;}
        const Rcpp::IntegerVector&   i() const {return d_i;}
        const Rcpp::IntegerVector&  nz() const {return d_nz;}
        const Rcpp::IntegerVector& nxt() const {return d_nxt;}
        const Rcpp::IntegerVector& prv() const {return d_prv;}
    };
    
    class dCHMsimpl : public CHMsimpl {
    protected:
        Rcpp::NumericVector   d_x;
    public:
        dCHMsimpl(Rcpp::S4&);
        const Rcpp::NumericVector& x() const {return d_x;}
        Rcpp::NumericVector&       x()       {return d_x;}
        Rcpp::NumericVector    Ldiag() const;

        Rcpp::NumericMatrix solve(const Rcpp::NumericMatrix&, int) const;
    };
    
    class CHMsuper : public CHMfactor {
    protected:
        const Rcpp::IntegerVector d_super;
        const Rcpp::IntegerVector d_pi;
        const Rcpp::IntegerVector d_px;
        const Rcpp::IntegerVector d_s;
    public:
        CHMsuper(Rcpp::S4&);
        const Rcpp::IntegerVector& super() const {return d_super;}
        const Rcpp::IntegerVector&    pi() const {return d_pi;}
        const Rcpp::IntegerVector&    px() const {return d_px;}
        const Rcpp::IntegerVector&     s() const {return d_s;}
    };
    
    class dCHMsuper : public CHMsuper {
    protected:
        Rcpp::NumericVector d_x;
    public:
        dCHMsuper(Rcpp::S4&);
        const Rcpp::NumericVector& x() const {return d_x;}
        Rcpp::NumericVector&       x()       {return d_x;}
        Rcpp::NumericVector    Ldiag() const;

        Rcpp::NumericMatrix solve(const Rcpp::NumericMatrix&, int) const;
    };
    
}
#endif // LME4_CHM_H
    
