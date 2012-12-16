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
#include <Matrix.h>

// add classes CsparseMatrix, 
// gCMatrix, scMatrix, tCMatrix, dgCMatrix, dsCMatrix, dtCMatrix
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
    const cholmod_sparse*         spp() const {return d_sp;}
    cholmod_sparse*               spp()       {return d_sp;}
    int                     n_factors() const {return d_factors.size();}
				// return cholmod struct pointers
    void               update_factors() {};
};

class CHMfactor {
protected:
    const Rcpp::IntegerVector d_colcount;
    const Rcpp::IntegerVector d_perm;
    const Rcpp::IntegerVector d_type;
    cholmod_factor*           d_fr;
public:
    CHMfactor(Rcpp::S4&);
    ~CHMfactor();
    const Rcpp::IntegerVector& colcount() const {return d_colcount;}
    const Rcpp::IntegerVector&     perm() const {return d_perm;}
    const Rcpp::IntegerVector&     type() const {return d_type;}
    cholmod_factor*                  fr()       {return d_fr;}
    const cholmod_factor*            fr() const {return d_fr;}
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
};

#endif // LME4_CHM_H
