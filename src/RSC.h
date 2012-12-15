#include <Rcpp.h>
#include <Matrix.h>

using namespace Rcpp;		// probably should remove this - it's bad form
class dsCMatrix {
protected:
    const IntegerVector d_Dim;
    const bool d_upper;
    IntegerVector d_colptr;
    IntegerVector d_rowval;
    const List d_factors;
    NumericVector d_nzval;
public:
    dsCMatrix(S4&);
				// extractor methods
    const IntegerVector&     Dim() const {return d_Dim;}
    int                     nrow() const {return d_Dim[0];}
    int                     ncol() const {return d_Dim[1];}
    int                      nnz() const {return d_rowval.size();}
    bool                   upper() const {return d_upper;}
    const IntegerVector&  colptr() const {return d_colptr;}
    IntegerVector&        colptr()       {return d_colptr;}
    const IntegerVector&  rowval() const {return d_rowval;}
    IntegerVector&        rowval()       {return d_rowval;}
    const NumericVector&   nzval() const {return d_nzval;}
    NumericVector&         nzval()       {return d_nzval;}
    int                n_factors() const {return d_factors.size();}
				// return cholmod struct pointers
    void          update_factors() {};
};

class CHM_SP_wrap {
protected:
    CHM_SP   d_sp;
    CHM_CM   d_cm;
    bool   d_mine;
public:
    CHM_SP_wrap(dsCMatrix&);
    CHM_SP_wrap(CHM_SP, CHM_CM);
    ~CHM_SP_wrap();
    CHM_SP       as_CHM_SP()       {return d_sp;}
    const CHM_SP as_CHM_SP() const {return d_sp;}
};

class CHMfactor {
protected:
    const IntegerVector d_colcount;
    const IntegerVector d_perm;
    const IntegerVector d_type;
public:
    CHMfactor(S4&);
    const IntegerVector &colcount() const {return d_colcount;}
    const IntegerVector     &perm() const {return d_perm;}
    const IntegerVector     &type() const {return d_type;}
};

class CHMsimpl : public CHMfactor {
protected:
    const IntegerVector d_p;
    const IntegerVector d_i;
    const IntegerVector d_nz;
    const IntegerVector d_nxt;
    const IntegerVector d_prv;
public:
    CHMsimpl(S4 &L);
    const IntegerVector   &p() const {return d_p;}
    const IntegerVector   &i() const {return d_i;}
    const IntegerVector  &nz() const {return d_nz;}
    const IntegerVector &nxt() const {return d_nxt;}
    const IntegerVector &prv() const {return d_prv;}
};

class dCHMsimpl : public CHMsimpl {
protected:
    NumericVector   d_x;
public:
    dCHMsimpl(S4 &L);
    const NumericVector& x() const {return d_x;}
    NumericVector&       x()       {return d_x;}
};

class CHMsuper : public CHMfactor {
protected:
    const IntegerVector d_super;
    const IntegerVector d_pi;
    const IntegerVector d_px;
    const IntegerVector d_s;
public:
    CHMsuper(S4 &L);
    const IntegerVector& super() const {return d_super;}
    const IntegerVector&    pi() const {return d_pi;}
    const IntegerVector&    px() const {return d_px;}
    const IntegerVector&     s() const {return d_s;}
};

class dCHMsuper : public CHMsuper {
protected:
    NumericVector d_x;
public:
    dCHMsuper(S4 &L);
    const NumericVector& x() const {return d_x;}
    NumericVector&       x()       {return d_x;}
};
    
class RSC { /**< const parts of mixed-effects predictor in regular sparse column format */
protected:
    const IntegerMatrix rv;	/**< rowvals matrix for Zt */
    const NumericMatrix xv;	/**< xvals matrix for ZtXt */
    const NumericVector lower; /**< lower bounds for covariance parameters */
    const int k;      /**< number of random effects per observation */
    const int kpp;    /**< number of rows in xv = k + p */
    const int n;      /**< number of observations */
    const int p;      /**< number of fixed-effects coefficients */
    const int q;      /**< total number of random effects */
public:
    RSC(const SEXP, const SEXP, const SEXP);
    NumericVector &apply_lambda(const NumericVector&, NumericVector&) const;
    void update_A(const NumericVector&, const NumericVector&, S4 &AA, NumericVector&)
	const;
};
