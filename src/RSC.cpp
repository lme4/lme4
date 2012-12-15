#include "RSC.h"

extern "C" {

    void R_cholmod_error(int status, const char *file, int line,
			 const char *message) {
	if(status < 0) {
	    Rcout << message << " at file " << file
		  << ", line " << line << std::endl;
	    stop("Cholmod error");
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

dsCMatrix::dsCMatrix(S4 &A)
    : d_Dim(A.slot("Dim")),
      d_upper(CharacterVector(A.slot("Dim"))[0] == "U"),
      d_colptr(A.slot("p")),
      d_rowval(A.slot("i")),
      d_factors(A.slot("factors")),
      d_nzval(A.slot("x")) {}
// still need to implement dsCMatrix::update_factors;

CHM_SP_wrap::CHM_SP_wrap(dsCMatrix& A)
    : d_sp(new cholmod_sparse),
      d_cm(new cholmod_common),
      d_mine(true) {
    d_sp = new cholmod_sparse;
    d_sp->nrow = A.nrow();
    d_sp->ncol = A.ncol();
    d_sp->nzmax = A.nnz();
    d_sp->p = (void*)&A.colptr()[0];
    d_sp->i = (void*)&A.rowval()[0];
    d_sp->nz = NULL;
    d_sp->x = (void*)&A.nzval()[0];
    d_sp->z = NULL;
    d_sp->stype = A.upper() ? 1 : -1;
    d_sp->itype = CHOLMOD_INT;
    d_sp->xtype = CHOLMOD_REAL;
    d_sp->dtype = CHOLMOD_DOUBLE;
    d_sp->sorted = 1;
    d_sp->packed = 1;
    if (!cholmod_start(d_cm)) stop("failure in cholmod_start");
}

CHM_SP_wrap::CHM_SP_wrap(CHM_SP A, CHM_CM c)
    : d_sp(A),
      d_cm(c),
      d_mine(false) {}

CHM_SP_wrap::~CHM_SP_wrap() {
    if (d_mine) {
	delete d_sp;
	delete d_cm;
    } else cholmod_free_sparse(&d_sp, d_cm);
}

	
CHMfactor::CHMfactor(S4 &L)
    : d_colcount(L.slot("colcount")),
      d_perm(L.slot("perm")),
      d_type(L.slot("type")) {}

CHMsimpl::CHMsimpl(S4 &L)
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

RSC::RSC(const SEXP rvSEXP, const SEXP xvSEXP, const SEXP lowerSEXP)
    : rv(rvSEXP),
      xv(xvSEXP),
      lower(lowerSEXP),
      k(rv.nrow()),
      kpp(xv.nrow()),
      n(xv.ncol()),
      p(kpp - k),
      q((*std::max_element(rv.begin(), rv.end())) + 1) {
    Rcout << "k = " << k << ", kpp = " << kpp << ", n = " << n
	  << ", p = " << p << ", q = " << q << std::endl;
    if (rv.ncol() != n || p < 0) stop("dimension mismatch of rv and xv");
    if (k != std::count_if(lower.begin(), lower.end(), ::R_finite))
	stop("dimension mismatch of rv and lower");
    if (*std::min_element(rv.begin(), rv.end()) != 0)
	stop("minimum row index must be 0");
}

NumericVector &RSC::apply_lambda(const NumericVector &theta,
				 NumericVector &dest) const {
    int dpos(-1);
    double *rr(0);		// initialize to NULL pointer
    for (int k = 0; k < kpp; ++k) {
	if (lower[k] == 0.) { // diagonal element of a factor
	    dest[++dpos] *= theta[k];
	    rr = &dest[k];
	}
	else dest[dpos] += theta[k] * *++rr;
    }
    return dest;
}

void RSC::update_A(const NumericVector &theta, const NumericVector &resid,
		  S4 &AA, NumericVector &ubeta) const {
    if (theta.size() != lower.size())
	stop("Dimension mismatch of theta and lower");
    if (resid.size() != n) stop("Dimension of resid should be n");
    if (ubeta.size() != q + p) stop("Dimension of ubeta should be q + p");
    if (!AA.is("dsCMatrix")) stop("A must be a \"dsCMatrix\" object");
    dsCMatrix A(AA);
    if (A.nrow() != q + p) stop("size of A must be q + p");
    const IntegerVector &rowval(A.rowval()), &colptr(A.colptr());
    NumericVector& nzval(A.nzval());
    NumericVector w(kpp);
				// initializations
    std::fill(nzval.begin(), nzval.end(), double(0)); // zero the contents of A
    std::fill(ubeta.begin(), ubeta.end(), double(0)); // and ubeta
    for (int i = 0; i < q; ++i) {  // initialize Z part of A to the identity
	int ll(colptr[i + 1] - 1); // index of last element in column i
	int ii(rowval[ll]);	   // should be the diagonal element
	if (ii != i) stop("A is not stored as the upper triangle");
	nzval[ii] = 1.;
    }
    for (int j = 0; j < n; ++j)	{ // iterate over columns of ZtXt
	double rj(resid[j]);
	std::copy(&(xv(0,j)), &(xv(0,j)) + kpp, w.begin()); // copy j'th column
	apply_lambda(theta, w);
	for (int i = 0; i < kpp; ++i) ubeta[rv(i,j)] += rj * w[i];
				// scan up the j'th column of ZtXt, which makes
				// it easier to evaluate the upper triangle of A
	for (int i = kpp - 1; i >= 0; --i) {
	    int ii(rv(i, j));	// row of ZtXt (column of m) 
	    int cpi(colptr[ii]), ll(colptr[ii + 1] - 1); // ll should be location of diagonal
	    if (rowval[ll] != ii) stop("u is not upper triangular");
	    nzval[ll] += w[i] * w[i];
	    for (int l = i - 1; l >= 0 && ll >= cpi; --l) { // off-diagonals in m
		int ii1(rv(l,j)); // row in ZtXt
		while (rowval[ll] > ii1) --ll; // up to the desired row
		if (rowval[ll] != ii1) stop("Pattern mismatch");
		nzval[ll] += w[i] * w[l];
	    }
	}
    }
}


//' Update for the penalized least squares problem
//' 
//' @param rv the matrix of row indices for Zt as a regular sparse column matrix
//' @param xv the non-zero values in ZtXt
//' @param theta the covariace parameter vector
//' @param lower lower bounds on the covariace parameter vector
//' @param resid current residual
//' @param m the sparse matrix to be updated.  Must have the correct pattern.
//' @param ZtXty the product of the updated LambdatZt and Xt with resid.
//' @examples
//' fm1 <- lmer(Yield ~ 1|Batch, Dyestuff, REML=FALSE)
//' rv <- with(Dyestuff, rbind(as.integer(Batch) - 1L, 6L))
//' xv <- matrix(1, nrow=nrow(rv), ncol=ncol(rv))
//' A  <- tcrossprod(sparseMatrix(i = as.vector(rv),
//'                               j = as.vector(col(rv) - 1L),
//'                               x = as.vector(xv), index1=FALSE))
//' ZtXty <- numeric(ncol(A))
//' with(Dyestuff, RSCupdate(rv, xv, getME(fm1,"theta"), 0., Yield, A, ZtXty))
//' L <- Cholesky(A, perm=FALSE, LDL=FALSE)
//' vv <- solve(L, ZtXty, system="A")
//' all.equal(vv[1:6], getME(fm1, "u"))
//' all.equal(vv[7], getME(fm1, "beta"))
//' LL <- as(L, "Matrix")
//' ## the transpose below is because a diagonal dtCMatrix is declared upper triangular
//' all.equal(t(LL[1:6,1:6]), as(getME(fm1, "L"), "Matrix"))
//' all.equal(LL[7,7], as.vector(getME(fm1, "RX")))
//' all.equal(LL[7,1:6], as.vector(getME(fm1, "RZX")))
//' @export
// [[Rcpp::export]]
void RSCupdate(const IntegerMatrix rv, const NumericMatrix xv,
	       const NumericVector theta, const NumericVector lower,
	       const NumericVector resid, S4 m, NumericVector ZtXty) {
    int k = std::count_if(lower.begin(), lower.end(), ::R_finite);
    int kpp(rv.nrow());		// k + p where p = length(fixef)
    int n(rv.ncol());		// number of observations
    int p(kpp - k);		// size of the fixed-effects vector
    int qpp(ZtXty.size());	// q + p (see below for def'n of q)
    int q(qpp - p);		// total number of random effects
    std::vector<double> w(kpp); // workspace for a copy of a column of xv

    if (xv.nrow() != kpp || xv.ncol() != n) stop("sizes of rv and xv must agree");
    if (resid.size() != n) stop("length(resid) != ncol(xv)");

    IntegerVector dims = m.slot("Dim"), rowind = m.slot("i"), colptr = m.slot("p");
    
    if (dims[0] != dims[1]) stop("m must be a square matrix");
    if (ZtXty.size() != dims[1]) stop("length(ZtXty) != ncol(m)");

    NumericVector nzvals = m.slot("x");
				// initializations
    std::fill(nzvals.begin(), nzvals.end(), double(0)); // zero the contents of m
    std::fill(ZtXty.begin(), ZtXty.end(), double(0));	// and ZtXty
    for (int i = 0; i < q; ++i) { // initialize Z part of m to the identity
	int ll(colptr[i + 1] - 1); // index of last element in column i
	int ii(rowind[ll]);	   // should be the diagonal element
	if (ii != i) stop("m is not stored as the upper triangle");
	nzvals[ii] = 1.; 
    }
    for (int j = 0; j < n; ++j)	{ // iterate over columns of ZtXt
	double rj(resid[j]);
	std::copy(&(xv(0,j)), &(xv(0,j)) + kpp, w.begin()); // copy j'th column
//	apply_lambda(theta, lower, w);
	for (int i = 0; i < kpp; ++i) ZtXty[rv(i,j)] += rj * w[i];
				// scan up the j'th column of ZtXt, which makes
				// it easier to evaluate the upper triangle of A
	for (int i = kpp - 1; i >= 0; --i) {
	    int ii(rv(i, j));	// row of ZtXt (column of m) 
	    int cpi(colptr[ii]), ll(colptr[ii + 1] - 1); // ll should be location of diagonal
	    if (rowind[ll] != ii) stop("u is not upper triangular");
	    nzvals[ll] += w[i] * w[i];
	    for (int l = i - 1; l >= 0 && ll >= cpi; --l) { // off-diagonals in m
		int ii1(rv(l,j)); // row in ZtXt
		while (rowind[ll] > ii1) --ll; // up to the desired row
		if (rowind[ll] != ii1) Rcpp::stop("Pattern mismatch");
		nzvals[ll] += w[i] * w[l];
	    }
	}
    }
}
