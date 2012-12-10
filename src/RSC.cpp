#include <Rcpp.h>

using namespace Rcpp;

static void apply_lambda(const NumericVector &theta,
			 const NumericVector &lower,
			 std::vector<double> &dest) {
    int dpos(-1);
    double *rr(0);		// initialize to NULL pointer
    for (int k = 0; k < theta.size(); ++k) {
	if (lower[k] == 0.) {	// diagonal element of a factor
	    dest[++dpos] *= theta[k];
	    rr = &dest[k];
	}
	else dest[dpos] += theta[k] * *++rr;
    }
}
//' Update for the penalized least squares problem
//' 
//' @param rv the matrix of row indices for the regular sparse column representation of ZtXt
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
	apply_lambda(theta, lower, w);
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
