#include <Rcpp.h>

using namespace Rcpp;
//' Update for the penalized least squares problem
//' 
//' @param rv the matrix of row indices for the regular sparse column representation of ZtXt
//' @param xv the non-zero values in ZtXt
//' @param theta vector of relative standard deviations
//' @param m the sparse matrix to be updated.  Must have the correct pattern.
//' @examples
//' rv <- rbind(as.integer(Dyestuff$Batch) - 1L, 6L)
//' xv <- matrix(1, nrow=nrow(rv), ncol=ncol(rv)
//' A <- sparseMatrix(i=as.integer(c(1,1,2,2,3,3,4,4,5,5,6,6,7)),
//'                   j=as.integer(c(1,7,2,7,3,7,4,7,5,7,6,7,7)),
//'                   x=rep.int(1.,13L), symmetric=TRUE,check=FALSE)
//' RSCupdate(rv, xv, 0.831, A)
//' A
// [[Rcpp::export]]
void RSCupdate(const IntegerMatrix rv, const NumericMatrix xv,
	       const NumericVector theta, S4 m) {
    int k(theta.size()), kpp(rv.nrow()), n(rv.ncol());
    int kppm1(kpp - 1);
    IntegerVector dims = m.slot("Dim");
    IntegerVector rowind = m.slot("i"), colptr = m.slot("p");
    NumericVector nzvals = m.slot("x");
    int nz(nzvals.size());
    if (dims[0] != dims[1])
	Rcpp::stop("m must be a square matrix");
    if (xv.nrow() != kpp || xv.ncol() != n)
	Rcpp::stop("sizes of rv and xv must agree");
    for (int i = 0; i < nz; ++i) nzvals[i] = 0.; // initialize m to zeros
    
    for (int j = 0; j < n; ++j)	{ // iterate over columns of ZtXt
	for (int i = kppm1; i >= 0; --i) { // reverse order of rows of the RSC rep
	    int ii(rv(i, j));		// row in ZtXt and column in m
	    int cpi(colptr[ii]),	// first rowind and nzvals index for this column
		ll(colptr[ii + 1] - 1); // should be location of diagonal
	    double rm(xv(i, j));	// non-zero value in ZtXt
	    if (rowind[ll] != ii) Rcpp::stop("u is not upper triangular");
	    if (i < k) {	 // if in the Z part
		rm *= theta[i];	 // scale by theta[i]
		nzvals[ll] = 1.; // initialize diagonal element to 1.
	    }
	    nzvals[ll] += rm * rm;
	    --ll;		// move up the column of m to next non-zero
	    for (int l = i - 1; l >= 0 && ll >= cpi; --l) { // off-diagonals in m
		int ii1(rv(l,j)); // row in ZtXt
		while (rowind[ll] > ii1) --ll;
		if (rowind[ll] != ii1) Rcpp::stop("Pattern mismatch");
		nzvals[ll] += (rm * xv(l,j) * (l < k ? theta[l] : 1.));
	    }
	}
    }
}
