// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

typedef Eigen::MatrixXi                      iMat;
typedef Eigen::MatrixXd                       Mat;
typedef Eigen::VectorXd                       Vec;
typedef Eigen::Map<iMat>                    MiMat;
typedef Eigen::Map<Mat>                      MMat;
typedef Eigen::Map<Vec>                      MVec;
typedef Eigen::MappedSparseMatrix<double>  MSpMat;

//' Update for the penalized least squares problem
//' 
//' @param rv the matrix of row indices for the regular sparse column representation of ZtXt
//' @param xv the non-zero values in ZtXt
//' @param theta vector of relative standard deviations
//' @param m the sparse matrix to be updated.  Must have the correct pattern.
// [[Rcpp::export]]
void RSCupdate(const MiMat rv, const MMat xv, const MVec theta, MSpMat m) {
    int k(theta.size()), kpp(rv.rows()), n(rv.cols()), nz(m.nonZeros()), qpp(m.cols());
    int kppm1(kpp - 1);
    if (m.rows() != qpp) Rcpp::stop("u must be a square matrix");
    if (xv.rows() != kpp || xv.cols() != n) Rcpp::stop("sizes of rv and xv must agree");
    double *nzvals(m.valuePtr());
    const int *colptr(m.outerIndexPtr()), *rowvals(m.innerIndexPtr());
    for (int i = 0; i < nz; ++i) nzvals[i] = 0.; // initialize m to zeros
    
    for (int j = 0; j < n; ++j)	{ // iterate over columns of ZtXt
	for (int i = kppm1; i >= 0; --i) { // reverse order of rows of the RSC rep
	    int ii(rv(i, j));		// row in ZtXt and column in m
	    int cpi(colptr[ii]),	// first rowvals and nzvals index for this column
		ll(colptr[ii + 1] - 1); // should be location of diagonal
	    double rm(xv(i, j));	// non-zero value in ZtXt
	    if (rowvals[ll] != ii) Rcpp::stop("u is not upper triangular");
	    if (i < k) {	 // if in the Z part
		rm *= theta[i];	 // scale by theta[i]
		nzvals[ll] = 1.; // initialize diagonal element to 1.
	    }
	    nzvals[ll] += rm * rm;
	    --ll;		// move up the column of m to next non-zero
	    for (int l = i - 1; l >= 0 && ll >= cpi; --l, --ll) { // off-diagonals in m
		int ii1(rv(l,j)); // row in ZtXt
		if (rowvals[ll] > ii1) continue;
		if (rowvals[ll] != ii1) Rcpp::stop("Pattern mismatch");
		nzvals[ll] += rm * nzvals[ll] * (l < k ? theta[l] : 1.);
	    }
	}
    }
}
