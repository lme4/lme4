#include <Rcpp.h>

using namespace Rcpp;

//typedef std::vector<NumericMatrix> matV;

static void apply_lambda(const double *src,
			 const std::vector<NumericVector> &vv,
			 const std::vector<int> &sz,
			 std::vector<double> &dest) {
    for (int i = 0; i < dest.size(); ++i) dest[i] = src[i];
    int pos(0);
    for (int k = 0; k < vv.size(); ++k) {
	switch(sz[k]) {
	case 1: 
	    dest[pos++] *= vv[k][0];
	    break;
	case 2:
	    dest[pos] = dest[pos] * vv[k][0] + vv[k][2] * dest[pos + 1];
	    dest[++pos] *= vv[k][3];
	    pos++;
	    break;
	default:
	    Rcpp::stop("Code for vector-valued random effects not yet written");
	}
    }
}
//' Update for the penalized least squares problem
//' 
//' @param rv the matrix of row indices for the regular sparse column representation of ZtXt
//' @param xv the non-zero values in ZtXt
//' @param vc list of relative variance-covariace factor matrices
//' @param m the sparse matrix to be updated.  Must have the correct pattern.
//' @examples
//' rv <- rbind(as.integer(Dyestuff$Batch) - 1L, 6L)
//' xv <- matrix(1, nrow=nrow(rv), ncol=ncol(rv))
//' A <- sparseMatrix(i=as.integer(c(1,1,2,2,3,3,4,4,5,5,6,6,7)),
//'                   j=as.integer(c(1,7,2,7,3,7,4,7,5,7,6,7,7)),
//'                   x=rep.int(1.,13L), symmetric=TRUE,check=FALSE)
//' RSCupdate(rv, xv, list(matrix(0.831,1L,1L)), A)
//' A
//' @export
// [[Rcpp::export]]
void RSCupdate(const IntegerMatrix rv, const NumericMatrix xv,
	       const List vc, S4 m) {
    int nr(vc.size()), kpp(rv.nrow()), n(rv.ncol());
    std::vector<double> ww(kpp); // will hold a copy of a column of xv
    int kppm1(kpp - 1), k = 0;
    std::vector<NumericVector> vv(nr);
    std::vector<int> sz(nr);
    for (int i = 0; i < nr; ++i) {
	vv[i] = NumericVector(vc[i]);
	k += (sz[i] = (int)std::sqrt(vv[i].size()));
    }

    IntegerVector dims = m.slot("Dim");
    IntegerVector rowind = m.slot("i"), colptr = m.slot("p");
    NumericVector nzvals = m.slot("x");
    int nz(nzvals.size());
    if (dims[0] != dims[1]) Rcpp::stop("m must be a square matrix");
    if (xv.nrow() != kpp || xv.ncol() != n)
	Rcpp::stop("sizes of rv and xv must agree");
    for (int i = 0; i < nz; ++i) nzvals[i] = 0.; // initialize m to zeros

    for (int j = 0; j < n; ++j)	{ // iterate over columns of ZtXt
	apply_lambda(&xv(0,j), vv, sz, ww);
	Rcout << "j = " << j << ": ww[0:1] = " << ww[0] <<", "<<ww[1]<<std::endl;
	for (int i = kppm1; i >= 0; --i) { // reverse order of rows of the RSC rep
	    int ii(rv(i, j));	// row in ZtXt and column in m
	    int cpi(colptr[ii]), // first rowind and nzvals index for this column
		ll(colptr[ii + 1] - 1); // should be location of diagonal
	    if (rowind[ll] != ii) Rcpp::stop("u is not upper triangular");
	    if (i < k) nzvals[ll] = 1.; // add Identity in the Z part
	    nzvals[ll] += ww[i] * ww[i];
	    --ll;		// move up the column of m to next non-zero
	    Rcout << "i = " << i << ", ii = " << ii << ", ll = " << ll
		  << ", cpi = " << cpi << std::endl;
	    for (int l = i - 1; l >= 0 && ll >= cpi; --l) { // off-diagonals in m
		int ii1(rv(l,j)); // row in ZtXt
		while (rowind[ll] > ii1) --ll; // up to the desired row
		if (rowind[ll] != ii1) Rcpp::stop("Pattern mismatch");
		nzvals[ll] += ww[i] * ww[l];
	    }
	}
    }
}
