#include "Gauss_Hermite.h"

#include <limits>
using namespace std;
using namespace Rcpp;

namespace lme4Eigen {

    /**
     * Generate zeros and weights of Hermite polynomial of order N, for
     * the AGQ method.
     *
     * Derived from Fortran code in package 'glmmML'
     *
     * @param x zeros of the polynomial, abscissas for AGQ
     * @param w weights used in AGQ
     */

    static void internal_ghq(NumericVector& x, NumericVector& w)
    {
	const double GHQ_EPS = 1e-15;
	const int GHQ_MAXIT = 40;
	int N = x.size();
	double Z = 0, HF = 0, HD = 0;
	double Z0, F0, F1, P, FD, Q, WP, GD, R, R1, R2;
	double HN = 1/double(N);
	std::vector<double> XX(N + 1), WW(N + 1);
	double *X = &XX[0], *W = &WW[0];
	
	for(int NR = 1; NR <= N / 2; NR++){
	    if(NR == 1)
		Z = -1.1611 + 1.46 * sqrt(double(N));
	    else
		Z -= HN * (N/2 + 1 - NR);
	    for (int IT = 0; IT <= GHQ_MAXIT; IT++) {
		Z0 = Z;
		F0 = 1.0;
		F1 = 2.0 * Z;
		for (int K = 2; K <= N; ++K){
		    HF = 2.0 * Z * F1 - 2.0 * double(K - 1.0) * F0;
		    HD = 2.0 * K * F1;
		    F0 = F1;
		    F1 = HF;
		}
		P = 1.0;
		for (int I = 1; I <= NR-1; ++I) P *= (Z - X[I]);
		FD = HF / P;
		Q = 0.0;
		for (int I = 1; I <= NR - 1; ++I){
		    WP = 1.0;
		    for(int J = 1; J <= NR - 1; ++J) if(J != I) WP *= ( Z - X[J] );
		    Q += WP;
		}
		GD = (HD-Q*FD)/P;
		Z -= (FD/GD);
		if (abs((Z - Z0) / Z) < GHQ_EPS) break;
	    }
	    
	    X[NR] = Z;
	    X[N+1-NR] = -Z;
	    R=1.0;
	    for (int K = 1; K <= N; ++K) R *= (2.0 * double(K));
	    W[N+1-NR] = W[NR] = 3.544907701811 * R / (HD*HD);
	}
	
	if( N % 2 ){
	    R1=1.0;
	    R2=1.0;
	    for (int J = 1; J <= N; ++J){
		R1=2.0*R1*J;
		if(J>=(N+1)/2) R2 *= J;
	    }
	    W[N/2+1]=0.88622692545276*R1/(R2*R2);
	    X[N/2+1]=0.0;
	}
	
	copy(XX.begin() + 1, XX.end(), x.begin());
	copy(WW.begin() + 1, WW.end(), w.begin());
    }

    GHQ::GHQ(int n)
	: d_xvals(n),
	  d_wts(n) {
	internal_ghq(d_xvals, d_wts);
    }

/** 
 * parchk checks parameters alpha and beta for classical weight functions. 
 * 
 * @param kind the rule
 * @param m the order of the highest moment to be calculated
 * @param alpha 
 * @param beta 
 */
    void parchk (int kind, int m, double alpha, double beta)     {
	double tmp;
	
	if (kind <= 0) throw std::runtime_error("parchk: kind <= 0");
	//  Check ALPHA for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
	if (3 <= kind && alpha <= -1.0)
	    throw std::runtime_error("parchk:  3 <= KIND and ALPHA <= -1.");
				//  Check BETA for Jacobi.
	if (kind == 4 && beta <= -1.0) throw std::runtime_error("parchk: KIND == 4 and BETA <= -1.0");
				//  Check ALPHA and BETA for rational.
	if (kind == 8) {
	    tmp = alpha + beta + m + 1.0;
	    if (0.0 <= tmp || tmp <= beta )
		throw std::runtime_error("parchk: kind == 8 but condition on alpha and beta fails");
	}
    }

/** 
 * Compute the Jacobi matrix for a quadrature rule.
 * 
 * This routine computes the diagonal AJ and sub-diagonal BJ
 * elements of the order M tridiagonal symmetric Jacobi matrix
 * associated with the polynomials orthogonal with respect to
 * the weight function specified by KIND.
 *
 * For weight functions 1-7, M elements are defined in BJ even
 * though only M-1 are needed.  For weight function 8, BJ(M) is
 * set to zero.
 *
 * The zero-th moment of the weight function is returned in ZEMU.
 * @param kind the rule
 * @param m the order of the Jacobi matrix
 * @param alpha 
 * @param beta 
 * @param aj diagonal of the Jacobi matrix
 * @param bj subdiagonal of the Jacobi matrix
 * 
 * @return the zero-th moment
 */
    double class_matrix (int kind, int m, double alpha, double beta, std::vector<double>& aj,
			 std::vector<double>& bj) {
	const double pi = 3.14159265358979323846264338327950;
	double a2b2;
	double ab;
	double aba;
	double abi;
	double abj;
	double abti;
	double apone;
	double zemu = 0.;

	parchk ( kind, 2 * m - 1, alpha, beta );

	if ( 500.0 *  std::numeric_limits<double>::epsilon() <
	     std::abs(std::pow(tgamma(0.5), 2) - pi))
	    throw std::runtime_error("Gamma function does not match machine parameters.");
	if ( kind == 1 ) {
	    ab = 0.0;
	    zemu = 2.0 / ( ab + 1.0 );
	    for (int i = 0; i < m; i++ ) aj[i] = 0.0;
	    for (int i = 1; i <= m; i++ ) {
		abi = i + ab * ( i % 2 );
		abj = 2 * i + ab;
		bj[i-1] = std::sqrt( abi * abi / ( abj * abj - 1.0 ) );
	    }
	} else if (kind == 2) {
	    zemu = pi;
	    for (int i = 0; i < m; i++ ) aj[i] = 0.0;
	    bj[0] = std::sqrt(0.5);
	    for (int i = 1; i < m; i++ ) bj[i] = 0.5;
	} else if ( kind == 3 ) {
	    ab = alpha * 2.0;
	    zemu = std::pow( 2.0, ab + 1.0 ) * std::pow( tgamma ( alpha + 1.0 ), 2 )
		/ tgamma ( ab + 2.0 );

	    for(int i = 0; i < m; i++ )
	    {
		aj[i] = 0.0;
	    }

	    bj[0] = std::sqrt( 1.0 / ( 2.0 * alpha + 3.0 ) );
	    for(int i = 2; i <= m; i++ )
	    {
		bj[i-1] = std::sqrt( i * ( i + ab ) / ( 4.0 * std::pow( i + alpha, 2 ) - 1.0 ) );
	    }
	}
	else if ( kind == 4 )
	{
	    ab = alpha + beta;
	    abi = 2.0 + ab;
	    zemu = std::pow( 2.0, ab + 1.0 ) * tgamma ( alpha + 1.0 ) 
		* tgamma ( beta + 1.0 ) / tgamma ( abi );
	    aj[0] = ( beta - alpha ) / abi;
	    bj[0] = std::sqrt( 4.0 * ( 1.0 + alpha ) * ( 1.0 + beta ) 
			   / ( ( abi + 1.0 ) * abi * abi ) );
	    a2b2 = beta * beta - alpha * alpha;

	    for(int i = 2; i <= m; i++ )
	    {
		abi = 2.0 * i + ab;
		aj[i-1] = a2b2 / ( ( abi - 2.0 ) * abi );
		abi = abi * abi;
		bj[i-1] = std::sqrt( 4.0 * i * ( i + alpha ) * ( i + beta ) * ( i + ab ) 
				 / ( ( abi - 1.0 ) * abi ) );
	    }
	}
	else if ( kind == 5 )
	{
	    zemu = tgamma ( alpha + 1.0 );

	    for(int i = 1; i <= m; i++ )
	    {
		aj[i-1] = 2.0 * i - 1.0 + alpha;
		bj[i-1] = std::sqrt( i * ( i + alpha ) );
	    }
	}
	else if ( kind == 6 )
	{
	    zemu = tgamma ( ( alpha + 1.0 ) / 2.0 );

	    for(int i = 0; i < m; i++ )
	    {
		aj[i] = 0.0;
	    }

	    for(int i = 1; i <= m; i++ )
	    {
		bj[i-1] = std::sqrt( ( i + alpha * ( i % 2 ) ) / 2.0 );
	    }
	}
	else if ( kind == 7 )
	{
	    ab = alpha;
	    zemu = 2.0 / ( ab + 1.0 );

	    for(int i = 0; i < m; i++ )
	    {
		aj[i] = 0.0;
	    }

	    for(int i = 1; i <= m; i++ )
	    {
		abi = i + ab * ( i % 2 );
		abj = 2 * i + ab;
		bj[i-1] = std::sqrt( abi * abi / ( abj * abj - 1.0 ) );
	    }
	}
	else if ( kind == 8 )
	{
	    ab = alpha + beta;
	    zemu = tgamma ( alpha + 1.0 ) * tgamma ( - ( ab + 1.0 ) ) 
		/ tgamma ( - beta );
	    apone = alpha + 1.0;
	    aba = ab * apone;
	    aj[0] = - apone / ( ab + 2.0 );
	    bj[0] = - aj[0] * ( beta + 1.0 ) / ( ab + 2.0 ) / ( ab + 3.0 );
	    for(int i = 2; i <= m; i++ )
	    {
		abti = ab + 2.0 * i;
		aj[i-1] = aba + 2.0 * ( ab + i ) * ( i - 1 );
		aj[i-1] = - aj[i-1] / abti / ( abti - 2.0 );
	    }

	    for(int i = 2; i <= m - 1; i++ )
	    {
		abti = ab + 2.0 * i;
		bj[i-1] = i * ( alpha + i ) / ( abti - 1.0 ) * ( beta + i ) 
		    / ( abti * abti ) * ( ab + i ) / ( abti + 1.0 );
	    }
	    bj[m-1] = 0.0;
	    for(int i = 0; i < m; i++ )
	    {
		bj[i] =  std::sqrt( bj[i] );
	    }
	}

	return zemu;
    }
//****************************************************************************80

    inline double r8_sign(const double& x) {return (x < 0) ? -1. : 1.;}

/** 
 * Diagonalize a symmetric tridiagonal matrix.
 * 
 * This routine is a slightly modified version of the EISPACK routine to 
 * perform the implicit QL algorithm on a symmetric tridiagonal matrix. 
 *
 * The authors thank the authors of EISPACK for permission to use this
 * routine. 
 *
 * Reference:
 *
 * Sylvan Elhay, Jaroslav Kautsky,
 * Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
 * Interpolatory Quadrature,
 * ACM Transactions on Mathematical Software,
 * Volume 13, Number 4, December 1987, pages 399-415.
 *
 * Roger Martin, James Wilkinson,
 * The Implicit QL Algorithm,
 * Numerische Mathematik,
 * Volume 12, Number 5, December 1968, pages 377-383.
 *
 * It has been modified to produce the product Q' * Z, where Z is an input 
 * vector and Q is the orthogonal matrix diagonalizing the input matrix.  
 * The changes consist (essentially) of applying the orthogonal transformations
 * directly to Z as they are generated.
 *
 * @param n the order of the matrix
 * @param d the diagonal entries of the matrix
 * @param e the subdiagonal entries of the matrix
 * @param z the value of Q' * Z where Q is the matrix that
 *          diagonalizes the input symmetric tridiagonal matrix 
 */
    void imtqlx (int n, std::vector<double>& d, std::vector<double>& e, std::vector<double>& z) {
	double b;
	double c;
	double f;
	double g;
	int i;
	int itn = 30;
	int m = 1;
	int mml;
	double p;
	double prec = std::numeric_limits<double>::epsilon();
	double r;
	double s;

	if (n == 1) return;

	e[n-1] = 0.0;

	for (int l = 1; l <= n; l++) {
	    int j = 0;
	    for ( ; ; ) {
		for (m = l; m <= n; m++ ) {
		    if ( m == n ) break; 

		    if (std::abs(e[m-1]) <= prec * (std::abs(d[m-1]) + std::abs(d[m]))) break;
		}
		p = d[l-1];
		if (m == l) break;
		if (itn <= j) throw std::runtime_error("imtqlx: Iteration limit exceeded");
		j++;
		g = ( d[l] - p ) / ( 2.0 * e[l-1] );
		r = std::sqrt(g * g + 1.0);
		g = d[m-1] - p + e[l-1] / (g + std::abs(r) * r8_sign(g));
		s = 1.0;
		c = 1.0;
		p = 0.0;
		mml = m - l;
		
		for(int ii = 1; ii <= mml; ii++ ) {
		    i = m - ii;
		    f = s * e[i-1];
		    b = c * e[i-1];
		    
		    if (std::abs(g) <= std::abs(f)) {
			c = g / f;
			r = std::sqrt(c * c + 1.0);
			e[i] = f * r;
			s = 1.0 / r;
			c = c * s;
		    } else {
			s = f / g;
			r = std::sqrt(s * s + 1.0);
			e[i] = g * r;
			c = 1.0 / r;
			s = s * c;
		    }
		    g = d[i] - p;
		    r = ( d[i-1] - g ) * s + 2.0 * c * b;
		    p = s * r;
		    d[i] = g + p;
		    g = c * r - b;
		    f = z[i];
		    z[i] = s * z[i-1] + c * f;
		    z[i-1] = c * z[i-1] - s * f;
		}
		d[l-1] = d[l-1] - p;
		e[l-1] = g;
		e[m-1] = 0.0;
	    }
	}
				//  Sorting.
	for(int ii = 2; ii <= m; ii++ ) {
	    int i = ii - 1;
	    int k = i;
	    double p = d[i-1];
	    
	    for (int j = ii; j <= n; j++ ) {
		if ( d[j-1] < p ) {
		    k = j;
		    p = d[j-1];
		}
	    }
	    
	    if (k != i) {
		d[k-1] = d[i-1];
		d[i-1] = p;
		p = z[i-1];
		z[i-1] = z[k-1];
		z[k-1] = p;
	    }
	}
    }


/** 
 * Scale a quadrature formula to a nonstandard interval.
 * 
 * @param nt number of knots
 * @param t the original knots
 * @param mlt the multiplicity of the knots
 * @param wts the weights
 * @param nwts the number of weights
 * @param ndx 
 * @param swts the scaled weights
 * @param st the scaled knots
 * @param kind 
 * @param alpha 
 * @param beta 
 * @param a 
 * @param b 
 */
    void scqf (int nt, std::vector<double>& t, std::vector<int> mlt, std::vector<double>& wts,
	       int nwts, std::vector<int>& ndx, 
	       std::vector<double>& swts, std::vector<double>& st, int kind,
	       double alpha, double beta, double a, 
	       double b) {
	double al;
	double be;
//	int i;
	int k;
	int l;
	double p;
	double shft = (a + b)/2.;
	double slp  = (b - a)/2.;
	double tmp;

	parchk(kind, 1, alpha, beta);

	switch (kind) {
	case 1:
	    al = 0.0;
	    be = 0.0;
	    if (std::abs( b - a ) <= std::numeric_limits<double>::epsilon())
		throw std::runtime_error("scqf: |B - A| too small");
	    break;
	case 2:
	    al = -0.5;
	    be = -0.5;
	    if (std::abs( b - a ) <= std::numeric_limits<double>::epsilon())
		throw std::runtime_error("scqf: |B - A| too small");
	    break;
	case 3:
	    al = alpha;
	    be = alpha;
	    if (std::abs( b - a ) <= std::numeric_limits<double>::epsilon())
		throw std::runtime_error("scqf: |B - A| too small");
	    break;
	case 4:
	    al = alpha;
	    be = beta;
	    if (std::abs( b - a ) <= std::numeric_limits<double>::epsilon())
		throw std::runtime_error("scqf: |B - A| too small");
	    break;
	case 5:
	    if ( b <= 0.0 ) throw std::runtime_error("scqf: b <= 0");
	    shft = a;
	    slp = 1.0 / b;
	    al = alpha;
	    be = 0.0;
	    break;
	case 6:
	    if ( b <= 0.0 ) throw std::runtime_error("scqf: b <= 0");
	    shft = a;
	    slp = 1.0 / std::sqrt( b );
	    al = alpha;
	    be = 0.0;
	    break;
	case 7:
	    al = alpha;
	    be = 0.0;
	    if (std::abs( b - a ) <= std::numeric_limits<double>::epsilon())
		throw std::runtime_error("scqf: |B - A| too small");
	    break;
	case 8:
	    if ( a + b <= 0.0 ) throw std::runtime_error("scqf: A + B <= 0.");
	    shft = a;
	    slp = a + b;
	    al = alpha;
	    be = beta;
	    break;
	case 9:
	    al = 0.5;
	    be = 0.5;
	    if (std::abs( b - a ) <= std::numeric_limits<double>::epsilon())
		throw std::runtime_error("scqf: |B - A| too small");
	    break;
	default:
	    throw std::runtime_error("unknown value of kind");
// FIXME: Make kind an enum
	}

	p = std::pow(slp, al + be + 1.0);

	for ( k = 0; k < nt; k++ ) {
	    st[k] = shft + slp * t[k];
	    l = std::abs(ndx[k]);
	    if ( l != 0 ) {
		tmp = p;
		for(int i = l - 1; i <= l - 1 + mlt[k] - 1; i++ ) {
		    swts[i] = wts[i] * tmp;
		    tmp = tmp * slp;
		}
	    }
	}
    }

/** 
 * Compute knots and weights of a Gauss Quadrature formula.
 * 
 * @param nt number of knots
 * @param aj the diagonal of the Jacobi matrix
 * @param bj the subdiagonal of the Jacobi matrix
 * @param zemu the zero-th moment of the weight function
 * @param t the knots
 * @param wts the weights
 */
    void sgqf (int nt, std::vector<double>& aj, std::vector<double>& bj,
	       double zemu, std::vector<double>& t, std::vector<double>& wts) {
	if (zemu <= 0.0) throw std::runtime_error("sgqf: zemu <= 0.");
				//  Set up vectors for IMTQLX.
	t = aj;
//	for(int i = 0; i < nt; i++ ) t[i] = aj[i];
	std::fill(wts.begin(), wts.end(), 0.);
	wts[0] = std::sqrt(zemu);
//	for(int i = 1; i < nt; i++ ) wts[i] = 0.0;
				//  Diagonalize the Jacobi matrix.
	imtqlx(nt, t, bj, wts);
//	wts *= wts;
	for(int i = 0; i < nt; i++ ) wts[i] = wts[i] * wts[i];
    }


/** 
 * Compute a Gauss quadrature formula with default A, B and simple knots.
 * 
 * This routine computes all the knots and weights of a Gauss quadrature
 * formula with a classical weight function with default values for A and B,
 * and only simple knots.
 *
 * @param nt number of knots
 * @param kind the rule
 * @param alpha 
 * @param beta 
 * @param t knots
 * @param wts weights
 */
    void cdgqf (int nt, int kind, double alpha, double beta, std::vector<double>& t, 
		std::vector<double>& wts)
    {
	std::vector<double> aj(nt), bj(nt);
	double zemu;

	parchk(kind, 2 * nt, alpha, beta);
				//  Get the Jacobi matrix and zero-th moment.
	zemu = class_matrix(kind, nt, alpha, beta, aj, bj);
				//  Compute the knots and weights.
	sgqf(nt, aj, bj, zemu, t, wts);
    }

/** 
 * Compute knots and weights of a Gauss quadrature formula.
 * 
 * @param nt number of knots
 * @param kind the rule
 * @param alpha 
 * @param beta 
 * @param a lower interval endpoint or location parameter
 * @param b upper interval endpoint or rate
 * @param t knots
 * @param wts weights
 */
    void cgqf (int nt, int kind, double alpha, double beta, double a, double b, 
	       std::vector<double>& t, std::vector<double>& wts)
    {
	std::vector<int> mlt(nt), ndx(nt);
				// Compute the Gauss quadrature formula for default values of A and B.
	cdgqf(nt, kind, alpha, beta, t, wts );
				// Prepare to scale the quadrature formula to other weight
				// function with valid A and B.
	std::fill(mlt.begin(), mlt.end(), 1);
	for (int i = 0; i < nt; i++) ndx[i] = i + 1;
	scqf(nt, t, mlt, wts, nt, ndx, wts, t, kind, alpha, beta, a, b);
    }

}
