#include "Gauss_Hermite.h"

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

    static void internal_ghq(NumericVector x, NumericVector w)
    {
	const double GHQ_EPS = 1e-15;
	const int GHQ_MAXIT = 40;
	int N = x.size(), NR, IT, I, K, J;
	double Z = 0, HF = 0, HD = 0;
	double Z0, F0, F1, P, FD, Q, WP, GD, R, R1, R2;
	double HN = 1/(double)N;
	NumericVector XX(N + 1), WW(N + 1);
	double *X = XX.begin(), *W = WW.begin();
	
	for(NR = 1; NR <= N / 2; NR++){
	    if(NR == 1)
		Z = -1.1611 + 1.46 * sqrt((double)N);
	    else
		Z -= HN * (N/2 + 1 - NR);
	    for (IT = 0; IT <= GHQ_MAXIT; IT++) {
		Z0 = Z;
		F0 = 1.0;
		F1 = 2.0 * Z;
		for(K = 2; K <= N; ++K){
		    HF = 2.0 * Z * F1 - 2.0 * (double)(K - 1.0) * F0;
		    HD = 2.0 * K * F1;
		    F0 = F1;
		    F1 = HF;
		}
		P = 1.0;
		for(I = 1; I <= NR-1; ++I){
		    P *= (Z - X[I]);
		}
		FD = HF / P;
		Q = 0.0;
		for(I = 1; I <= NR - 1; ++I){
		    WP = 1.0;
		    for(J = 1; J <= NR - 1; ++J){
			if(J != I) WP *= ( Z - X[J] );
		    }
		    Q += WP;
		}
		GD = (HD-Q*FD)/P;
		Z -= (FD/GD);
		if (abs((Z - Z0) / Z) < GHQ_EPS) break;
	    }
	    
	    X[NR] = Z;
	    X[N+1-NR] = -Z;
	    R=1.0;
	    for(K = 1; K <= N; ++K){
		R *= (2.0 * (double)K );
	    }
	    W[N+1-NR] = W[NR] = 3.544907701811 * R / (HD*HD);
	}
	
	if( N % 2 ){
	    R1=1.0;
	    R2=1.0;
	    for(J = 1; J <= N; ++J){
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
}














