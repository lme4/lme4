// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// predModule.h: predictor module using Eigen
//
// Copyright (C)       2011 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#ifndef LME4_PREDMODULE_H
#define LME4_PREDMODULE_H

#include "eigen.h"

using Rcpp::IntegerVector;
using Rcpp::NumericVector;
using Rcpp::S4;
using std::invalid_argument;
using std::runtime_error;

namespace lme4Eigen {

    class ddiMatrix : public DiagType {
    public:
	typedef typename DiagType::DiagonalVectorType  VType;
	ddiMatrix(const S4& xp)
	    : DiagType(Rcpp::as<VectorXd>(xp.slot("x"))) {}
	
	VType               diag()       {return diagonal();}
	double*              xPt()       {return diag().begin();}
	int                  nnz() const {return cols();}
    };
    
    class MdgCMatrix : public MSpMatrixXd { // mapped sparse Matrix, uses R's storage
    protected:
	S4    d_xp;
    public:
	MdgCMatrix(const S4& xp)
	    : MSpMatrixXd(Rcpp::as<MSpMatrixXd>(xp)), d_xp(xp) {}

//	const double*    xPt() const {return _valuePtr();}
//	int              nnz() const {return nonZeros();}
    };

    class dgCMatrix : public SpMatrixXd { // sparse Matrix, copies contents of R object
    public:
	dgCMatrix(const S4& xp)
	    : SpMatrixXd(Rcpp::as<SpMatrixXd>(xp)) {}

	double*          xPt()       {return _valuePtr();}
        VectorXd        diag() const {
	    VectorXd     ans(outerSize());
	    for (int j = 0; j < outerSize(); ++j) {
		SpMatrixXd::InnerIterator it(*this, j);
		if (it.index() != j)  // because *this is ColMajor and Lower
		    throw runtime_error("first element of column in lower is not a diagonal");
		ans[j] = it.value();
	    }
	    return ans;
	}
    };

    class modelMatrix {
    protected:
	const IntegerVector                d_assign;
	const Rcpp::List            d_contrasts;
    public:
	modelMatrix(const S4& xp)
	    : d_assign(   xp.slot("assign")),
	      d_contrasts(xp.slot("contrasts")) {}

	const IntegerVector         &assign() const {return d_assign;}
	const Rcpp::List  &contrasts() const {return d_contrasts;}
    };

    class ddenseModelMatrix : public MMatrixXd, public modelMatrix {
    protected:
	S4             d_xp;
    public:
	typedef MatrixXd           MatrixType;
	ddenseModelMatrix(const S4& xp)
	    : MMatrixXd(NumericVector(xp.slot("x")).begin(),
			::Rf_asInteger(xp.slot("Dim")),
			IntegerVector(xp.slot("Dim"))[1]),
	      modelMatrix(xp), d_xp(xp) {}
    };

    class dsparseModelMatrix : public MdgCMatrix, public modelMatrix {
    protected:
	S4             d_xp;
    public:
	typedef SpMatrixXd         MatrixType;
	dsparseModelMatrix(const S4& xp)
	    : MdgCMatrix(xp), modelMatrix(xp), d_xp(xp) {}
    };

#if 0
    namespace traits {
	template<typename Xtype> struct merPred_XTraits;
	template<> struct merPred_XTraits<ddenseModelMatrix> {
	    enum {
		RowsAtCompileTime    = ddenseModelMatrix::RowsAtCompileTime,
		ColsAtCompileTime    = ddenseModelMatrix::ColsAtCompileTime,
		MaxColsAtCompileTime = ddenseModelMatrix::MaxColsAtCompileTime
	    };
	    typedef typename ddenseModelMatrix::Index                           Index;
	    typedef typename ddenseModelMatrix::Scalar                          Scalar;
//	    typedef Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime> MatrixType;
//	    typedef Eigen::SelfAdjointView<MatrixType, Eigen::Lower>            LowerSymType;
	    // typedef dsyMatrixU                                                  UpperSymType;
	    // typedef Eigen::TriangularView<MatrixType, Eigen::Lower>             LowerTriType;
	    // typedef Eigen::TriangularView<MatrixType, Eigen::Upper>             UpperTriType;
	    // typedef Eigen::LLT<MatrixType, Eigen::Upper>                        LLtUpperType;
	    // typedef Eigen::LLT<MatrixType, Eigen::Lower>                        LLtLowerType;
	};

	template<> struct merPred_XTraits<dsparseModelMatrix> {
	    typedef typename dsparseModelMatrix::Index                          Index;
	    typedef typename dsparseModelMatrix::Scalar                         Scalar;
	    // typedef Eigen::SparseMatrix<Scalar, Eigen::ColMajor, Index>         MatrixType;
	    // typedef Eigen::SparseSelfAdjointView<MatrixType, Eigen::Lower>      LowerSymType;
	    // typedef Eigen::SparseSelfAdjointView<MatrixType, Eigen::Upper>      UpperSymType;
	    // typedef Eigen::SparseTriangularView<MatrixType, Eigen::Lower>       LowerTriType;
	    // typedef Eigen::SparseTriangularView<MatrixType, Eigen::Upper>       UpperTriType;
	    // typedef Eigen::SimplicialLLt<MatrixType, Eigen::Lower>              LLtLowerType;
	    // typedef Eigen::SimplicialLLt<MatrixType, Eigen::Upper>              LLtUpperType;
	};
    }
#endif

    class merPredD {
    public:
	typedef ddenseModelMatrix                                                  XType;
	typedef typename XType::Scalar                                             Scalar;
	typedef typename XType::Index                                              Index;
	typedef typename Eigen::Matrix<Scalar, XType::ColsAtCompileTime, 1>        VectorType;
	typedef typename Eigen::CholmodDecomposition<Eigen::SparseMatrix<double> > ChmDecomp;
    protected:
	XType                 d_X;
	MdgCMatrix            d_Zt;
	NumericVector   d_theta;
	IntegerVector   d_Lind;
	Index                 d_n, d_p, d_q, d_nnz;
	dgCMatrix             d_Lambdat;
	Scalar                d_ldL2, d_ldRX2;
	MatrixXd              d_RZX, d_V, d_VtV;
	VectorType            d_Vtr, d_Utr, d_delb, d_delu, d_beta0, d_u0;
	SpMatrixXd            d_Ut;//, d_I;
	ChmDecomp             d_L;
	bool                  d_isDiagLam;
	LLTType               d_RX;
    public:
	merPredD(S4, S4, S4, IntegerVector, NumericVector)        throw (invalid_argument,
								         runtime_error);
	IntegerVector          Pvec() const ;
	MatrixXd                 RX() const {MatrixXd ans = d_RX.matrixU(); return ans;}
	MatrixXd                RXi() const {return d_RX.matrixU().solve(MatrixXd::Identity(d_p,d_p));}
	MatrixXd                VtV() const {return d_VtV;}
	MatrixXd               unsc() const {MatrixXd rxi(RXi()); return rxi * rxi.adjoint();}
	SpMatrixXd        LambdatUt() const; 
	VectorType           RXdiag() const {return d_RX.matrixLLT().diagonal();}
	VectorType                b(const Scalar& fac) const {return d_Lambdat.adjoint()*u(fac);}
	VectorType             beta(const Scalar& fac) const {return d_beta0 + fac * d_delb;}
	VectorType          linPred(const Scalar& fac) const {
	    return d_X * beta(fac) + d_Zt.adjoint() * b(fac);}
	VectorType                u(const Scalar& fac) const {return d_u0 + fac * d_delu;}
	Scalar                 ldL2() const {return d_ldL2;}
	Scalar                ldRX2() const {return d_ldRX2;}
	Scalar                 sqrL(const Scalar& fac) const {return u(fac).squaredNorm();}
	int                    info() const {return d_L.info();}
	const MatrixXd&         RZX() const {return d_RZX;}
	const NumericVector&  theta() const {return d_theta;}
	const SpMatrixXd&   Lambdat() const {return d_Lambdat;}
	const MSpMatrixXd&       Zt() const {return d_Zt;}
	const VectorType&      delb() const {return d_delb;}
	const VectorType&      delu() const {return d_delu;}
	const VectorType&     beta0() const {return d_beta0;}
	const VectorType&        u0() const {return d_u0;}
	void            installPars(const Scalar& fac) {d_u0 = u(fac); d_beta0 = beta(fac);}
	void               setBeta0(const VectorType&)            throw (invalid_argument);
	void               setTheta(const NumericVector&)         throw (invalid_argument,
							                 runtime_error);
	void                  setU0(const VectorType& newU0)      throw (invalid_argument);
	void                  solve();
	void                 solveU() {d_delu = d_L.solve(d_Utr);}
	void           updateDecomp();
	void                updateL()                             throw (runtime_error);
	void              updateRes(const VectorType&)            throw (invalid_argument);
	void             updateXwts(const VectorType&)            throw (invalid_argument);

    };

}

extern "C" {
    SEXP merPredDCreate(SEXP, SEXP, SEXP, SEXP, SEXP); // constructor (returns external pointer)
    
    SEXP merPredDsetTheta(SEXP, SEXP); // setters
    SEXP merPredDsetBeta0(SEXP, SEXP);
    SEXP merPredDsetU0(SEXP, SEXP);

    SEXP merPredDLambdat(SEXP);	// getters
    SEXP merPredDL(SEXP);
    SEXP merPredDPvec(SEXP);
    SEXP merPredDRX(SEXP);
    SEXP merPredDRXdiag(SEXP);
    SEXP merPredDRZX(SEXP);
    SEXP merPredDVtV(SEXP);
    SEXP merPredDZt(SEXP);
    SEXP merPredDbeta0(SEXP);
    SEXP merPredDdelb(SEXP);
    SEXP merPredDdelu(SEXP);
    SEXP merPredDldL2(SEXP);
    SEXP merPredDldRX2(SEXP);
    SEXP merPredDtheta(SEXP);
    SEXP merPredDu0(SEXP);
    SEXP merPredDunsc(SEXP);

    SEXP merPredDlinPred(SEXP, SEXP); // methods
    SEXP merPredDinstallPars(SEXP, SEXP);
    SEXP merPredDsqrL(SEXP,SEXP);
}
#endif // LME4_EIGEN_H
