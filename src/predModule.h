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

namespace lme4Eigen {

    class ddiMatrix : public DiagType {
    public:
	typedef typename DiagType::DiagonalVectorType  VType;
	ddiMatrix(const Rcpp::S4& xp)
	    : DiagType(Rcpp::as<VectorXd>(xp.slot("x"))) {}
	
	VType               diag()       {return diagonal();}
	double*              xPt()       {return diag().begin();}
	int                  nnz() const {return cols();}
    };
    
    class MdgCMatrix : public MSpMatrixXd { // mapped sparse Matrix, uses R's storage
    protected:
	Rcpp::S4    d_xp;
    public:
	MdgCMatrix(const Rcpp::S4& xp)
	    : MSpMatrixXd(Rcpp::as<MSpMatrixXd>(xp)), d_xp(xp) {}

//	const double*    xPt() const {return _valuePtr();}
//	int              nnz() const {return nonZeros();}
    };

    class dgCMatrix : public SpMatrixXd { // sparse Matrix, copies contents of R object
    public:
	dgCMatrix(const Rcpp::S4& xp)
	    : SpMatrixXd(Rcpp::as<SpMatrixXd>(xp)) {}

	double*          xPt()       {return _valuePtr();}
        VectorXd        diag() const {
	    VectorXd     ans(outerSize());
	    for (int j = 0; j < outerSize(); ++j) {
		SpMatrixXd::InnerIterator it(*this, j);
		if (it.index() != j)  // because *this is ColMajor and Lower
		    throw std::runtime_error("first element of column in lower is not a diagonal");
		ans[j] = it.value();
	    }
	    return ans;
	}
    };

    class modelMatrix {
    protected:
	const Rcpp::IntegerVector                d_assign;
	const Rcpp::List            d_contrasts;
    public:
	modelMatrix(const Rcpp::S4& xp)
	    : d_assign(   xp.slot("assign")),
	      d_contrasts(xp.slot("contrasts")) {}

	const Rcpp::IntegerVector         &assign() const {return d_assign;}
	const Rcpp::List  &contrasts() const {return d_contrasts;}
    };

    class ddenseModelMatrix : public MMatrixXd, public modelMatrix {
    protected:
	Rcpp::S4             d_xp;
    public:
	typedef MatrixXd           MatrixType;
	ddenseModelMatrix(const Rcpp::S4& xp)
	    : MMatrixXd(Rcpp::NumericVector(xp.slot("x")).begin(),
			::Rf_asInteger(xp.slot("Dim")),
			Rcpp::IntegerVector(xp.slot("Dim"))[1]),
	      modelMatrix(xp), d_xp(xp) {}
    };

    class dsparseModelMatrix : public MdgCMatrix, public modelMatrix {
    protected:
	Rcpp::S4             d_xp;
    public:
	typedef SpMatrixXd         MatrixType;
	dsparseModelMatrix(const Rcpp::S4& xp)
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
	typedef ddenseModelMatrix                                               XType;
	typedef typename XType::Scalar                                          Scalar;
	typedef typename XType::Index                                           Index;
	typedef typename Eigen::Matrix<Scalar, XType::ColsAtCompileTime, 1>     VectorType;
	typedef typename Eigen::Array<Scalar, XType::ColsAtCompileTime, 1>      Array1DType;
    protected:
	XType                 d_X;
	MdgCMatrix            d_Z;
	Rcpp::NumericVector   d_theta;
	Rcpp::IntegerVector   d_Lind;
	Index                 d_n, d_p, d_q;
	dgCMatrix             d_Lambda;
	Scalar                d_ldL2, d_ldRX2;
	MatrixXd              d_RZX, d_V;
	VectorType            d_Vtr, d_Utr, d_delb, d_delu, d_beta0, d_u0;
	SpMatrixXd            d_U, d_I;
	SpLDLt                d_L;
	LLTType               d_RX;
    public:
	merPredD(Rcpp::S4, Rcpp::S4, Rcpp::S4,
		   Rcpp::IntegerVector, Rcpp::NumericVector)         throw (std::invalid_argument,
									    std::runtime_error);

	MatrixXd                      RX() const {MatrixXd ans = d_RX.matrixU(); return ans;}
	MatrixXd                     RXi() const {return d_RX.matrixU().solve(MatrixXd::Identity(d_p,d_p));}
	MatrixXd                    unsc() const {MatrixXd rxi(RXi()); return rxi * rxi.adjoint();}
	VectorType               Dvector() const {return d_L.vectorD();}
	VectorType               linPred(const double& fac) const {
	    return d_X * (d_beta0 + fac * d_delb) + d_Z * (d_Lambda * (d_u0 + fac * d_delu));
	}
	VectorType                     b(const Scalar& fac) const {
	    return d_Lambda * (d_u0 + fac * d_delu);}
	double                      ldL2() const {return d_ldL2;}
	double                     ldRX2() const {return d_ldRX2;}
	double                      sqrL(const double& fac) const {
	    return (d_u0 + fac * d_delu).squaredNorm();
	}
	int                         info() const {return d_L.info();}
	const MatrixXd&              RZX() const {return d_RZX;}
	const Rcpp::NumericVector& theta() const {return d_theta;}
	const PermutationType&         P() const {return d_L.permutationP();}
	const SpMatrixXd&         Lambda() const {return d_Lambda;}
	const MSpMatrixXd&             Z() const {return d_Z;}
	const VectorXi&             Pvec() const {return d_L.permutationP().indices();}
	const VectorType&           delb() const {return d_delb;}
	const VectorType&           delu() const {return d_delu;}
	const VectorType&          beta0() const {return d_beta0;}
	const VectorType&             u0() const {return d_u0;}
	void                 installPars(const double& fac) {
	    d_u0    = d_u0 + fac * d_delu;
	    d_beta0 = d_beta0 + fac * d_delb;
	}
	void                    setBeta0(const VectorType& nBeta)    throw (std::invalid_argument) {
	    if (d_beta0.size() != nBeta.size())
		throw std::invalid_argument("setBeta0: dimension mismatch");
	    d_beta0 = nBeta;	// does this cause a copy?
	}
	void                    setTheta(const Rcpp::NumericVector&) throw (std::invalid_argument,
									    std::runtime_error);
	void                       setU0(const VectorType& newU0)    throw (std::invalid_argument) {
	    if (d_u0.size() != newU0.size())
		throw std::invalid_argument("setU0: dimension mismatch");
	    d_u0 = newU0;
	}
	void                       solve();
	void                      solveU() {d_delu = d_L.solve(d_Utr);}
	void                     updateL(const SpMatrixXd&);
	void                   updateRes(const VectorType&,
					 const VectorType&)          throw (std::invalid_argument);
    };

}

extern "C" {
    SEXP merPredDCreate(SEXP, SEXP, SEXP, SEXP, SEXP);
    
    SEXP merPredDsetTheta(SEXP, SEXP);
    SEXP merPredDsetBeta0(SEXP, SEXP);
    SEXP merPredDsetU0(SEXP, SEXP);

    SEXP merPredDLambda(SEXP);
    SEXP merPredDPvec(SEXP);
    SEXP merPredDRX(SEXP);
    SEXP merPredDRZX(SEXP);
    SEXP merPredDZ(SEXP);
    SEXP merPredDbeta0(SEXP);
    SEXP merPredDdelb(SEXP);
    SEXP merPredDdelu(SEXP);
    SEXP merPredDldL2(SEXP);
    SEXP merPredDldRX2(SEXP);
    SEXP merPredDtheta(SEXP);
    SEXP merPredDu0(SEXP);
    SEXP merPredDunsc(SEXP);

    SEXP merPredDlinPred(SEXP, SEXP);
    SEXP merPredDinstallPars(SEXP, SEXP);
    SEXP merPredDsqrL(SEXP,SEXP);
}
#endif // LME4_EIGEN_H
