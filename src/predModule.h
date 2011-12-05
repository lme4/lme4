// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// predModule.h: predictor module using Eigen
//
// Copyright (C)       2011 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#ifndef LME4_PREDMODULE_H
#define LME4_PREDMODULE_H

#include <RcppEigen.h>

namespace lme4Eigen {
    using Eigen::ArrayXi;
    using Eigen::CholmodDecomposition;
    using Eigen::Dynamic;
    using Eigen::LLT;
    using Eigen::Lower;
    using Eigen::Map;
    using Eigen::Matrix;
    using Eigen::MatrixXd;
    using Eigen::SelfAdjointView;
    using Eigen::SparseMatrix;
    using Eigen::Upper;
    using Eigen::VectorXd;
    using Eigen::VectorXi;

    using Rcpp::IntegerVector;
    using Rcpp::List;
    using Rcpp::NumericVector;
    using Rcpp::S4;
    using Rcpp::as;

    using std::runtime_error;
    
    class dgCMatrix : public SparseMatrix<double> { // sparse Matrix, copies contents of R object
    public:
	dgCMatrix(const S4& xp)
	    : SparseMatrix<double>(as<SparseMatrix<double> >(xp)) {}

	double*          xPt()       {return _valuePtr();}
        VectorXd        diag() const {
	    VectorXd     ans(outerSize());
	    for (int j = 0; j < outerSize(); ++j) {
		SparseMatrix<double>::InnerIterator it(*this, j);
		if (it.index() != j)  // because *this is ColMajor and Lower
		    throw runtime_error("first element of column in lower is not a diagonal");
		ans[j] = it.value();
	    }
	    return ans;
	}
    };

    class modelMatrix {
    protected:
	const IntegerVector          d_assign;
	const List                   d_contrasts;
    public:
	modelMatrix(const S4& xp)
	    : d_assign(   xp.slot("assign")),
	      d_contrasts(xp.slot("contrasts")) {}

	const IntegerVector         &assign() const {return d_assign;}
	const List               &contrasts() const {return d_contrasts;}
    };

    class ddenseModelMatrix : public Map<MatrixXd>, public modelMatrix {
    protected:
	S4             d_xp;
    public:
	typedef MatrixXd           MatrixType;
	ddenseModelMatrix(const S4& xp)
	    : Map<MatrixXd>(NumericVector(xp.slot("x")).begin(),
			    ::Rf_asInteger(xp.slot("Dim")),
			    IntegerVector(xp.slot("Dim"))[1]),
	      modelMatrix(xp), d_xp(xp) {}
    };

#if 0
    class dsparseModelMatrix : public MappedSparseMatrix<double>, public modelMatrix {
    protected:
	S4             d_xp;
    public:
	typedef SparseMatrix<double> MatrixType;
	dsparseModelMatrix(const S4& xp)
	    : MappedSparseMatrix<double>(as<MappedSparseMatrix<double> >(xp)),
	      modelMatrix(xp), d_xp(xp) {}
    };

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
//	    typedef Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime> MatrixType;
//	    typedef SelfAdjointView<MatrixType, Lower>            LowerSymType;
	    // typedef dsyMatrixU                                                  UpperSymType;
	    // typedef TriangularView<MatrixType, Lower>             LowerTriType;
	    // typedef TriangularView<MatrixType, Upper>             UpperTriType;
	    // typedef Eigen::LLT<MatrixType, Upper>                        LLtUpperType;
	    // typedef Eigen::LLT<MatrixType, Lower>                        LLtLowerType;
	};

	template<> struct merPred_XTraits<dsparseModelMatrix> {
	    typedef typename dsparseModelMatrix::Index                          Index;
	    typedef typename dsparseModelMatrix::Scalar                         Scalar;
	    // typedef SparseMatrix<Scalar, Eigen::ColMajor, Index>         MatrixType;
	    // typedef Eigen::SparseSelfAdjointView<MatrixType, Lower>      LowerSymType;
	    // typedef Eigen::SparseSelfAdjointView<MatrixType, Upper>      UpperSymType;
	    // typedef Eigen::SparseTriangularView<MatrixType, Lower>       LowerTriType;
	    // typedef Eigen::SparseTriangularView<MatrixType, Upper>       UpperTriType;
	    // typedef Eigen::SimplicialLLt<MatrixType, Lower>              LLtLowerType;
	    // typedef Eigen::SimplicialLLt<MatrixType, Upper>              LLtUpperType;
	};
    }
#endif

    class merPredD {
    public:
	typedef ddenseModelMatrix                            XType;
	typedef XType::Scalar                                Scalar;
	typedef XType::Index                                 Index;
	typedef CholmodDecomposition<SparseMatrix<double> >  ChmDecomp;
	typedef SparseMatrix<double>                         SpMatrixd;
	typedef Eigen::MappedSparseMatrix<double>            MSpMatrixd;
    protected:
	XType          d_X;
	MSpMatrixd     d_Zt;
	Map<VectorXd>  d_theta;
	Map<VectorXi>  d_Lind;
	Index          d_n, d_nnz, d_p, d_q;
	MSpMatrixd     d_Lambdat;
	Scalar         d_CcNumer, d_ldL2, d_ldRX2;
	MatrixXd       d_RZX, d_V, d_VtV;
	VectorXd       d_Vtr, d_Utr, d_delb, d_delu;
	Map<VectorXd>  d_beta0, d_u0;
	SpMatrixd      d_Ut, d_LamtUt;
	ChmDecomp      d_L;
	LLT<MatrixXd>  d_RX;
	bool           d_LamtUtRestructure;
    public:
	merPredD(S4, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

	IntegerVector          Pvec() const;

	MatrixXd                 RX() const;
	MatrixXd                RXi() const;
	MatrixXd               unsc() const;

	VectorXd             RXdiag() const;
	VectorXd                  b(const Scalar& f) const;
	VectorXd               beta(const Scalar& f) const;
	VectorXd            linPred(const Scalar& f) const;
	VectorXd                  u(const Scalar& f) const;

	Scalar              CcNumer() const {return d_CcNumer;}
	Scalar                 ldL2() const {return d_ldL2;}
	Scalar                ldRX2() const {return d_ldRX2;}
	Scalar                solve();
	Scalar               solveU();
	Scalar                 sqrL(const Scalar& f) const;

	const ChmDecomp&          L() const {return d_L;}

	const MatrixXd&           V() const {return d_V;}
	const MatrixXd&         VtV() const {return d_VtV;}
	const MatrixXd&         RZX() const {return d_RZX;}

	const Map<VectorXd>&  theta() const {return d_theta;}

	const MSpMatrixd&   Lambdat() const {return d_Lambdat;}
	const SpMatrixd&     LamtUt() const {return d_LamtUt;}
	const SpMatrixd&         Ut() const {return d_Ut;}

	const MSpMatrixd&        Zt() const {return d_Zt;}

	const VectorXd&         Utr() const {return d_Utr;}
	const VectorXd&         Vtr() const {return d_Vtr;}
	const VectorXd&        delb() const {return d_delb;}
	const VectorXd&        delu() const {return d_delu;}
	const Map<VectorXd>&  beta0() const {return d_beta0;}
	const Map<VectorXd>&     u0() const {return d_u0;}

	int                    info() const {return d_L.info();}

	void            installPars(const Scalar& f);
	void               setBeta0(const VectorXd&);
	void                   setS(int);
	void               setTheta(const VectorXd&);
	void                  setU0(const VectorXd&);
	void           updateDecomp();
	void                updateL();
	void           updateLamtUt();
	void              updateRes(const VectorXd&);
	void             updateXwts(const VectorXd&);
    };
}

#endif // LME4_PREDMODULE_H
