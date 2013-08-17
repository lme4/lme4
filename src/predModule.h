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

namespace lme4 {

    using Eigen::ArrayXd;
    using Eigen::LLT;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Eigen::VectorXi;

    class merPredD {
    public:
	typedef Eigen::Map<MatrixXd>                      MMap;
	typedef Eigen::Map<VectorXd>                      MVec;
	typedef Eigen::Map<VectorXi>                      MiVec;
	typedef MatrixXd::Scalar                          Scalar;
	typedef MatrixXd::Index                           Index;
	typedef Eigen::SparseMatrix<double>               SpMatrixd;
	typedef Eigen::CholmodDecomposition<SpMatrixd>    ChmDecomp;
	typedef Eigen::MappedSparseMatrix<double>         MSpMatrixd;
    protected:
	MMap          d_X, d_RZX, d_V, d_VtV;
	MSpMatrixd    d_Zt, d_Ut, d_LamtUt, d_Lambdat;
	MVec          d_theta, d_Vtr, d_Utr, d_Xwts, d_beta0, d_delb, d_delu, d_u0;
	MiVec         d_Lind;
	Index         d_N, d_p, d_q;
	Scalar        d_CcNumer, d_ldL2, d_ldRX2;
	ChmDecomp     d_L;
	LLT<MatrixXd> d_RX;
    public:
	merPredD(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, 
		 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

	VectorXi             Pvec() const;

	MatrixXd               RX() const;
	MatrixXd              RXi() const;
	MatrixXd             unsc() const;

	VectorXd           RXdiag() const;
	VectorXd                b(const Scalar& f) const;
	VectorXd             beta(const Scalar& f) const;
	VectorXd          linPred(const Scalar& f) const;
	VectorXd                u(const Scalar& f) const;

	Rcpp::List        condVar(const Rcpp::Environment&) const;

	Scalar            CcNumer() const {return d_CcNumer;}
	Scalar               ldL2() const {return d_ldL2;}
	Scalar              ldRX2() const {return d_ldRX2;}
	Scalar              solve();
	Scalar             solveU();
	Scalar               sqrL(const Scalar& f) const;

	const ChmDecomp&        L() const {return d_L;}

	const MMap&             V() const {return d_V;}
	const MMap&           VtV() const {return d_VtV;}
	const MMap&           RZX() const {return d_RZX;}

	const MSpMatrixd& Lambdat() const {return d_Lambdat;}
	const MSpMatrixd&  LamtUt() const {return d_LamtUt;}
	const MSpMatrixd&      Ut() const {return d_Ut;}
	const MSpMatrixd&      Zt() const {return d_Zt;}

	const MVec&           Utr() const {return d_Utr;}
	const MVec&           Vtr() const {return d_Vtr;}
	const MVec&          delb() const {return d_delb;}
	const MVec&          delu() const {return d_delu;}
	const MVec&         beta0() const {return d_beta0;}
	const MVec&         theta() const {return d_theta;}
	const MVec&            u0() const {return d_u0;}

	int                  info() const {return d_L.info();}

	void          installPars(const Scalar& f);
	void          MCMC_beta_u(const Scalar& sigma);
	void             setBeta0(const VectorXd&);
	void              setDelb(const VectorXd&);
	void              setDelu(const VectorXd&);
	void             setTheta(const VectorXd&);
	void                setU0(const VectorXd&);
	void         updateDecomp();
	void              updateL();
	void         updateLamtUt();
	void            updateRes(const VectorXd&);
	void           updateXwts(const  ArrayXd&);
    };
}

#endif // LME4_PREDMODULE_H
