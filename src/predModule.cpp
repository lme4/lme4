//
// predModule.cpp: implementation of predictor module using Eigen
//
// Copyright (C) 2011-2012 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#include "predModule.h"

namespace lme4 {
    using    Rcpp::as;

    using     std::invalid_argument;
    using     std::runtime_error;

    using   Eigen::ArrayXd;

    typedef Eigen::Map<MatrixXd>  MMat;
    typedef Eigen::Map<VectorXd>  MVec;
    typedef Eigen::Map<VectorXi>  MiVec;

    merPredD::merPredD(SEXP X, SEXP Lambdat, SEXP LamtUt, SEXP Lind,
		       SEXP RZX, SEXP Ut, SEXP Utr, SEXP V, SEXP VtV,
		       SEXP Vtr, SEXP Xwts, SEXP Zt, SEXP beta0,
		       SEXP delb, SEXP delu, SEXP theta, SEXP u0)
	: d_X(       as<MMat>(X)),
	  d_RZX(     as<MMat>(RZX)),
	  d_V(       as<MMat>(V)),
	  d_VtV(     as<MMat>(VtV)),
          d_Zt(      as<MSpMatrixd>(Zt)),
	  d_Ut(      as<MSpMatrixd>(Ut)),
	  d_LamtUt(  as<MSpMatrixd>(LamtUt)),
	  d_Lambdat( as<MSpMatrixd>(Lambdat)),
	  d_theta(   as<MVec>(theta)),
	  d_Vtr(     as<MVec>(Vtr)),
	  d_Utr(     as<MVec>(Utr)),
	  d_Xwts(    as<MVec>(Xwts)),
	  d_beta0(   as<MVec>(beta0)),
	  d_delb(    as<MVec>(delb)),
	  d_delu(    as<MVec>(delu)),
	  d_u0(      as<MVec>(u0)),
	  d_Lind(    as<MiVec>(Lind)),
	  d_N(       d_X.rows()),
	  d_p(       d_X.cols()),
	  d_q(       d_Zt.rows()),
	  d_RX(      d_p)
    {				// Check consistency of dimensions
	if (d_N != d_Zt.cols())
	    throw invalid_argument("Z dimension mismatch");
	if (d_Lind.size() != d_Lambdat.nonZeros())
	    throw invalid_argument("size of Lind does not match nonzeros in Lambda");
	// checking of the range of Lind is now done in R code for reference class
				// initialize beta0, u0, delb, delu and VtV
	d_VtV.setZero().selfadjointView<Eigen::Upper>().rankUpdate(d_V.adjoint());
	d_RX.compute(d_VtV);	// ensure d_RX is initialized even in the 0-column X case

	setTheta(d_theta);	    // starting values into Lambda
        d_L.cholmod().final_ll = 1; // force an LL' decomposition
	updateLamtUt();
	d_L.analyzePattern(d_LamtUt); // perform symbolic analysis
        if (d_L.info() != Eigen::Success)
	    throw runtime_error("CholeskyDecomposition.analyzePattern failed");
    }

    void merPredD::updateLamtUt() {
	// This complicated code bypasses problems caused by Eigen's
	// sparse/sparse matrix multiplication pruning zeros.  The
	// Cholesky decomposition croaks if the structure of d_LamtUt changes.
	MVec(d_LamtUt.valuePtr(), d_LamtUt.nonZeros()).setZero();
	for (Index j = 0; j < d_Ut.outerSize(); ++j) {
	    for(MSpMatrixd::InnerIterator rhsIt(d_Ut, j); rhsIt; ++rhsIt) {
		Scalar                        y(rhsIt.value());
		Index                         k(rhsIt.index());
		MSpMatrixd::InnerIterator prdIt(d_LamtUt, j);
		for (MSpMatrixd::InnerIterator lhsIt(d_Lambdat, k); lhsIt; ++lhsIt) {
		    Index                     i = lhsIt.index();
		    while (prdIt && prdIt.index() != i) ++prdIt;
		    if (!prdIt) throw runtime_error("logic error in updateLamtUt");
		    prdIt.valueRef()           += lhsIt.value() * y;
		}
	    }
	}
    }

    VectorXd merPredD::b(const double& f) const {return d_Lambdat.adjoint() * u(f);}

    VectorXd merPredD::beta(const double& f) const {return d_beta0 + f * d_delb;}

    VectorXd merPredD::linPred(const double& f) const {
	return d_X * beta(f) + d_Zt.adjoint() * b(f);
    }

    Rcpp::List merPredD::condVar(const Rcpp::Environment& rho) const {
	const Rcpp::List ll(as<Rcpp::List>(rho["flist"])), trmlst(as<Rcpp::List>(rho["terms"]));
	const int nf(ll.size());
	const MiVec nc(as<MiVec>(rho["ncols"])), nl(as<MiVec>(rho["nlevs"])),
	    nct(as<MiVec>(rho["nctot"])), off(as<MiVec>(rho["offsets"]));

	Rcpp::List ans(nf);
	ans.names() = clone(as<Rcpp::CharacterVector>(ll.names()));

	const SpMatrixd d_Lambda(d_Lambdat.adjoint());
	for (int i = 0; i < nf; i++) {
	    int         ncti(nct[i]), nli(nl[i]);
	    Rcpp::NumericVector ansi(ncti * ncti * nli);
	    ansi.attr("dim") = Rcpp::IntegerVector::create(ncti, ncti, nli);
	    ans[i] = ansi;
	    const MiVec trms(as<MiVec>(trmlst(i)));
	    if (trms.size() == 1) { // simple case
		int offset = off[trms[0] - 1];
		for (int j = 0; j < nli; ++j) {
		    MatrixXd Lv(d_Lambda.innerVectors(offset + j * ncti, ncti));
		    d_L.solveInPlace(Lv, CHOLMOD_A);
		    MatrixXd rr(MatrixXd(ncti, ncti).setZero().
				selfadjointView<Eigen::Lower>().rankUpdate(Lv.adjoint()));
		    std::copy(rr.data(), rr.data() + rr.size(), &ansi[j * ncti * ncti]);
		}
	    } else {
		throw std::runtime_error("multiple terms per factor not yet written");
	    }
	}
	return ans;
    }
	
    VectorXd merPredD::u(const double& f) const {return d_u0 + f * d_delu;}

    merPredD::Scalar merPredD::sqrL(const double& f) const {return u(f).squaredNorm();}

    void merPredD::updateL() {
	updateLamtUt();
	// More complicated code to handle the case of zeros in
	// potentially nonzero positions.  The factorize_p method is
	// for a SparseMatrix<double>, not a MappedSparseMatrix<double>.
	SpMatrixd  m(d_LamtUt.rows(), d_LamtUt.cols());
	m.resizeNonZeros(d_LamtUt.nonZeros());
	std::copy(d_LamtUt.valuePtr(),
		  d_LamtUt.valuePtr() + d_LamtUt.nonZeros(),
		  m.valuePtr()); 
	std::copy(d_LamtUt.innerIndexPtr(),
		  d_LamtUt.innerIndexPtr() + d_LamtUt.nonZeros(),
		  m.innerIndexPtr()); 
	std::copy(d_LamtUt.outerIndexPtr(),
		  d_LamtUt.outerIndexPtr() + d_LamtUt.cols() + 1,
		  m.outerIndexPtr()); 
	d_L.factorize_p(m, Eigen::ArrayXi(), 1.);
	d_ldL2 = ::M_chm_factor_ldetL2(d_L.factor());
    }

    void merPredD::setTheta(const VectorXd& theta) {
	if (theta.size() != d_theta.size())
	    throw invalid_argument("theta size mismatch");
				// update theta
	std::copy(theta.data(), theta.data() + theta.size(), d_theta.data());
				// update Lambdat
	int    *lipt = d_Lind.data();
	double *LamX = d_Lambdat.valuePtr(), *thpt = d_theta.data();
	for (int i = 0; i < d_Lind.size(); ++i) {
	    LamX[i] = thpt[lipt[i] - 1];
	}
    }

    merPredD::Scalar merPredD::solve() {
	d_delu          = d_Utr - d_u0;
	d_L.solveInPlace(d_delu, CHOLMOD_P);
	d_L.solveInPlace(d_delu, CHOLMOD_L);    // d_delu now contains cu
	d_CcNumer       = d_delu.squaredNorm(); // numerator of convergence criterion

	d_delb          = d_RX.matrixL().solve(d_Vtr - d_RZX.adjoint() * d_delu);
	d_CcNumer      += d_delb.squaredNorm(); // increment CcNumer
	d_RX.matrixU().solveInPlace(d_delb);

	d_delu         -= d_RZX * d_delb;
	d_L.solveInPlace(d_delu, CHOLMOD_Lt);
	d_L.solveInPlace(d_delu, CHOLMOD_Pt);
	return d_CcNumer;
    }

    merPredD::Scalar merPredD::solveU() {
	d_delb.setZero(); // in calculation of linPred delb should be zero after solveU
	d_delu          = d_Utr - d_u0;
	d_L.solveInPlace(d_delu, CHOLMOD_P);
	d_L.solveInPlace(d_delu, CHOLMOD_L);    // d_delu now contains cu
	d_CcNumer       = d_delu.squaredNorm(); // numerator of convergence criterion
	d_L.solveInPlace(d_delu, CHOLMOD_Lt);
	d_L.solveInPlace(d_delu, CHOLMOD_Pt);
	return d_CcNumer;
    }

    void merPredD::updateXwts(const ArrayXd& sqrtXwt) {
	if (d_Xwts.size() != sqrtXwt.size())
	    throw invalid_argument("updateXwts: dimension mismatch");
	std::copy(sqrtXwt.data(), sqrtXwt.data() + sqrtXwt.size(), d_Xwts.data());
	if (sqrtXwt.size() == d_V.rows()) { // W is diagonal
	    d_V              = d_Xwts.asDiagonal() * d_X;
	    for (int j = 0; j < d_N; ++j)
		for (MSpMatrixd::InnerIterator Utj(d_Ut, j), Ztj(d_Zt, j);
		     Utj && Ztj; ++Utj, ++Ztj)
		    Utj.valueRef() = Ztj.value() * d_Xwts.data()[j];
	} else {
	    SpMatrixd      W(d_V.rows(), sqrtXwt.size());
	    const double *pt = sqrtXwt.data();
	    W.reserve(sqrtXwt.size());
	    for (Index j = 0; j < W.cols(); ++j, ++pt) {
		W.startVec(j);
		W.insertBack(j % d_V.rows(), j) = *pt;
	    }
	    W.finalize();
	    d_V              = W * d_X;
	    SpMatrixd      Ut(d_Zt * W.adjoint());
	    if (Ut.cols() != d_Ut.cols())
		throw std::runtime_error("Size mismatch in updateXwts");

	    // More complex code to handle the pruning of zeros
	    MVec(d_Ut.valuePtr(), d_Ut.nonZeros()).setZero();
	    for (int j = 0; j < d_Ut.outerSize(); ++j) {
		MSpMatrixd::InnerIterator lhsIt(d_Ut, j);
		for (SpMatrixd::InnerIterator  rhsIt(Ut, j); rhsIt; ++rhsIt, ++lhsIt) {
		    Index                         k(rhsIt.index());
		    while (lhsIt && lhsIt.index() != k) ++lhsIt;
		    if (lhsIt.index() != k)
			throw std::runtime_error("Pattern mismatch in updateXwts");
		    lhsIt.valueRef() = rhsIt.value();
		}
	    }
	}
	d_VtV.setZero().selfadjointView<Eigen::Upper>().rankUpdate(d_V.adjoint());
	updateL();
    }

    void merPredD::updateDecomp() { // update L, RZX and RX
	updateL();
	d_RZX         = d_LamtUt * d_V;
	if (d_p > 0) {
	    d_L.solveInPlace(d_RZX, CHOLMOD_P);
	    d_L.solveInPlace(d_RZX, CHOLMOD_L);

	    MatrixXd      VtVdown(d_VtV);
	    d_RX.compute(VtVdown.selfadjointView<Eigen::Upper>().rankUpdate(d_RZX.adjoint(), -1));
	    if (d_RX.info() != Eigen::Success)
		::Rf_error("Downdated VtV is not positive definite");
	    d_ldRX2         = 2. * d_RX.matrixLLT().diagonal().array().abs().log().sum();
	}
    }

    void merPredD::updateRes(const VectorXd& wtres) {
	if (d_V.rows() != wtres.size())
	    throw invalid_argument("updateRes: dimension mismatch");
	d_Vtr           = d_V.adjoint() * wtres;
	d_Utr           = d_LamtUt * wtres;
    }

    void merPredD::installPars(const Scalar& f) {
	d_u0 = u(f);
	d_beta0 = beta(f);
	d_delb.setZero();
	d_delu.setZero();
    }

    void merPredD::setBeta0(const VectorXd& nBeta) {
	if (nBeta.size() != d_p) throw invalid_argument("setBeta0: dimension mismatch");
	std::copy(nBeta.data(), nBeta.data() + d_p, d_beta0.data());
    }

    void merPredD::setDelb(const VectorXd& newDelb) {
	if (newDelb.size() != d_p)
	    throw invalid_argument("setDelb: dimension mismatch");
	std::copy(newDelb.data(), newDelb.data() + d_p, d_delb.data());
    }

    void merPredD::setDelu(const VectorXd& newDelu) {
	if (newDelu.size() != d_q)
	    throw invalid_argument("setDelu: dimension mismatch");
	std::copy(newDelu.data(), newDelu.data() + d_q, d_delu.data());
    }

    void merPredD::setU0(const VectorXd& newU0) {
	if (newU0.size() != d_q) throw invalid_argument("setU0: dimension mismatch");
	std::copy(newU0.data(), newU0.data() + d_q, d_u0.data());
    }

    template <typename T>
    struct Norm_Rand : std::unary_function<T, T> {
	const T operator()(const T& x) const {return ::norm_rand();}
    };

    inline static VectorXd Random_Normal(int size, double sigma) {
	return ArrayXd(size).unaryExpr(Norm_Rand<double>()) * sigma;
    }

    void merPredD::MCMC_beta_u(const Scalar& sigma) {
	VectorXd del2(d_RX.matrixU().solve(Random_Normal(d_p, sigma)));
	d_delb += del2;
	VectorXd del1(Random_Normal(d_q, sigma) - d_RZX * del2);
	d_L.solveInPlace(del1, CHOLMOD_Lt);
	d_delu += del1;
    }

    VectorXi merPredD::Pvec() const {
	int*                 ppt((int*)d_L.factor()->Perm);
	VectorXi             ans(d_q);
	std::copy(ppt, ppt + d_q, ans.data());
	return ans;
    }

    MatrixXd merPredD::RX() const {
	return d_RX.matrixU();
    }

    MatrixXd merPredD::RXi() const {
	return d_RX.matrixU().solve(MatrixXd::Identity(d_p,d_p));
    }

    MatrixXd merPredD::unsc() const {
	return MatrixXd(MatrixXd(d_p, d_p).setZero().
			selfadjointView<Eigen::Lower>().
			rankUpdate(RXi()));
    }

    VectorXd merPredD::RXdiag() const {
	return d_RX.matrixLLT().diagonal();
    }
}
