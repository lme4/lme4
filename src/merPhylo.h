//
// merPhylo.h: phylogeny object
//
// Copyright (C)       2014 Steve Walker
//
// This file is part of flexLambda/lme4.

#ifndef LME4_MERPHYLO_H
#define LME4_MERPHYLO_H

#include <RcppEigen.h>

namespace lme4 {

    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Eigen::VectorXi;

    class merPhylo {
    public:
	typedef Eigen::Map<MatrixXd> MMap;
	typedef Eigen::Map<VectorXd> MVec;
	typedef Eigen::Map<VectorXi> MiVec;

    protected:
	MiVec d_edgeAncestor, d_edgeDescendent;
	MVec  d_edgeLength;
	int   d_Nnode;
    public:
	merPhylo(SEXP,SEXP,SEXP,SEXP);
	void  updateEdgeLength(MVec);
    };
}

#endif // LME4_MERPHYLO_H
