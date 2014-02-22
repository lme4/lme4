//
// merPhylo.cpp: phylogeny object
//
// Copyright (C) 2014 Steve Walker
//
// This file is part of flexLambda/lme4

#include "predModule.h"
#include "merPhylo.h"

namespace lme4 {
    using   Rcpp::as;
    typedef Eigen::Map<MatrixXd> MMat;
    typedef Eigen::Map<VectorXd> MVec;
    typedef Eigen::Map<VectorXi> MiVec;
  
  merPhylo::merPhylo(SEXP edgeAncestor, SEXP edgeDescendent, SEXP edgeLength, SEXP Nnode)
    : d_edgeAncestor(  as<MiVec>(edgeAncestor)),
      d_edgeDescendent(as<MiVec>(edgeDescendent)),
      d_edgeLength(    as<MVec>(edgeLength)),
      d_Nnode(         as<int>(Nnode)) {
  }

  void merPhylo::updateEdgeLength(MVec newEdgeLength){
    d_edgeLength = newEdgeLength;
  }
}
