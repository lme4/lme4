// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// eigen.h: Rcpp/Eigen declarations for lme4
//
// Copyright (C)       2011 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#ifndef LME4_EIGEN_H
#define LME4_EIGEN_H

#include <RcppEigen.h>
				// Basic matrix, vector and array types for double
typedef Eigen::MatrixXd                                    MatrixXd;
typedef Eigen::VectorXd                                    VectorXd;
typedef Eigen::ArrayXd                                     ArrayXd;
typedef Eigen::ArrayXXd                                    ArrayXXd;
				// integer
typedef Eigen::MatrixXi                                    MatrixXi;
typedef Eigen::VectorXi                                    VectorXi;
typedef Eigen::ArrayXi                                     ArrayXi;
typedef Eigen::ArrayXXi                                    ArrayXXi;
				// complex
typedef Eigen::MatrixXcd                                   MatrixXcd;
typedef Eigen::VectorXcd                                   VectorXcd;
typedef Eigen::ArrayXcd                                    ArrayXcd;
typedef Eigen::ArrayXXcd                                   ArrayXXcd;
				// these should be defined separately for each base matrix type
typedef MatrixXd::Index                                    Index;
typedef MatrixXd::Scalar                                   Scalar;
typedef MatrixXd::RealScalar                               RealScalar;
				// Mapped matrix and vector types (share storage)
typedef Eigen::Map<MatrixXd>                               MMatrixXd;
typedef Eigen::Map<VectorXd>                               MVectorXd;
typedef Eigen::Map<ArrayXd>                                MArrayXd;
typedef Eigen::Map<ArrayXXd>                               MArrayXXd;
typedef Eigen::Map<MatrixXi>                               MMatrixXi;
typedef Eigen::Map<VectorXi>                               MVectorXi;
typedef Eigen::Map<ArrayXi>                                MArrayXi;
typedef Eigen::Map<ArrayXXi>                               MArrayXXi;
				// Views
typedef Eigen::TriangularView<MatrixXd, Eigen::Upper>      UpperTri;
typedef Eigen::TriangularView<MatrixXd, Eigen::Lower>      LowerTri;
typedef Eigen::SelfAdjointView<MatrixXd, Eigen::Upper>     UpperSym;
typedef Eigen::SelfAdjointView<MatrixXd, Eigen::Lower>     LowerSym;
				// Decomposition types
typedef Eigen::LLT<MatrixXd>                               LLTType;
typedef Eigen::LDLT<MatrixXd>                              LDLTType;
typedef Eigen::ColPivHouseholderQR<MatrixXd>               PivQRType;
typedef Eigen::HouseholderQR<MatrixXd>                     QRType;
typedef Eigen::JacobiSVD<MatrixXd>                         SVDType;
typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic>      DiagType;
typedef Eigen::SelfAdjointEigenSolver<MatrixXd>            EigvType;
				// Types derived from decompositions
typedef PivQRType::PermutationType                         PermutationType;
typedef LDLTType::TranspositionType                        TranspositionType;
typedef PermutationType::IndicesType                       IndicesType;
				// Sparse types
typedef Eigen::SparseVector<double>                              SpVectorXd;
typedef Eigen::SparseMatrix<double>                              SpMatrixXd;
typedef Eigen::MappedSparseMatrix<double>                        MSpMatrixXd;
typedef Eigen::SparseTriangularView<SpMatrixXd, Eigen::Upper>    SpUpperTri;
typedef Eigen::SparseTriangularView<SpMatrixXd, Eigen::Lower>    SpLowerTri;
typedef Eigen::SparseTriangularView<SpMatrixXd, Eigen::UnitUpper>SpUUpperTri;
typedef Eigen::SparseTriangularView<SpMatrixXd, Eigen::UnitLower>SpULowerTri;
typedef Eigen::SparseTriangularView<SpMatrixXd, Eigen::Lower>    SpLowerTri;
typedef Eigen::SparseSelfAdjointView<SpMatrixXd, Eigen::Lower>   SpLowerSym;
typedef Eigen::SparseSelfAdjointView<SpMatrixXd, Eigen::Upper>   SpUpperSym;
typedef Eigen::SimplicialLDLt<SpMatrixXd,Eigen::Lower>           SpLDLt;
typedef Eigen::SimplicialLLt<SpMatrixXd,Eigen::Lower>            SpLLt;

#endif
