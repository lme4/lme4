// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// mcmcsamp.h: Markov-chain Monte Carlo sample class using Eigen
//
// Copyright (C)       2012 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.
#ifndef LME4_MCMCSAMP_H
#define LME4_MCMCSAMP_H

#include "predModule.h"
#include "respModule.h"

namespace lme4 {
    class mcmcsamp {
    public:
	typedef Eigen::ArrayXd   Ar1;
	typedef Eigen::Map<Ar1> MAr1;
	typedef Eigen::VectorXd  Vec;
	typedef Eigen::Map<Vec> MVec;
	typedef Eigen::ArrayXXd  Ar2;
	typedef Eigen::Map<Ar2> MAr2;
	typedef Eigen::MatrixXd  Mat;
	typedef Eigen::Map<Mat> MMat;
    protected:
	// lme4::merPredD *d_pred;
	// lme4::lmResp   *d_resp;
	MVec                    d_dev;
	MMat                    d_fixef;
	MVec                    d_sigma;
        MMat                    d_ranef;
    public:
				// all the work is done in the constructor
	mcmcsamp(lme4::merPredD *pred, lme4::lmResp *resp,
		 SEXP dev, SEXP fixef, SEXP sigma, SEXP ranef);
    };
}
    
#endif /* LME4_GLMFAMILY_H */

