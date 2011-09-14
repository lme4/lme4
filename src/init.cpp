#include "predModule.h"
#include "respModule.h"
#include "lmer.h"
#include <R_ext/Rdynload.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {

    CALLDEF(lmerDeviance, 3),

    CALLDEF(glm_Create, 2),	// generate external pointer

    CALLDEF(glm_setN, 2),	// setters

    CALLDEF(glm_devResid, 1),	// getters
    CALLDEF(glm_eta, 1),
    CALLDEF(glm_family, 1),
    CALLDEF(glm_link, 1),
    CALLDEF(glm_muEta, 1),
    CALLDEF(glm_n, 1),
    CALLDEF(glm_resDev, 1),
    CALLDEF(glm_sqrtWrkWt, 1),
    CALLDEF(glm_variance, 1),
    CALLDEF(glm_wrkResids, 1),
    CALLDEF(glm_wrkResp, 1),

    CALLDEF(glm_Laplace, 4),	// methods
    CALLDEF(glm_updateMu, 2),
    CALLDEF(glm_updateWts, 1),

    CALLDEF(lm_setOffset, 2),	// setters
    CALLDEF(lm_setWeights, 2),

    CALLDEF(lm_mu, 1),		// getters
    CALLDEF(lm_offset, 1),
    CALLDEF(lm_sqrtXwt, 1),
    CALLDEF(lm_sqrtrwt, 1),
    CALLDEF(lm_weights, 1),
    CALLDEF(lm_wrss, 1),
    CALLDEF(lm_wtres, 1),
    CALLDEF(lm_y, 1),

    CALLDEF(lm_updateMu, 2),	// method

    CALLDEF(lmer_Create, 1),	// generate external pointer

    CALLDEF(lmer_setREML, 2),   // setter

    CALLDEF(lmer_REML, 1),	// getter

    CALLDEF(lmer_Laplace, 4),   // method

    CALLDEF(merPredDCreate, 5),	// generate external pointer

    CALLDEF(merPredDsetTheta, 2), // setters
    CALLDEF(merPredDsetBeta0, 2),
    CALLDEF(merPredDsetU0, 2),

    CALLDEF(merPredDCcNumer, 1), // getters
    CALLDEF(merPredDL, 1),
    CALLDEF(merPredDLambdat, 1), 
    CALLDEF(merPredDLamtUt, 1),
    CALLDEF(merPredDPvec, 1),
    CALLDEF(merPredDRX, 1),
    CALLDEF(merPredDRXdiag, 1),
    CALLDEF(merPredDRZX, 1),
    CALLDEF(merPredDUt, 1),
    CALLDEF(merPredDV, 1),
    CALLDEF(merPredDVtV, 1),
    CALLDEF(merPredDZt, 1),
    CALLDEF(merPredDbeta0, 1),
    CALLDEF(merPredDdelb, 1),
    CALLDEF(merPredDdelu, 1),
    CALLDEF(merPredDldL2, 1),
    CALLDEF(merPredDldRX2, 1),
    CALLDEF(merPredDtheta, 1),
    CALLDEF(merPredDu0, 1),
    CALLDEF(merPredDunsc, 1),

    CALLDEF(merPredDb, 2),	// methods
    CALLDEF(merPredDbeta, 2),
    CALLDEF(merPredDlinPred, 2),
    CALLDEF(merPredDinstallPars, 2),
    CALLDEF(merPredDsolve, 1),
    CALLDEF(merPredDsolveU, 1),
    CALLDEF(merPredDsqrL, 2),
    CALLDEF(merPredDu, 2),
    CALLDEF(merPredDupdateDecomp, 1),
    CALLDEF(merPredDupdateRes, 2),
    CALLDEF(merPredDupdateXwts, 2),
    {NULL, NULL, 0}
};

/** Initializer for lme4Eigen, called upon loading the package.
 *
 *  Register routines that can be called directly from R.
 *  Initialize CHOLMOD and require the LL' form of the factorization.
 *  Install the symbols to be used by functions in the package.
 */
extern "C"
void R_init_lme4Eigen(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, (Rboolean)FALSE);
}

