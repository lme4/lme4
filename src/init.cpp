#include "predModule.h"
#include "respModule.h"
#include "lmer.h"
#include <R_ext/Rdynload.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {

    CALLDEF(lmerDeviance, 3),

    CALLDEF(lmerRespCreate, 1),	 // generate external pointer

    CALLDEF(lmerRespsetREML, 2), // setters

    CALLDEF(lmerRespREML, 1),	 // getters

    CALLDEF(lmerRespLaplace, 4), // methods
    CALLDEF(lmerRespupdateMu, 2),

    CALLDEF(merPredDCreate, 5),	  // generate external pointer

    CALLDEF(merPredDsetTheta, 2), //setters
    CALLDEF(merPredDsetBeta0, 2),
    CALLDEF(merPredDsetU0, 2),

    CALLDEF(merPredDLambdat, 1), //getters
    CALLDEF(merPredDPvec, 1),
    CALLDEF(merPredDRX, 1),
    CALLDEF(merPredDRXdiag, 1),
    CALLDEF(merPredDRZX, 1),
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

    CALLDEF(merPredDlinPred, 2), // methods
    CALLDEF(merPredDinstallPars, 2),
    CALLDEF(merPredDsqrL, 2),	 

    CALLDEF(modRespsetOffset, 2),
    CALLDEF(modRespsetWeights, 2),

    CALLDEF(modRespmu, 1),
    CALLDEF(modRespoffset, 1),
    CALLDEF(modRespweights, 1),
    CALLDEF(modRespwrss, 1),
    CALLDEF(modRespwtres, 1),
    CALLDEF(modRespy, 1),
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

