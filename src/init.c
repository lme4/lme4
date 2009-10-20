#include "lmer.h"
#include <R_ext/Rdynload.h>
#include "Matrix.h"
#include "Syms.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {

    CALLDEF(lme4_ghq, 1),
    CALLDEF(mer_MCMCsamp, 2),
    CALLDEF(mer_ST_chol, 1),
    CALLDEF(mer_ST_getPars, 1),
    CALLDEF(mer_ST_initialize, 3),
    CALLDEF(mer_ST_setPars, 2),
    CALLDEF(mer_create_L, 1),
    CALLDEF(mer_optimize, 1),
    CALLDEF(mer_postVar, 2),
    CALLDEF(mer_update_L, 1),
    CALLDEF(mer_update_RX, 1),
    CALLDEF(mer_update_dev, 1),
    CALLDEF(mer_update_projection, 1),
    CALLDEF(mer_update_ranef, 1),
    CALLDEF(mer_update_mu, 1),
    CALLDEF(mer_update_u, 1),
    CALLDEF(mer_validate, 1),

    CALLDEF(merMCMC_validate, 1),
    CALLDEF(merMCMC_VarCorr, 2),

    CALLDEF(spR_optimize, 2),
    CALLDEF(spR_update_mu, 1),

  {NULL, NULL, 0}
};

/** cholmod_common struct local to the lme4 package */
cholmod_common c;

/** Need our own CHOLMOD error handler */
void attribute_hidden
lme4_R_cholmod_error(int status, const char *file, int line, const char *message)
{
    if(status < 0) {
#ifdef Failure_in_matching_Matrix
/* This fails unexpectedly with
 *  function 'cholmod_l_defaults' not provided by package 'Matrix'
 * from ../tests/lmer-1.R 's  (l.68)  lmer(y ~ habitat + (1|habitat*lagoon)
 */
	M_cholmod_defaults(&c);/* <--- restore defaults,
				* as we will not be able to .. */
	c.final_ll = 1;	    /* LL' form of simplicial factorization */
#endif

	error("Cholmod error '%s' at file:%s, line %d", message, file, line);
    }
    else
	warning("Cholmod warning '%s' at file:%s, line %d",
		message, file, line);
}

/** Initializer for lme4, called upon loading the package.
 *
 *  Register routines that can be called directly from R.
 *  Initialize CHOLMOD and require the LL' form of the factorization.
 *  Install the symbols to be used by functions in the package.
 */
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
void R_init_lme4(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);


    M_R_cholmod_start(&c);
    c.final_ll = 1;	    /* LL' form of simplicial factorization */

    /* need own error handler, that resets  final_ll (after *_defaults()) : */
    c.error_handler = lme4_R_cholmod_error;

    lme4_ASym = install("A");
    lme4_CmSym = install("Cm");
    lme4_CxSym = install("Cx");
    lme4_DimSym = install("Dim");
    lme4_GpSym = install("Gp");
    lme4_LSym = install("L");
    lme4_RXSym = install("RX");
    lme4_RZXSym = install("RZX");
    lme4_STSym = install("ST");
    lme4_VSym = install("V");
    lme4_XSym = install("X");
    lme4_XstSym = install("Xst");
    lme4_ZtSym = install("Zt");
    lme4_devianceSym = install("deviance");
    lme4_dimsSym = install("dims");
    lme4_envSym = install("env");
    lme4_etaSym = install("eta");
    lme4_fixefSym = install("fixef");
    lme4_flistSym = install("flist");
    lme4_ghwSym = install("ghw");
    lme4_ghxSym = install("ghx");
    lme4_gradientSym = install("gradient");
    lme4_iSym = install("i");
    lme4_ncSym = install("nc");
    lme4_nlmodelSym = install("nlmodel");
    lme4_muEtaSym = install("muEta");
    lme4_muSym = install("mu");
    lme4_offsetSym = install("offset");
    lme4_pSym = install("p");
    lme4_permSym = install("perm");
    lme4_pWtSym = install("pWt");
    lme4_ranefSym = install("ranef");
    lme4_residSym = install("resid");
    lme4_sigmaSym = install("sigma");
    lme4_sqrtrWtSym = install("sqrtrWt");
    lme4_sqrtXWtSym = install("sqrtXWt");
    lme4_uSym = install("u");
    lme4_varSym = install("var");
    lme4_xSym = install("x");
    lme4_ySym = install("y");
}

/** Finalizer for lme4 called upon unloading the package.
 *
 */
void R_unload_lme4(DllInfo *dll){
    M_cholmod_finish(&c);
}
