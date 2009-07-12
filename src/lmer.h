#ifndef LME4_LMER_H
#define LME4_LMER_H

#include <R.h>
/* Rdefines.h includes Rinternals.h (for SEXP, REAL, etc.) and defines
 * GET_SLOT, MAKE_CLASS, NEW_OBJECT, SET_SLOT, etc. */
#include <Rdefines.h>

/* SEXP lme4_rWishart(SEXP ns, SEXP dfp, SEXP scal); */
SEXP lme4_ghq(SEXP np);

SEXP mer_MCMCsamp(SEXP x, SEXP fm);
SEXP mer_ST_chol(SEXP x);
SEXP mer_ST_getPars(SEXP x);
SEXP mer_ST_initialize(SEXP ST, SEXP Gp, SEXP Zt);
SEXP mer_ST_setPars(SEXP x, SEXP pars);
SEXP mer_create_L(SEXP CmP);
SEXP mer_optimize(SEXP x);
SEXP mer_postVar(SEXP x, SEXP which);
SEXP mer_update_L(SEXP x);
SEXP mer_update_RX(SEXP x);
SEXP mer_update_dev(SEXP x);
SEXP mer_update_projection(SEXP x);
SEXP mer_update_ranef(SEXP x);
SEXP mer_update_mu(SEXP x);
SEXP mer_update_u(SEXP x);
SEXP mer_validate(SEXP x);
SEXP merMCMC_validate(SEXP x);
SEXP merMCMC_VarCorr(SEXP x, SEXP typ);

SEXP spR_optimize(SEXP x, SEXP verbP);
SEXP spR_update_mu(SEXP x);

#endif /* LME4_LMER_H */
