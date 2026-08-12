/* schizophrenia: imps79 ~ TxDrug*Week + (1|id), Gamma(link=log).
   TxDrug is a plain 0/1 numeric in lme4::schizophrenia, NOT a factor
   (confirmed via str(schizophrenia)$TxDrug -> integer) -- so it's a
   continuous covariate here too, no CLASS statement needed at all
   (unlike epil2's trt). R beta order: (Intercept), TxDrug, Week,
   TxDrug:Week -- nbeta=4, term order matches R's Base*trt-style
   expansion since main effects already precede the interaction in the
   formula as written. */

%include "fitlib.sas";

%fit_glimmix(
  example=schizophrenia,
  method=rspl,
  classvars=,
  modelrhs=%str(TxDrug Week TxDrug*Week),
  nbeta=4,
  randomstmt=%str(random intercept / subject=id solution;)
);

%fit_glimmix(
  example=schizophrenia,
  method=laplace,
  classvars=,
  modelrhs=%str(TxDrug Week TxDrug*Week),
  nbeta=4,
  randomstmt=%str(random intercept / subject=id solution;)
);
