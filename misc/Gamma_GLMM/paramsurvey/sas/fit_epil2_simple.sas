/* epil2 (simple): y ~ trt + (1|subject), Gamma(link=log).
   R beta order: (Intercept), trtprogabide -- nbeta=2. */

%include "fitlib.sas";

%fit_glimmix(
  example=epil2_simple,
  method=rspl,
  classvars=trt,
  modelrhs=%str(trt),
  nbeta=2,
  randomstmt=%str(random intercept / subject=subject solution;)
);

%fit_glimmix(
  example=epil2_simple,
  method=laplace,
  classvars=trt,
  modelrhs=%str(trt),
  nbeta=2,
  randomstmt=%str(random intercept / subject=subject solution;)
);
