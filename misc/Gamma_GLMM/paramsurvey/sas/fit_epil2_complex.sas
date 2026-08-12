/* epil2 (complex): y ~ Base*trt + Age + Visit + (Visit|subject), Gamma(link=log).
   R beta order: (Intercept), Base, trtprogabide, Age, Visit,
   Base:trtprogabide -- nbeta=6. R's terms() moves the interaction to the
   end even though "Base*trt" is written first in the formula; SAS won't
   do that automatically, so the MODEL statement below is written in
   that exact final order by hand (Base trt Age Visit Base*trt), not in
   the more natural-looking "Base*trt Age Visit" order. */

%include "fitlib.sas";

%fit_glimmix(
  example=epil2_complex,
  method=rspl,
  classvars=trt,
  modelrhs=%str(Base trt Age Visit Base*trt),
  nbeta=6,
  randomstmt=%str(random intercept Visit / subject=subject type=un solution;)
);

%fit_glimmix(
  example=epil2_complex,
  method=laplace,
  classvars=trt,
  modelrhs=%str(Base trt Age Visit Base*trt),
  nbeta=6,
  randomstmt=%str(random intercept Visit / subject=subject type=un solution;)
);
