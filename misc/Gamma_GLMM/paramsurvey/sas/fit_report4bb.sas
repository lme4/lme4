/* whale_crate/report4bb: crate ~ (1|location) + (1|fyear), Gamma(link=log).
   Intercept-only fixed effects (nbeta=1); two INDEPENDENT scalar random
   intercepts (no correlation between them -- separate RANDOM statements,
   not a single TYPE=UN block). sd1=sd(location), sd2=sd(fyear), corr=NA
   (matches the R side's convention for two independent grouping
   factors -- see toolkit.R's extract_sdcorr()). */

%include "fitlib.sas";

%fit_glimmix(
  example=report4bb,
  method=rspl,
  classvars=,
  modelrhs=%str(),
  nbeta=1,
  randomstmt=%str(
    random intercept / subject=location solution;
    random intercept / subject=fyear solution;
  )
);

%fit_glimmix(
  example=report4bb,
  method=laplace,
  classvars=,
  modelrhs=%str(),
  nbeta=1,
  randomstmt=%str(
    random intercept / subject=location solution;
    random intercept / subject=fyear solution;
  )
);
