/* Shared macro for the paramsurvey SAS arm: fit PROC GLIMMIX to all B
   replicates of one example (via BY-group processing on a `rep` column,
   see ../12_export_csv_for_sas.R) under one estimation METHOD, and write
   a results CSV in the same column shape 10_ingest_julia_results.R /
   13_ingest_sas_results.R already expect: i, status, singular, msg,
   time_sec, sd1, sd2, corr, phi, negll, beta1..betaN (positional, same
   convention as the Julia arm -- SAS's own effect-name strings, e.g.
   "trt" vs R's "trtprogabide", don't string-match R's factor-contrast
   names any more than Julia's coefnames() do).

   NOT YET RUN AGAINST REAL SAS -- written from PROC GLIMMIX/ODS OUTPUT
   documentation, not verified interactively. Things flagged "VERIFY"
   below are the parts most likely to need adjustment against real ODS
   OUTPUT table contents (exact CovParm label text, ConvergenceStatus
   column names, Fit Statistics row labels can all vary a little by SAS
   release/options).

   Timing caveat: BY-group processing fits all B replicates inside one
   PROC GLIMMIX step, so there's no per-replicate wall-clock split the
   way R's system.time()-per-fit or Julia's per-file timing gives.
   time_sec here is the *average* over the whole batch (total elapsed /
   B, identical for every row), not a real per-replicate measurement.

   Dispersion caveat: METHOD=RSPL's Fit Statistics are on the linearized
   pseudo-data scale and are NOT a real -2*logLik comparable to the other
   methods -- negll is left missing for method=rspl rather than reporting
   a misleadingly-comparable-looking number. Only method=laplace's negll
   is filled in. */

%macro fit_glimmix(
    example=,           /* e.g. epil2_simple -- must match <example>.csv */
    method=,            /* rspl | laplace (PROC GLIMMIX METHOD= value) */
    classvars=,         /* space-separated CLASS variables, or empty */
    modelrhs=,          /* fixed-effect RHS, terms in the SAME ORDER as
                            the R formula's beta names, e.g.
                            "Base trt Age Visit Base*trt" -- SAS does NOT
                            auto-reorder interactions to the end the way
                            R's terms() does, so get this order right by
                            hand rather than relying on SAS to match R */
    nbeta=,             /* number of fixed-effect coefficients expected,
                            INCLUDING the intercept */
    randomstmt=,        /* one or more full RANDOM statements (with
                            trailing semicolons), e.g.
                            %str(random intercept Visit / subject=subject type=un solution;) */
    indir=%str(data),   /* relative to CWD -- these defaults assume SAS is
                            launched with the `sas/` directory itself as
                            the working directory (matching %include
                            "fitlib.sas" above, which has no path prefix
                            and so also needs CWD=sas/) */
    outdir=%str(.)
  );

  %local t0 t1 nrep;

  proc import datafile="&indir./&example..csv" out=work._raw dbms=csv replace;
    guessingrows=max;
  run;

  proc sort data=work._raw; by rep; run;

  %let t0 = %sysfunc(datetime());

  proc glimmix data=work._raw method=&method;
    by rep;
    %if %length(&classvars) %then %do;
      class &classvars / param=ref ref=first;
    %end;
    model &modelrhs / dist=gamma link=log solution;
    &randomstmt
    ods output ParameterEstimates=work._fe
               CovParms=work._cov
               ConvergenceStatus=work._conv
               FitStatistics=work._fit;
  run;
  quit;

  %let t1 = %sysfunc(datetime());

  /* ---- status/convergence: one row per rep ----
     VERIFY: ConvergenceStatus's exact column names (Status/Reason) and
     the meaning of Status codes against real ODS OUTPUT -- Status=0 is
     documented as "converged with no issues"; anything else is treated
     as a warning here (SAS doesn't cleanly distinguish "warning" from
     "error" the way glmer's tryCatch-based status column does). */
  data work._status(keep=rep status msg);
    set work._conv;
    length status $10 msg $200;
    if Status = 0 then do; status = "clean"; msg = ""; end;
    else do; status = "warning"; msg = Reason; end;
  run;

  proc sql noprint;
    select count(*) into :nrep from work._status;
  quit;

  /* ---- fixed effects: reshape long (one row per rep*effect, in MODEL-
     statement order) to wide (one row per rep, beta1..beta&nbeta), by
     position within each BY group -- mirrors the positional matching
     10_ingest_julia_results.R already uses for Julia's coefnames(). ---- */
  data work._fe_wide(keep=rep beta1-beta%eval(&nbeta));
    set work._fe;
    by rep;
    array betas{&nbeta} beta1-beta%eval(&nbeta);
    retain betas idx;
    if first.rep then do;
      call missing(of betas{*});
      idx = 0;
    end;
    idx + 1;
    if idx <= &nbeta then betas{idx} = estimate;
    if last.rep then output;
  run;

  /* ---- RE sd/corr + dispersion (Scale): CovParms is long (one row per
     rep*covparm). VERIFY the exact CovParm label text produced for your
     RANDOM statement (e.g. "Intercept" for a plain VC random intercept,
     "UN(1,1)"/"UN(2,1)"/"UN(2,2)" for a TYPE=UN 2-effect block, "Scale"
     for the Gamma dispersion parameter) -- the pattern-matching below
     assumes the documented default labels, not verified interactively.
     `singular` is left at 0 (placeholder): GLIMMIX flags a non-positive-
     definite/boundary G matrix in the log, not cleanly in CovParms
     itself, so a real check needs the log, not this table. */
  data work._cov_wide(keep=rep sd1 sd2 corr phi singular);
    set work._cov;
    by rep;
    retain v11 v21 v22 vsingle vsingle2 scale;
    if first.rep then do;
      v11 = .; v21 = .; v22 = .; vsingle = .; vsingle2 = .; scale = .;
    end;
    if upcase(CovParm) = "SCALE" then scale = Estimate;
    else if upcase(CovParm) = "UN(1,1)" then v11 = Estimate;
    else if upcase(CovParm) = "UN(2,1)" then v21 = Estimate;
    else if upcase(CovParm) = "UN(2,2)" then v22 = Estimate;
    else if upcase(CovParm) = "INTERCEPT" and vsingle = . then vsingle = Estimate;
    else if upcase(CovParm) = "INTERCEPT" and vsingle ne . then vsingle2 = Estimate; /* 2nd independent RE term, e.g. report4bb's fyear */
    if last.rep then do;
      singular = 0;
      if v11 ne . and v22 ne . then do;
        sd1 = sqrt(v11); sd2 = sqrt(v22); corr = v21 / (sd1 * sd2);
      end;
      else if vsingle ne . and vsingle2 ne . then do;
        sd1 = sqrt(vsingle); sd2 = sqrt(vsingle2); corr = .; /* independent REs, no correlation param */
      end;
      else if vsingle ne . then do;
        sd1 = sqrt(vsingle); sd2 = .; corr = .;
      end;
      phi = scale;
      output;
    end;
  run;

  /* ---- negll: only meaningful for METHOD=LAPLACE (real -2*logLik on
     the marginal-likelihood scale); left missing for RSPL, see caveat
     above. VERIFY the Fit Statistics row label against real output --
     documented as "-2 Log Likelihood" but confirm the exact Descr text
     ODS OUTPUT emits for your SAS release. */
  %if %upcase(&method) = LAPLACE %then %do;
    data work._negll(keep=rep negll);
      set work._fit;
      where upcase(Descr) = "-2 LOG LIKELIHOOD";
      rename Value = negll;
    run;
  %end;
  %else %do;
    data work._negll(keep=rep negll);
      set work._status(keep=rep);
      negll = .;
    run;
  %end;

  /* ---- assemble + write: time_sec is the whole-batch average (see
     caveat above), identical for every replicate. ---- */
  data work._final(keep=i status singular msg time_sec sd1 sd2 corr phi negll beta1-beta%eval(&nbeta));
    merge work._status work._fe_wide work._cov_wide work._negll;
    by rep;
    i = rep;
    time_sec = (&t1 - &t0) / &nrep;
  run;

  proc export data=work._final
      outfile="&outdir./results_&example._&method..csv"
      dbms=csv replace;
  run;

  %put NOTE: fit_glimmix &example / &method done, &nrep replicates, wrote &outdir./results_&example._&method..csv;

%mend fit_glimmix;
