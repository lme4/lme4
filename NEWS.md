# CHANGES IN lme4 VERSION 1.1-37.9000 (2025-05-06)

## USER-VISIBLE CHANGES

- lme4 now depends completely on the reformulas package for formula processing etc.; the corresponding functions (expandDoubleVerts, subbars, findbars, nobars, etc.) are no longer implemented within/exported from the package

# CHANGES IN lme4 VERSION 1.1-37 (2025-03-26)

## NEW FEATURES

- simulate.merMod and .simulateFun gain the cluster.rand argument to allow specifying non-Normal cluster random effects in simulations.

- simulation.Rmd, and derived .md and .html files in misc/notes/, trace the flow of the simulation just before this change. Likely of interest only to developers.

## USER-VISIBLE CHANGES

- headline of print and summary output for glmer fits now labels the minimum of the objective function (correctly) as "-2*log(L)" rather than "deviance" (inspired by glmmTMB GH #1156)

# CHANGES IN lme4 VERSION 1.1-36

## NEW FEATURES

- The full (joint) conditional covariance matrix of fixed effects and conditional modes is now available via vcov(fitted_model, full = TRUE)

- Random effects formulas are now processed by code from the reformulas package, meaning that some extended random effects formulas (such as (1|f*g) for crossed random effects) will now work

## USER-VISIBLE CHANGES

- simulate with newparams no longer prints a message if the individual parameter vectors are unnamed

## BUG FIXES

- floating-point errors leading to slightly negative deviances could lead to NaN deviance residuals (Wolfgang Viechtbauer, GH #812)

# CHANGES IN lme4 VERSION 1.1-35.5 (2024-07-03)

## USER-VISIBLE CHANGES

- in predict, the synonyms ReForm, REForm, and REform of re.form are now hard-deprecated, giving an error (after 10 years of soft deprecation)

## OTHER CHANGES

- minor adjustments in test tolerances

- correct Matrix dependency to >= 1.2-3

# CHANGES IN lme4 VERSION 1.1-35.4 (2024-06-19)

## BUG FIXES

- predict(., re.form=...) works in a wider range of cases (GH #691)

- Gamma simulation now uses correct shape parameter (GH #782)

- Avoid triggering RcppEigen UBSAN bug(?) in the case of profiling fixed effects in a merMod object with a single fixed effect parameter (GH #794: lots of help from Dirk Eddelbuettel and Mikael Jagan)

- fix bug in plot methods (cbind'ing zero-length objects)

# CHANGES IN lme4 VERSION 1.1-35.3 (2024-04-16)

## BUG FIXES

- bug-fix for ASAN/memory access problem in CholmodDecomposition (Mikael Jagan)

# CHANGES IN lme4 VERSION 1.1-35.2 (2024-03-28)

## BUG FIXES

- simulate works (again) with re.form=NULL when NA values are present in the data (GH #737, @frousseu)

## USER-VISIBLE CHANGES

- updated tests of upstream Matrix version; should now warn only on ABI incompatibility, not on package version mismatch alone

# CHANGES IN lme4 VERSION 1.1-35.1 (2023-11-05)

## USER-VISIBLE CHANGES

- lFormula and glFormula once again _do_ allow matrix-valued responses (for use in downstream packages like galamm)

# CHANGES IN lme4 VERSION 1.1-35 (2023-11-03)

## NEW FEATURES

- predict.merMod now has a se.fit method, which computes the standard errors of the predictions, conditional on the estimated theta (variance-covariance) parameters

## USER-VISIBLE CHANGES

- using lmer with a matrix-valued response now throws a more informative error message, directing the user to ?refit

# CHANGES IN lme4 VERSION 1.1-34 (2023-07-04)

## BUG FIXES

- summary(<merMod>) now records if correlation was specified explicitly and to what; and its print() method takes it into account; notably summary(<merMod>, correlation=TRUE) will by default print the correlation (matrix of the fixed effects) fixing GH #725

## NEW FEATURES

- refit gains a newweights argument

# CHANGES IN lme4 VERSION 1.1-33 (2023-04-25)

## BUG FIXES

- a boundary check could fail occasionally when large data produced an NA value in a computed gradient; now warns instead (GH #719, Mathias Ambuehl)

- allFit now works better when optimx and dfoptim packages are not installed (GH #724)

- refit reset internal degrees of freedom component incorrectly for REML fits (resulted in incorrect reported REML criteria, but otherwise harmless: side effect of GH #678)

## NEW FEATURES

- dotplot and qqmath methods gain a level argument to set the width of confidence intervals

- dotplot method is now more flexible, using ".v" options (lty.v, col.line.v, lwd.v) to set appearance of vertical lines (Iago Giné Vázquez)

- refit gains a newweights argument (GH #678)

# CHANGES IN lme4 VERSION 1.1-32 (2023-03-14)

## USER-VISIBLE CHANGES

- formatVC() gets a new optional argument corr indicating if correlations or covariances should be used for vector random effects; this corresponds to print(<merMod>, ranef.corr = ...) .  By default, it is FALSE for comp = "Variance", fixing (GH #707).

- qqmath.merMod adds a (useless) data argument for S3 compatibility. Going forward, the id and idLabels arguments should always be specified by name.  We have added code to try to detect/warn when this is not done.

## BUG FIXES

- nobars now retains the environment of its formula argument (GH #713, Mikael Jagan)

# CHANGES IN lme4 VERSION 1.1-31 (2022-11-01)

## BUG FIXES

- confint(fm, <single string>) now works (after years of being broken) again.

- simulating from binomial model with a factor response, when the simulated response contains only a single factor level, now works (Daniel Kennedy)

# CHANGES IN lme4 VERSION 1.1-30 (2022-07-08)

## USER-VISIBLE CHANGES

- nl (term names) component added to output list of mkReTrms (GH #679)

- eliminate partial-matching of eta (for etastart) (GH #686: not actually "user-visible" unless getOption("warnPartialMatchDollar") is TRUE)

- summary method doesn't break for GLMMs other than binomial/Poisson when merDeriv's vcov.glmerMod method is attached (GH #688)

## BUG FIXES

- better handling of simulate(., re.form = NULL) when model frame contains derived components (e.g. offset(), log(x)) (<https://github.com/florianhartig/DHARMa/issues/335>)

- bootMer works with glmmTMB again (broken in 1.1-29)

- maxfun argument to allFit controls max function evaluations for every optimizer type (GH#685)

# CHANGES IN lme4 VERSION 1.1-29 (2022-04-07)

## USER-VISIBLE CHANGES

- prediction with new levels (when not allowed) returns a more informative error message (displays a list of unobserved levels)

## BUG FIXES

- glmer.nb now works when lme4 is not loaded (GH #658, @brgew)

- tests for singularity (check.conv.singular) now run independently of derivative computation (e.g., when calc.derivs=FALSE) (GH #660, @palday)

- influence.merMod now works when data were originally specified as a tibble

- fixed bug in cooks.distance method for influence.merMod (i.e., objects created via influence(fitted_model)) (John Fox) (GH #672)

- predict works for formulas containing . when newdata is specified (GH #653)

- bootMer now correctly inherits control settings from original fit

# CHANGES IN lme4 VERSION 1.1-28 (2022-02-04)

## USER-VISIBLE CHANGES

- construction of interacting factors (e.g. when f1:f2 or f1/f2 occur in random effects terms) is now more efficient for partially crossed designs (doesn't try to create all combinations of f1 and f2) (GH #635 and #636)

- mkNewReTrms is exported

- singular-fit message now refers to help("isSingular") rather than ?isSingular

## TESTS

- fix all.equal(p1,p2,p3) and similar expect_equal() thinkos

- fix some tests only run when lme4:::testLevel() > 1; adapt tests for upcoming Matrix 1.4-1 which has names(diag(<sparse>))

## BUG FIXES

- reOnly preserves environment (GH #654, Mikael Jagan)

- backward-compatibility hooks changed to evaluate at run-time (i.e., in .onLoad()) rather than at build time (GH #649)

- lmList no longer warns when data is a tibble (GH #645)

# CHANGES IN lme4 VERSION 1.1-27.1 (2021-06-22)

## USER-VISIBLE CHANGES

- influence.merMod allows user-specified starting parameters

- cleaned up performance vignette

## BUG FIXES

- cooks.distance now works with objects computed by influence method

- influence.merMod now works with glmer models using nAGQ=0

- predict (with new data) and simulate methods now work for models with >100 levels in a random effect grouping variable (GH #631)

# CHANGES IN lme4 VERSION 1.1-27 (2021-05-15)

## USER-VISIBLE CHANGES

- improvements from Lionel Henry (via https://github.com/lme4/lme4/pull/587) to fix corner cases in data checking; also resolves GH #601 (allFit scoping)

- getME(., "lower") now has names (request of GH #609)

- improved detection of NaN in internal calculations (typically due to underflow/overflow or out-of-bounds linear predictors from non-constraining link functions such as identity-link Gamma models)

- influence.merMod allows parallel computation

- the statmod package is no longer required unless attempting to simulate results from a model with an inverse Gaussian response

## BUG FIXES

- long formulas work better in anova headings (GH #611)

# CHANGES IN lme4 VERSION 1.1-26 (2020-11-30)

## BUG FIXES

- predict, model.frame(.,fixed.only=TRUE) work with variable names containing spaces (GH #605)

- simulate works when original response variable was logical

- densityplot handles partly broken profiles more robustly

## NEW FEATURES

- thpr method for densityplot() (for plotting profiles scaled as densities) gets new arguments

# CHANGES IN lme4 VERSION 1.1-25 (2020-10-23)

- Set more tests to run only if environment variable LME4_TEST_LEVEL>1

# CHANGES IN lme4 VERSION 1.1-24

## USER-VISIBLE CHANGES

- anova() now returns a p-value of NA if the df difference between two models is 0 (implying they are equivalent models) (GH#583, @MetaEntropy)

- speedup in coef() for large models, by skipping conditional variance calculation (Alexander Bauer)

- simulate.formula machinery has changed slightly, for compatibility with the ergm package (Pavel Krivitsky)

- informational messages about (non-)convergence improved (GH #599)

- improved error messages for 0 non-NA cases in data (GH #533)

## NEW FEATURES

- getME(.,"devfun") now works for glmer objects.  Additionally, profile/confint for GLMMs no longer depend on objects in the fitting environment remaining unchanged (GH #589). This change also affects likelihood profiling machinery; results of glmer profiling/CIs may not match results from previous versions exactly.

## BUG FIXES

- improved handling/documentation of glmer.nb controls (GH #556)

- predict works better for gamm4 objects (GH #575)

- resolved some long-standing UBSAN issues (GH #561)

# CHANGES IN lme4 VERSION 1.1-23 (2020-03-06)

- Some PROTECT/UNPROTECT fixes

# CHANGES IN lme4 VERSION 1.1-22

## USER-VISIBLE CHANGES

- prediction now works better for factors with many levels (GH#467, solution by @sihoward)

- minor changes to argument order in [g]lmerControl; default tolerance for convergence checks increased from 0.001 to 0.002 for glmerControl (now consistent with lmerControl)

- lmer(*, family="<fam>") is no longer valid; it had been deprecated since 2013-06.

- lmer(), glmer(), and nlmer() no longer have a formal ... argument.  This defunctifies the use of a sparseX = . argument and will reveal some user errors, where extraneous arguments were previously disregarded.

- In isSingular(x, tol), the default tolerance (tol) has been increased from 1e-5 to 1e-4, the default of check.conv.singular in g?lmerControl().

- for clarity and consistency with base R methods, some column names of anova() output are changed: "Df" becomes "npar", "Chi Df" becomes "Df" (GH #528)

- simulate() now works with inverse-Gaussian models (GH #284 revisited, @nahorp/Florian Hartig)

- single-model mode of anova() now warns about unused arguments in ...  (e.g. type="III")

- default tolerances for nloptwrap/BOBYQA optimizer tightened (xtol_abs and ftol_abs were 1e-6, now 1e-8). (To revert to former tolerances, use control=lmerControl(optimizer="nloptwrap", optCtrl=list(xtol_abs=1e-6, ftol_abs=1e-6)).)

## BUG FIXES

- improved checking for missing data (@lionel-)

- internal checkZrank() should be able to deal with (Matrix package) rankMatrix() returning NA.

- allFit(fm) now works for a model that had an explicit control = lmerControl(..) call.

- internal getStart() now works when model's start was specified as a list, and when called from drop1() on a submodel, fixing GH #521.

- internal function mkdevfun now works even if there is an extraneous getCall function defined in the global environment (GH #535)

- allFit() works even if a variable with symbol i is used somewhere in the original model call (GH #538, reported by Don Cohen); generally more robust

- glmer.nb works even if an alternative version of negative.binomial (other than the one from MASS) is loaded in the workspace (e.g. by the GLMMadaptive package) (GH#516)

- level argument is now honoured by confint(..., type="boot", level=...) (GH #543)

# CHANGES IN lme4 VERSION 1.1-21 (2019-03-05)

## USER-VISIBLE CHANGES

- bootMer now traps and stores messages, warnings, and errors

- bootMer returns an object of class c("bootMer","boot"); new print and confint methods for class bootMer

- small changes to wording of singular-fit messages

# CHANGES IN lme4 VERSION 1.1-20 (2019-02-04)

## USER-VISIBLE CHANGES

- default value for condVar (whether to return conditional variances as part of the ranef.merMod object) is now TRUE

- changed default optimizer to "nloptwrap" (BOBYQA implementation from the nloptr package) for lmer models. To revert to the old default, use control=lmerControl(optimizer="bobyqa")

## BUG FIXES

- adapted tests to work with R-devel's more consistent formula(model.frame(.)) behavior.

# CHANGES IN lme4 VERSION 1.1-19 (2018-11-10)

## NEW FEATURES

- influence measure code from car rolled in (see ?influence.merMod)

- mkReTrm gets new arguments reorder.terms, reorder.vars to control arrangement of RE terms and individual effects with RE terms within model structures

- adding material from the RePsychLing package (on GitHub; see Bates et al 2015 arXiv:1506.04967) to show orthogonal variance components.

- new utility isSingular() function for detecting singular fits

- allFit function/methods have been moved to the main package, rather than being included in an auxiliary source file; computations can (in principle) be done in parallel

## USER-VISIBLE CHANGES

- by default a message is now printed for singular fits (i.e., fits with linear combinations of variance components that are exactly zero)

- as.data.frame.merMod finds conditional variance information stored either as attr(.,"postVar") or attr(.,"condVar") (for glmmTMB compatibility)

- change to defaults of [g]lmerControl to print a message when fits are singular

- post-fitting convergence checks based on estimated gradient and Hessian (see troubleshooting) are no longer performed for (nearly-)singular fits (see isSingular)

# CHANGES IN lme4 VERSION 1.1-18-1 (2018-08-17)

- This is a minor release; the only change is to roll back (unexport) the influence.merMod method, pending resolution of conflicts with the car package

# CHANGES IN lme4 VERSION 1.1-18 (2018-08-16)

## USER-VISIBLE CHANGES

- ranef(.,condVar=TRUE) now works when there are multiple random effects terms per factor

## NEW FEATURES

- rstudent and influence methods are available for merMod objects

- devfun2 function (for generating a deviance function that works on the standard deviation/correlation scale) is now exported

## BUG FIXES

- lmList now obeys its pool argument (instead of always using what currently is the default, GH #476)

# CHANGES IN lme4 VERSION 1.1-17 (2018-04-03)

- This is a maintenance release only (fixes CRAN problems with cross-platform tests and examples)

# CHANGES IN lme4 VERSION 1.1-16 (2018-03-28)

## BUG FIXES

- lmList no longer ignores the subset argument (John Fox)

- fixed several minor issues with predicting when (1) grouping variables have different levels from original model (e.g. missing levels/factor levels not explicitly specified in newdata) or (2) re.form is a subset of the original RE formula and some (unused) grouping variables are omitted from newdata (GH #452, #457)

## USER-VISIBLE CHANGES

- lmList tries harder to collect errors and pass them on as warnings

- documented as.function method (given a merMod object, returns a function that computes the deviance/REML criterion for specified parameters)

- print method for summary.merMod objects no longer collapses small values of the t-statistic to zero

# CHANGES IN lme4 VERSION 1.1-15 (2017-12-21)

## BUG FIXES

- model.frame(., fixed.only=TRUE) now handles models with "non-syntactic" (e.g. space-containing/backtick-delimited) variables in the formula.

- confint(<merMod>) now works again for the default method "profile".

## USER-VISIBLE CHANGES

- exported dotplot.ranef.mer

# CHANGES IN lme4 VERSION 1.1-14 (2017-09-27)

## NEW FEATURES

- added transf argument to dotplot.ranef.mer to allow back-transformation (Ferenci Tamás, GH #134)

- added as.data.frame.ranef.mer convenience method

- user can specify initial value for overdispersion parameter in glmer.nb (Timothy Lau, GH #423)

## BUG FIXES

- fix bug where NAs in fitting data were carried over into predictions on new data (!) (lmwang9527, GH #420)

- fix bug with long terms in models with || notation

- nlmer now respects user-specified lower/upper bounds (GH #432)

- confint.thpr (confint method applied to an already-computed profile now respects "theta_"/"beta_" specifications to return all random-effect or all fixed-effect confidence intervals, respectively.

## DOCUMENTATION IMPROVEMENTS

- document need to export packages and objects to workers when using bootMer with snow

## USER-VISIBLE CHANGES

- improved warning message when using lmerControl() with glmer (GH #415)

- avoid deparsing big data frames when checking data (GH #410)

- pass verbose options to nloptr optimizers when using nloptwrap (previously ignored, with a warning)

- the fl (factor list) component of mkReTrms objects is now returned as a list rather than a data frame

# CHANGES IN lme4 VERSION 1.1-13 (2017-04-18)

## NEW FEATURES

- added prof.scale argument to profile.merMod, documented caveats about using varianceProf/logProf transformation methods for correlation parameters

## BUG FIXES

- suppressed spurious contrast-dropping warning (GH #414)

- fixed bug in confint.lmList4 (GH #26)

- fixed bug when FUN returned an unnamed vector in confint(.,FUN=FUN,method="boot")

- fixed small bug relating to nAGQ0initStep=FALSE

## CRAN-COMPATIBILITY UPDATES

- fixed time stamps on compiled versions of vignettes

# CHANGES IN lme4 VERSION 1.1-12 (2016-04-15)

## USER-VISIBLE CHANGES

- reduced default print precision of fixed-effect correlation matrix in summary.merMod (related to GH #300)

## BUG FIXES

- fixed bug in _de novo_ Gamma-response simulations

# CHANGES IN lme4 VERSION 1.1-11 (2016-02-11)

## USER-VISIBLE CHANGES

- change VarCorr method signature (for compatibility with upstream nlme changes)

## BUG FIXES

- several glmer.nb bugs fixed (generally not changing results, but causing warnings and errors e.g.  during bootstrapping)

- fixes to some lmList bugs (Github #320)

- minor documentation, vignette updates

- minor fix to plot.merMod with id specified

- bootMer now handles separate offset term properly (Github #250)

# CHANGES IN lme4 VERSION 1.1-10 (2015-10-05)

## USER-VISIBLE CHANGES

- updated CITATION file.

## NEW FEATURES

- We export set of about a dozen printing utility functions which are used in our print methods.

- bootMer now allows the use of re.form.

## BUG FIXES

- fixed reordering bug in names of getME(.,"Ztlist") (terms are reordered in decreasing order of the number of levels of the grouping variable, but names were not being reordered)

- fixed issue with simulation when complex forms (such as nested random effects terms) are included in the model (Github #335)

# CHANGES IN lme4 VERSION 1.1-9 (2015-08-20)

## USER-VISIBLE CHANGES

- explicit maxit arguments for various functions (refit, mkGlmerDevfun, ...)

## NEW FEATURES

- terms and formula methods now have random.only options

- getME gains a glmer.nb.theta option.  It is now (an S3) generic with an "merMod" method in lme4 and potentially other methods in dependent packages.

- simulate now works for glmer.nb models (Github #284: idea from @aosmith16)

## BUG FIXES

- prediction and simulation now work when random-effects terms have data-dependent bases (e.g., poly(.) or ns(.) terms) (Github #313, Edgar Gonzalez)

- logLik for glmer.nb models now includes the overdispersion parameter in the parameter count (df attribute)

- lmList handles offsets and weights better

- lots of fixes to glmer.nb (Github #176, #266, #287, #318). *Please note that glmer.nb is still somewhat unstable/under construction.*

## CRAN-COMPATIBILITY UPDATES

- import functions from base packages to pass CRAN checks

- tweak to failing tests on Windows

# CHANGES IN lme4 VERSION 1.1-8 (2015-06-22)

## NEW FEATURES

- getME gains a "Tlist" option (returns a vector of template matrices from which the blocks of Lambda are generated)

- hatvalues method returns the diagonal of the hat matrix of LMMs

- nlminbwrap convenience function allows use of nlminb without going through the optimx package

- as.data.frame.VarCorr.merMod gains an order option that allows the results to be sorted with variances first and covariances last (default) or in lower-triangle order

- allow more flexibility in scales for xyplot.thpr method (John Maindonald)

- models with only random effects of the form 1|f have better starting values for lmer optimization (Gabor Grothendieck)

- glmer now allows a logical vector as the response for binomial models

- anova will now do (sequential) likelihood ratio tests for two or more models including both merMod and glm or lm models (at present, only for GLMMs fitted with the Laplace approximation)

## USER-VISIBLE CHANGES

- deviance() now returns the deviance, rather than half the negative log-likelihood, for GLMMs fitted with Laplace (the behaviour for LMMs and GLMMs fitted with nAGQ>1 has not changed)

- convergence warning and diagnostic test issues are now reported in print and summary methods

- update now (attempts to) re-evaluate the original fit in the environment of its formula (as is done with drop1)

- refit of a nonlinear mixed model fit now throws an error, but this will hopefully change in future releases (related to bug fixes for Github #231)

- lmList now returns objects of class lmList4, to avoid overwriting lmList methods from the recommended nlme package

- names of random effects parameters in confint changed (modified for consistency across methods); oldNames=TRUE (default) gives ".sig01"-style names, oldNames=FALSE gives "sd_(Intercept)|Subject"-style names

- confint(.,method="Wald") result now contains rows for random effects parameters (values set to NA) as well as for fixed-effect parameters

## BUG FIXES

- simulate and predict now work more consistently with different-length data, differing factor levels, and NA values (Github #153, #197, #246, #275)

- refit now works correctly for glmer fits (Github #231)

- fixed bug in family.merMod; non-default links were not retrieved correctly (Alessandro Moscatelli)

- fixed bootMer bug for type=="parametric", use.u=TRUE (Mark Lai)

- gradient scaling for convergence checks now uses the Cholesky factor of the Hessian; while it is more correct, this will lead to some additional (probably false-positive) convergence warnings

- As with lm(), users now get an error for non-finite (Inf, NA, or NaN) values in the response unless na.action is set to exclude or omit them (Github #310)

# CHANGES IN lme4 VERSION 1.1-7 (2014-07-19)

## NEW FEATURES

- the nloptr package is now imported; a wrapper function (nloptwrap) is provided so that lmerControl(optimizer="nloptwrap") is all that's necessary to use nloptr optimizers in the nonlinear optimization stage (the default algorithm is NLopt's implementation of BOBYQA: see ?nloptwrap for examples)

- preliminary implementation of checks for scaling of model matrix columns (see check.scaleX in ?lmerControl)

- beta is now allowed as a synonym for fixef when specifying starting parameters (Github #194)

## USER-VISIBLE CHANGES

- the use of deviance to return the REML criterion is now deprecated; users should use REMLcrit() instead (Github #211)

- changed the default value of check.nobs.vs.rankZ to "ignore" (Github #214)

## BUG FIXES

- change gradient testing from absolute to relative

- fix confint(.,method="boot") to allow/work properly with boot.type values other than "perc" (reported by Alan Zaslavsky)

- allow plot() to work when data are specified in a different environment (reported by Dieter Menne)

- predict and simulate work for matrix-valued predictors (Github #201)

- other simulate bugs (Github #212)

- predict no longer warns spuriously when original response was a factor (Github #205)

- fix memory access issues (Github #200)

# CHANGES IN lme4 VERSION 1.1-6 (2014-04-13)

## BUG FIXES

- change drop1 example to prevent use of old/incompatible pbkrtest versions, for 2.15.3 compatibility

- explicitly require(mlmRev) for tests to prevent cyclic dependency

- bump RcppEigen Imports: requirement from >0.3.1.2.3 to >=0.3.2.0; Rcpp dependency to >= 0.10.5

# CHANGES IN lme4 VERSION 1.1-5 (2014-03-14)

## BUG FIXES

- improved NA handling in simulate and refit

- made internal handling of weights/offset arguments slightly more robust (Github #191)

- handle non-positive-definite estimated fixed effect variance-covariance matrices slightly more generally/robustly (fall back on RX approximation, with a warning, if finite-difference Hessian is non-PD; return NA matrix if RX approximation is also bad)

## MINOR USER-VISIBLE CHANGES

- Added output specifying when Gauss-Hermite quadrature was used to fit the model, and specifying number of GHQ points (Github #190)

# CHANGES IN lme4 VERSION 1.1-4

## BUG FIXES

- Models with prior weights returned an incorrect sigma and deviance (Github issue #155). The deviance bug was only a practical issue in model comparisons, not with inferences given a particular model. Both bugs are now fixed.

- Profiling failed in some cases for models with vector random effects (Github issue #172)

- Standard errors of fixed effects are now computed from the approximate Hessian by default (see the use.hessian argument in vcov.merMod); this gives better (correct) answers when the estimates of the random- and fixed-effect parameters are correlated (Github #47)

## MAJOR USER-VISIBLE CHANGES

- The default optimizer for lmer fits has been switched from "Nelder_Mead" to "bobyqa" because we have generally found the latter to be more reliable.  To switch back to the old behaviour, use control=lmerControl(optimizer="Nelder_Mead").

- Better handling of rank-deficient/overparameterized fixed-effect model matrices; see check.rankX option to [g]lmerControl.  The default value is "message+drop.cols", which automatically drops redundant columns and issues a message (not a warning). (Github #144)

## MINOR USER-VISIBLE CHANGES

- slight changes in convergence checking; tolerances can be specified where appropriate, and some default tolerances have changed (e.g., check.conv.grad)

- improved warning messages about rank-deficiency in X and Z etc. (warnings now try to indicate whether the unidentifiability is in the fixed- or random-effects part of the model)

- predict and simulate now prefer re.form as the argument to specify which random effects to condition on, but allow ReForm, REForm, or REform, giving a message (not a warning) that they are deprecated (addresses Github #170)

- small fixes for printing consistency in models with no fixed effects

- we previously exported a fortify function identical to the one found in ggplot2 in order to be able to define a fortify.merMod S3 method without inducing a dependency on ggplot2.  This has now been unexported to avoid masking ggplot2's own fortify methods; if you want to add diagnostic information to the results of a model, use fortify.merMod explicitly.

- simulate.formula now checks for names associated with the theta and beta parameter vectors. If missing, it prints a message (not a warning); otherwise, it re-orders the parameter vectors to match the internal representation.

- preliminary implementation of a check.scaleX argument in [g]lmerControl that warns about scaling if some columns of the fixed-effect model matrix have large standard deviations (relative to 1, or to each other)

# CHANGES IN lme4 VERSION 1.1-3

## NEW FEATURES

- The gradient and Hessian are now computed via finite differencing after the nonlinear fit is done, and the results are used for additional convergence tests. Control of the behaviour is available through the check.conv.* options in [g]lmerControl. Singular fits (fits with estimated variances of zero or correlations of +/- 1) can also be tested for, although the current default value of the check.conv.singular option is "ignore"; this may be changed to "warning" in the future. The results are stored in @optinfo$derivs.  (Github issue #120; based on code by Rune Christensen.)

- The simulate method will now work to generate simulations "from scratch" by providing a model formula, a data frame holding the predictor variables, and a list containing the values of the model parameters: see ?simulate.merMod. (Github issue #115)

- VarCorr.merMod objects now have an as.data.frame method, converting the list of matrices to a more convenient form for reporting and post-processing. (Github issue #129)

## MINOR USER-VISIBLE CHANGES

- results of fitted(), predict(), and residuals() now have names in all cases (previously results were unnamed, or named only when predicting from new data)

- the anova method now has a refit argument that controls whether objects of class lmerMod should be refitted with ML before producing the anova table.  (Github issues #141, #165; contributed by Henrik Singmann.)

- the print method for VarCorr objects now has a formatter argument for finer control of standard deviation and variance formats

- the optinfo slot now stores slightly more information, including the number of function evaluations ($feval).

- dotplot.ranef.mer now adds titles to sub-plots by default, like qqmath.ranef.mer

## BUG FIXES

- fitted now respects na.action settings (Github issue #149)

- confint(.,method="boot") now works when there are NA values in the original data set (Github issue #158)

- previously, the code stored the results (parameter values, residuals, etc.) based on the _last_ set of parameters evaluated, rather than the optimal parameters.  These were not always the same, but were almost always very close, but some previous results will change slightly (Github issue #166)

# CHANGES IN lme4 VERSION 1.1-0

## MINOR USER-VISIBLE CHANGES

- when using the default method="profile", confint now returns appropriate upper/lower bounds (-1/1 for correlations, 0/Inf for standard deviations) rather than NA when appropriate

## BUG FIXES

- in a previous development version, ranef returned incorrect conditional variances (github issue #148). this is now fixed

# CHANGES IN lme4 VERSION 1.0-6 (2014-02-02)

## BUG FIXES

- prediction now works when new data have fewer factor levels than are present in the original data (Github issue #143, reported by Rune Haubo)

- the existence of a variable "new" in the global environment would mess lme4 up: reported at <http://stackoverflow.com/questions/19801070/error-message-glmer-using-r-what-must-be-a-character-string-or-a-function>

# CHANGES IN lme4 VERSION 1.0-5 (2013-10-24)

## USER-VISIBLE CHANGES

- confint.merMod and vcov.merMod are now exported, for downstream package-author convenience

- the package now depends on Matrix >=1.1-0 and RcppEigen >=0.3.1.2.3

- new rename.response option for refit (see BUG FIXES section)

## BUG FIXES

- eliminated redundant messages about suppressed fixed-effect correlation matrices when p>20

- most inverse-link functions are now bounded where appropriate by .Machine$double.eps, allowing fitting of GLMMs with extreme parameter values

- merMod objects created with refit did not work with update: optional rename.response option added to refit.merMod, to allow this (but the default is still FALSE, for back-compatibility) (reported by A. Kuznetsova)

- fixed buglet preventing on-the-fly creation of index variables, e.g. y~1+(1|rownames(data)) (reported by J. Dushoff)

- predict now works properly for glmer models with basis-creating terms (e.g. poly, ns)

- step sizes determined from fixed effect coefficient standard errors after first state of glmer fitting are now bounded, allowing some additional models to be fitted

# CHANGES IN lme4 VERSION 1.0-4 (2013-09-08)

## BUG FIXES

- refit() now works, again, with lists of length 1, so that e.g. refit(.,simulate(.)) works.  (Reported by Gustaf Granath)

- getME(.,"ST") was returning a list containing the Cholesky factorizations that get repeated in Lambda. But this was inconsistent with what ST represents in lme4.0. This inconsistency has now been fixed and getME(.,"ST") is now consistent with the definition of the ST matrix in lme4.0. See https://github.com/lme4/lme4/issues/111 for more detail. Thanks to Vince Dorie.

- Corrected order of unpacking of standard deviation/correlation components, which affected results from confint(.,method="boot"). (Reported by Reinhold Kliegl)

- fixed a copying bug that made refitML() modify the original model

# CHANGES IN lme4 VERSION 1.0-1 (2013-08-17)

## MINOR USER-VISIBLE CHANGES

- check.numobs.* and check.numlev.* in (g)lmerControl have been changed (from recent development versions) to check.nobs.* and check.nlev.* respectively, and the default values of check.nlev.gtreq.5 and check.nobs.vs.rankZ have been changed to "ignore" and "warningSmall" respectively

- in (g)lmerControl, arguments to the optimizer should be passed as a list called optCtrl, rather than specified as additional (ungrouped) arguments

- the postVar argument to ranef has been changed to the (more sensible) condVar ("posterior variance" was a misnomer, "conditional variance" - short for "variance of the conditional mode" - is preferred)

- the REform argument to predict has been changed to ReForm for consistency

- the tnames function, briefly exported, has been unexported

- getME(.,"cnms") added

- print method for merMod objects is now more terse, and different from summary.merMod

- the objective method for the respMod reference class now takes an optional sigma.sq parameter (defaulting to NULL) to allow calculation of the objective function with a residual variance different from the profiled value (Vince Dorie)

# CHANGES IN lme4 VERSION 1.0-0 (2013-08-01)

## MAJOR USER-VISIBLE CHANGES

- Because the internal computational machinery has changed, results from the newest version of lme4 will not be numerically identical to those from previous versions.  For reasonably well- defined fits, they will be extremely close (within numerical tolerances of 1e-4 or so), but for unstable or poorly-defined fits the results may change, and very unstable fits may fail when they (apparently) succeeded with previous versions. Similarly, some fits may be slower with the new version, although on average the new version should be faster and more stable. More numerical tuning options are now available (see below); non-default settings may restore the speed and/or ability to fit a particular model without an error. If you notice significant or disturbing changes when fitting a model with the new version of lme4, _please notify the maintainers_.

- VarCorr returns its results in the same format as before (as a list of variance-covariance matrices with correlation and stddev attributes, plus a sc attribute giving the residual standard deviation/scale parameter when appropriate), but prints them in a different (nicer) way.

- By default residuals gives deviance (rather than Pearson) residuals when applied to glmer fits (a side effect of matching glm behaviour more closely).

- As another side effect of matching glm behaviour, reported log-likelihoods from glmer models are no longer consistent with those from pre-1.0 lme4, but _are_ consistent with glm; see glmer examples.

## MINOR USER-VISIBLE CHANGES

- More use is made of S3 rather than S4 classes and methods: one side effect is that the nlme and lme4 packages are now much more compatible; methods such as fixef no longer conflict.

- The internal optimizer has changed. [gn]lmer now has an optimizer argument; "Nelder_Mead" is the default for [n]lmer, while a combination of "bobyqa" (an alternative derivative-free method) and "Nelder_Mead" is the default for glmer. To use the nlminb optimizer as in the old version of lme4, you can use optimizer="optimx" with control=list(method="nlminb") (you will need the optimx package to be installed and loaded). See lmerControl for details.

- Families in GLMMs are no longer restricted to built-in/hard- coded families; any family described in family, or following that design, is usable (although there are some hard-coded families, which will be faster).

- [gn]lmer now produces objects of class merMod rather than class mer as before.

- the structure of the Zt (transposed random effect design matrix) as returned by getME(.,"Zt"), and the corresponding order of the random effects vector (getME(.,"u")) have changed. To retrieve Zt in the old format, use do.call(Matrix::rBind,getME(.,"Ztlist")).

- the package checks input more thoroughly for non-identifiable or otherwise problematic cases: see lmerControl for fine control of the test behaviour.

## NEW FEATURES

- A general-purpose getME accessor method allows extraction of a wide variety of components of a mixed-model fit. getME also allows a vector of objects to be returned as a list of mixed-model components. This has been backported to be compatible with older versions of lme4 that still produce mer objects rather than merMod objects. However, backporting is incomplete; some objects are only extractable in newer versions of lme4.

- Optimization information (convergence codes, warnings, etc.) is now stored in an @optinfo slot.

- bootMer provides a framework for obtaining parameter confidence intervals by parametric bootstrapping.

- plot.merMod provides diagnostic plotting methods similar to those from the nlme package (although missing augPred).

- A predict.merMod method gives predictions; it allows an effect-specific choice of conditional prediction or prediction at the population level (i.e., with random effects set to zero).

- Likelihood profiling for lmer and glmer results (see link{profile-methods}).

- Confidence intervals by likelihood profiling (default), parametric bootstrap, or Wald approximation (fixed effects only): see confint.merMod

- nAGQ=0, an option to do fast (but inaccurate) fitting of GLMMs.

- Using devFunOnly=TRUE allows the user to extract a deviance function for the model, allowing further diagnostics/customization of model results.

- The internal structure of [gn]lmer is now more modular, allowing finer control of the different steps of argument checking; construction of design matrices and data structures; parameter estimation; and construction of the final merMod object (see ?modular).

- the formula, model.frame, and terms methods return full versions (including random effect terms and input variables) by default, but a fixed.only argument allows access to the fixed effect submodel.

## EXPERIMENTAL FEATURES

- glmer.nb provides an embryonic negative binomial fitting capability.

## STILL NON-EXISTENT FEATURES

- Adaptive Gaussian quadrature (AGQ) is not available for multiple and/or non-scalar random effects.

- Posterior variances of conditional models for non-scalar random effects.

- Standard errors for predict.merMod results.

- Automatic MCMC sampling based on the fit turns out to be very difficult to implement in a way that is really broadly reliable and robust; mcmcsamp will not be implemented in the near future. See pvalues for alternatives.

- "R-side" structures (within-block correlation and heteroscedasticity) are not on the current timetable.

## BUG FIXES

- In a development version, prior weights were not being used properly in the calculation of the residual standard deviation, but this has been fixed.  Thanks to Simon Wood for pointing this out.

- In a development version, the step-halving component of the penalized iteratively reweighted least squares algorithm was not working, but this is now fixed.

- In a development version, square RZX matrices would lead to a pwrssUpdate did not converge in 30 iterations error. This has been fixed by adding an extra column of zeros to RZX.

## DEPRECATED AND DEFUNCT

- Previous versions of lme4 provided the mcmcsamp function, which efficiently generated a Markov chain Monte Carlo sample from the posterior distribution of the parameters, assuming flat (scaled likelihood) priors. Due to difficulty in constructing a version of mcmcsamp that was reliable even in cases where the estimated random effect variances were near zero (e.g. <https://stat.ethz.ch/pipermail/r-sig-mixed-models/2009q4/003115.html>), mcmcsamp has been withdrawn (or more precisely, not updated to work with lme4 versions >=1.0).

- Calling glmer with the default gaussian family redirects to lmer, but this is deprecated (in the future glmer(...,family="gaussian") may fit a LMM using the penalized iteratively reweighted least squares algorithm). Please call lmer directly.

- Calling lmer with a family argument redirects to glmer; this is deprecated. Please call glmer directly.

# CHANGES IN lme4 VERSION 0.999375-16 (2008-06-23)

## MAJOR USER-VISIBLE CHANGES

- The underlying algorithms and representations for all the mixed-effects models fit by this package have changed - for the better, we hope. The class "mer" is a common mixed-effects model representation for linear, generalized linear, nonlinear and generalized nonlinear mixed-effects models.

- ECME iterations are no longer used at all, nor are analytic gradients. Components named 'niterEM', 'EMverbose', or 'gradient' can be included in the 'control' argument to lmer(), glmer() or nlmer() but have no effect.

- PQL iterations are no longer used in glmer() and nlmer().  Only the Laplace approximation is currently available. AGQ, for certain classes of GLMMs or NLMMs, is being added.

- The 'method' argument to lmer(), glmer() or nlmer() is deprecated. Use the 'REML = FALSE' in lmer() to obtain ML estimates. Selection of AGQ in glmer() and nlmer() will be controlled by the argument 'nAGQ', when completed.

## NEW FEATURES

- The representation of mixed-effects models has been dramatically changed to allow for smooth evaluation of the objective as the variance-covariance matrices for the random effects approach singularity. Beta testers found this representation to be more robust and usually faster than previous versions of lme4.

- The mcmcsamp function uses a new sampling method for the variance-covariance parameters that allows recovery from singularity. The update is not based on a sample from the Wishart distribution. It uses a redundant parameter representation and a linear least squares update.

- CAUTION: Currently the results from mcmcsamp look peculiar and are probably incorrect. I hope it is just a matter of my omitting a scaling factor but I have seen patterns such as the parameter estimate for some variance-covariance parameters being the maximum value in the chain, which is highly unlikely.

- The 'verbose' argument to lmer(), glmer() and nlmer() can be used instead of 'control = list(msVerbose = TRUE)'.

