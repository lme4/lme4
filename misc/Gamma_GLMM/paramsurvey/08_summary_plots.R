## Original six-method (glmmTMB, joint-phi, PIRLS/digamma, PIRLS/moment,
## lme4current, lme4 2.0-6) x five-dataset (epil2 simple/complex/phi>1,
## whale_crate, schizophrenia) summary plots -- preserved as-is (same
## methods/datasets/output filenames as always) as a thin call into the
## shared, parameterized plotlib.R. For a different method/dataset
## subset (e.g. comparing against the Julia MixedModels.jl arm), write a
## new driver script calling make_summary_plots() with its own
## `out_prefix` instead of editing this one -- see
## 11_summary_plots_julia.R for an example.

source(here::here("misc/Gamma_GLMM/paramsurvey/plotlib.R"))

make_summary_plots(
  examples = c("epil2_simple", "epil2_phigt1", "epil2_complex", "report4bb", "schizophrenia"),
  methods  = c("glmmTMB", "jointphi", "pirlsdigamma", "pirlsmoment", "lme4current", "lme4old"),
  negll_ref_method = "glmmTMB",
  negll_diff_exclude = "lme4old",
  negll_diff_outliers = rbind(
    data.frame(example = "epil2_simple", method = "jointphi", threshold = 300),
    ## replicate 417: non-convergence (singular Hessian warning), negll
    ## collapses to a ~1e10 sentinel-like value, not a real likelihood
    data.frame(example = "epil2_phigt1", method = "lme4current", threshold = 1000)
  ),
  out_prefix = ""
)

## "new methods only" pointrange plot -- same five datasets, drops
## PIRLS/fixed-phi (CRAN)/lme4 2.0-6 (the pre-fix baseline, not a fair
## phi-estimator alternative) so the remaining five methods' spread is
## easier to read without the old-CRAN outlier stretching each y-axis.
make_stderr_plot(
  examples = c("epil2_simple", "epil2_phigt1", "epil2_complex", "report4bb", "schizophrenia"),
  methods  = c("glmmTMB", "jointphi", "pirlsdigamma", "pirlsmoment", "lme4current"),
  outfile  = here::here("misc/Gamma_GLMM/paramsurvey/param_summary_stderr_newonly.png")
)
