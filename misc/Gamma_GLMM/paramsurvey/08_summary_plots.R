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
  negll_diff_outliers = data.frame(example = "epil2_simple", method = "jointphi", threshold = 300),
  out_prefix = ""
)
