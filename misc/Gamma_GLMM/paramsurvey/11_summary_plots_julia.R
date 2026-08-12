## Comparison plots for the Julia arm (MixedModels.jl, pa/dispersion-
## again) against the three R methods it's most directly comparable to
## (glmmTMB reference, joint-phi, and lme4's current C++ moment-phi fix
## -- the R-level PIRLS-digamma/moment prototypes and old lme4 2.0-6
## aren't part of this comparison). Uses a distinct `out_prefix` so
## these never collide with 08_summary_plots.R's original PNGs.

source(here::here("misc/Gamma_GLMM/paramsurvey/plotlib.R"))

make_summary_plots(
  examples = c("epil2_simple", "report4bb", "schizophrenia", "epil2_complex"),
  methods  = c("glmmTMB", "jointphi", "lme4current", "juliaMixedModels"),
  negll_ref_method = "glmmTMB",
  out_prefix = "julia_compare_"
)
