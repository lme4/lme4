## Prep step for the crossed-RE extension of the log-link-gaussian survey.
## Simulates ONE full crossed factorial dataset (16 x 16 levels x 3 reps,
## n=768) with two random intercepts (Rail1, Rail2) per replicate, then
## draws eight fixed subsamples of that same simulated data (same
## underlying RE realizations, only the visible rows differ), fit with
## the crossed model ~1+(1|Rail1)+(1|Rail2) except the two "oneway*" cases
## (single-RE control, no crossing, fit with ~1+(1|Rail1)):
##   - "structured" / "structured5" / "structured6": first k levels of
##     each factor (k=4/5/6), full crossing preserved (k^2 x 3 rows)
##   - "random" / "random75" / "random108": uniform subsamples of the
##     full 768 rows (n=48/75/108), crossing not guaranteed
##   - "oneway": all 16 levels of Rail1, only the first level of Rail2
##     (16 x 1 x 3 = 48) -- q is unambiguously nlevels(Rail1)=16
##   - "oneway10": first 10 levels of Rail1, only the first level of
##     Rail2 (10 x 1 x 3 = 30) -- same control idea at a smaller q=10
## These sizes let structured/random pairs be compared at matched n
## (48, 75, 108) while the crossed-vs-oneway contrast at similar q
## (16 vs the crossed cases' q1+q2 in the same ballpark) tests whether
## the effective degrees-of-freedom correction for two crossed RE terms
## is just q1+q2 or something else (see README's GH #936 section).
## RE SDs (1 on the log scale) and a tiny residual SD (0.01) are
## deliberately chosen to keep fits away from the singular boundary (an
## earlier attempt with RE sd=3 gave travel values spanning 7+ orders of
## magnitude and broke glmer's PIRLS Cholesky step outright --
## "Downdated VtV is not positive definite" -- on every replicate, even
## though glmmTMB/mgcv tolerated it; sd=1 keeps the dynamic range sane).

suppressMessages(library(glmmTMB))
## RTMB backend doesn't implement Gamma yet (as of this glmmTMB dev
## version); force the legacy backend so gaussian and Gamma runs are
## on the same numerical footing
glmmTMB:::useRTMB(FALSE)

outdir <- here::here("misc", "Gamma_GLMM", "paramsurvey_loggaussian")

## args: [1] B (default 10), [2] family name (default "gaussian" -- also
## supports "Gamma", both with link="log"; the q_eff mechanism is a
## dispersion-scale property that should be family-generic, this is the
## flexible hook for testing that on other families without duplicating
## the whole pipeline)
args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args) >= 1) as.integer(args[[1]]) else 10
family_name <- if (length(args) >= 2) args[[2]] else "gaussian"
family <- switch(family_name,
                  gaussian = gaussian(link = "log"),
                  Gamma = Gamma(link = "log"),
                  stop("unrecognized family_name: ", family_name))

form <- travel ~ 1 + (1 | Rail1) + (1 | Rail2)

beta_pretty <- 4.1
sd_pretty <- 1
sigma_pretty <- 0.01

design_seed <- 20260807
set.seed(design_seed)

full_design <- expand.grid(Rail1 = factor(1:16), Rail2 = factor(1:16), rep = 1:3)
cat(sprintf("full design: %d rows, %d x %d levels (Rail1 x Rail2)\n",
            nrow(full_design), nlevels(full_design$Rail1), nlevels(full_design$Rail2)))

## fixed row-index subsamples, drawn once (same subsample reused for all B
## replicates -- only the response is resimulated, per 01_prep_rail.R's
## convention). Evaluated in this order so the pre-existing "random" case
## keeps drawing the same first sample() call as before (reproducible
## across the earlier 3-case run).
case_idx <- list(
  structured  = with(full_design, which(Rail1 %in% levels(Rail1)[1:4] & Rail2 %in% levels(Rail2)[1:4])),
  random      = sample(nrow(full_design), 48),
  oneway      = with(full_design, which(Rail2 == levels(Rail2)[1])),
  structured5 = with(full_design, which(Rail1 %in% levels(Rail1)[1:5] & Rail2 %in% levels(Rail2)[1:5])),
  random75    = sample(nrow(full_design), 75),
  structured6 = with(full_design, which(Rail1 %in% levels(Rail1)[1:6] & Rail2 %in% levels(Rail2)[1:6])),
  random108   = sample(nrow(full_design), 108),
  oneway10    = with(full_design, which(Rail1 %in% levels(Rail1)[1:10] & Rail2 == levels(Rail2)[1]))
)

designs <- lapply(case_idx, function(idx) droplevels(full_design[idx, ]))

for (nm in names(designs)) {
  d <- designs[[nm]]
  cat(sprintf("%-12s: %3d rows, %2d x %2d levels (Rail1 x Rail2)\n",
              nm, nrow(d), nlevels(d$Rail1), nlevels(d$Rail2)))
}

sim_seed <- 20260807

## glmmTMB's betadisp is on the log scale of a family-specific dispersion
## parameter: for gaussian that's sigma itself, but for Gamma it's the
## shape parameter (= 1/sigma^2, sigma being the CV) -- verified
## empirically (fit a known-shape simulated Gamma sample, checked
## exp(betadisp) against the true shape) since the two families don't
## share a convention here.
betadisp_pretty <- switch(family_name,
                           gaussian = log(sigma_pretty),
                           Gamma = -2 * log(sigma_pretty))
newparams <- list(beta = beta_pretty, betadisp = betadisp_pretty, theta = log(rep(sd_pretty, 2)))

## simulate B response vectors on the FULL design once -- every case above
## is a subsample of these same simulated realizations
sim_y_full <- simulate_new(~ 1 + (1 | Rail1) + (1 | Rail2), nsim = B, seed = sim_seed,
                            family = family, newdata = full_design, newparams = newparams)

mk_sim_data <- function(idx, design) lapply(sim_y_full, function(y) { d <- design; d$travel <- y[idx]; d })

sim_data_list <- Map(mk_sim_data, case_idx, designs)

for (nm in names(sim_data_list)) {
  sd_list <- sim_data_list[[nm]]
  cat(sprintf("\n=== sanity check on simulated responses (%s) ===\n", nm))
  for (b in seq_len(B)) cat(sprintf("rep %2d: n=%d range=[%.3f, %.2f]\n",
                                      b, nrow(sd_list[[b]]), min(sd_list[[b]]$travel), max(sd_list[[b]]$travel)))
}

cases <- Map(function(design, sim_data) list(design = design, sim_data = sim_data), designs, sim_data_list)

## ---- random-slope control case: correlated (1+x|Rail1) ----
## Reuses "oneway10"'s 30-row design (first 10 Rail1 levels, Rail2
## fixed/irrelevant here) and adds a small-variance covariate x. A
## single RE term (K=1), but with 2 columns/level (intercept, slope)
## instead of 1 -- tests whether q_eff = sum(q_k) - (K-1) generalizes to
## non-scalar terms, since here it predicts q_eff = q = 2*10 = 20 (no
## extra -1, because K=1). True intercept/slope SD = 1, correlation =
## 0.1.
design_randomslope10 <- designs$oneway10
design_randomslope10$x <- rnorm(nrow(design_randomslope10), 0, 0.1)
cat(sprintf("%-12s: %3d rows, %2d levels (Rail1), x added\n",
            "randomslope10", nrow(design_randomslope10), nlevels(design_randomslope10$Rail1)))

corr_rs_pretty <- 0.1
## glmmTMB's "us" (unstructured) covariance parameterization is NOT the
## same convention as lme4's theta (a raw relative-Cholesky-factor of
## the covariance itself): the first n elements are log-SDs, the
## remaining n(n-1)/2 are a *scaled* Cholesky term of the correlation
## matrix, not the correlation directly -- see
## https://cran.r-project.org/web/packages/glmmTMB/vignettes/covstruct.html#mappings.
## For sd1=sd2=1, correlation=rho: theta = c(0, 0, rho/sqrt(1-rho^2)).
theta_rs <- c(log(sd_pretty), log(sd_pretty), corr_rs_pretty / sqrt(1 - corr_rs_pretty^2))
newparams_rs <- list(beta = beta_pretty, betadisp = betadisp_pretty, theta = theta_rs)

sim_y_randomslope10 <- simulate_new(~ 1 + (1 + x | Rail1), nsim = B, seed = sim_seed,
                                     family = family, newdata = design_randomslope10,
                                     newparams = newparams_rs)
sim_data_randomslope10 <- lapply(sim_y_randomslope10, function(y) { d <- design_randomslope10; d$travel <- y; d })

cat("\n=== sanity check on simulated responses (randomslope10) ===\n")
for (b in seq_len(B)) cat(sprintf("rep %2d: n=%d range=[%.3f, %.2f]\n",
                                    b, nrow(sim_data_randomslope10[[b]]), min(sim_data_randomslope10[[b]]$travel), max(sim_data_randomslope10[[b]]$travel)))

cases$randomslope10 <- list(design = design_randomslope10, sim_data = sim_data_randomslope10)

## ---- nested-RE control case: ~1 + (1|Rail1/Rail2) ----
## Same 75-row 5x5x3 design as "structured5" (full crossing), but fit
## as NESTED (Rail1, Rail1:Rail2) rather than crossed (Rail1, Rail2):
## q1=nlevels(Rail1)=5, q2=nlevels(Rail1:Rail2)=25, K=2 terms. Tests
## whether q_eff=sum(q_k)-(K-1) depends on crossing vs. nesting -- the
## identifiability argument behind the -1 didn't depend on how the
## terms relate to each other, only that each is its own intercept-
## bearing term, so the same formula (q_eff=5+25-1=29) is predicted to
## hold here too. Same true SD=1 for both terms as the original crossed
## case.
sim_y_nested5 <- simulate_new(~ 1 + (1 | Rail1 / Rail2), nsim = B, seed = sim_seed,
                               family = family, newdata = designs$structured5, newparams = newparams)
sim_data_nested5 <- lapply(sim_y_nested5, function(y) { d <- designs$structured5; d$travel <- y; d })

cat("\n=== sanity check on simulated responses (nested5) ===\n")
for (b in seq_len(B)) cat(sprintf("rep %2d: n=%d range=[%.3f, %.2f]\n",
                                    b, nrow(sim_data_nested5[[b]]), min(sim_data_nested5[[b]]$travel), max(sim_data_nested5[[b]]$travel)))

cases$nested5 <- list(design = designs$structured5, sim_data = sim_data_nested5)

## ---- crossed random-slopes case: (1+x|Rail1) + (1+x|Rail2) ----
## Same 75-row 5x5x3 design plus a covariate x, but now TWO independent
## correlated random-slope terms instead of one -- combines the
## crossed-term (-1) mechanism with the multi-column-term (spherical-
## dimension) mechanism at once: q1=2*5=10, q2=2*5=10, K=2, predicted
## q_eff=10+10-1=19. Same per-term parameters as "randomslope10" (SD=1,
## SD=1, corr=0.1) for both terms -- theta is simply that term's theta
## repeated once per term (the two terms are independent of each other
## by construction, so there's no cross-term parameter to add).
design_crossedslopes5 <- designs$structured5
design_crossedslopes5$x <- rnorm(nrow(design_crossedslopes5), 0, 0.1)
newparams_cs <- list(beta = beta_pretty, betadisp = betadisp_pretty, theta = rep(theta_rs, 2))

sim_y_crossedslopes5 <- simulate_new(~ 1 + (1 + x | Rail1) + (1 + x | Rail2), nsim = B, seed = sim_seed,
                                      family = family, newdata = design_crossedslopes5,
                                      newparams = newparams_cs)
sim_data_crossedslopes5 <- lapply(sim_y_crossedslopes5, function(y) { d <- design_crossedslopes5; d$travel <- y; d })

cat("\n=== sanity check on simulated responses (crossedslopes5) ===\n")
for (b in seq_len(B)) cat(sprintf("rep %2d: n=%d range=[%.3f, %.2f]\n",
                                    b, nrow(sim_data_crossedslopes5[[b]]), min(sim_data_crossedslopes5[[b]]$travel), max(sim_data_crossedslopes5[[b]]$travel)))

cases$crossedslopes5 <- list(design = design_crossedslopes5, sim_data = sim_data_crossedslopes5)

outfile <- file.path(outdir, sprintf("crossed_simdata_%s.rds", family_name))
saveRDS(list(name = "crossed", formula = form, family = family, family_name = family_name,
             pretty = list(beta = beta_pretty, sd1 = sd_pretty, sd2 = sd_pretty, sigma = sigma_pretty,
                            corr = corr_rs_pretty),
             cases = cases,
             B = B, seed = sim_seed, design_seed = design_seed),
        outfile)
cat("\nsaved to", outfile, "\n")
