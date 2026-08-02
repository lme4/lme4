## Prep step for the "epil2 (complex)" example: fit the full model from
## ?glmmTMB::epil2 (y ~ Base*trt + Age + Visit + (Visit|subject)) to the
## real data with a regularizing ranef prior (normal(0,3)) to get a
## non-singular reference fit, round its parameters to "pretty" values
## near the real estimates, and simulate B new datasets from those pretty
## parameters using the real design (same covariates/grouping, new
## response only). The prior is used ONLY here, to pin down realistic
## non-degenerate true parameters -- it is NOT used when fitting any of
## the four methods to the simulated data later (per-replicate fits are
## expected to include a substantial fraction of singular results, which
## is itself part of what this survey is measuring).

suppressMessages(library(glmmTMB))
source("/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/param_survey/toolkit.R")

outdir <- "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/param_survey"

data(epil2, package = "glmmTMB")
epil2$subject <- factor(epil2$subject)
epil2_red <- subset(epil2, y > 0)   ## Gamma requires y>0; drops 23/236 obs (precedent: misc/Gamma_GLMM/nAGQ_epil2.R)

form <- y ~ Base * trt + Age + Visit + (Visit | subject)
family <- Gamma(link = "log")
prior <- data.frame(prior = "normal(0,3)", class = "ranef")

fit_ref <- glmmTMB(form, data = epil2_red, family = family, priors = prior)
cat("=== reference fit (with regularizing prior) ===\n")
print(summary(fit_ref))
stopifnot(!performance::check_singularity(fit_ref))

## round the reference fit's estimates to "pretty" values
beta_fit <- fixef(fit_ref)$cond
cat("\nfitted beta:\n"); print(beta_fit)
beta_pretty <- c(`(Intercept)` = -0.7, Base = 0.85, trtprogabide = -0.5,
                  Age = 0.3, Visit = -0.1, `Base:trtprogabide` = 0.15)
stopifnot(identical(names(beta_pretty), names(beta_fit)))

vc_fit <- VarCorr(fit_ref)$cond$subject
cat("\nfitted RE sd/corr:\n"); print(attr(vc_fit, "stddev")); print(attr(vc_fit, "correlation"))
sd_pretty <- c(`(Intercept)` = 0.4, Visit = 0.15)
rho_pretty <- -0.1

phi_fit <- sigma(fit_ref)^2
cat("\nfitted dispersion (phi):", phi_fit, "\n")
phi_pretty <- 0.2   ## matches the disp=0.2 point already characterized in README #3

cat("\n=== pretty (rounded) true parameters ===\n")
print(beta_pretty); cat("sd:", sd_pretty, " corr:", rho_pretty, " phi:", phi_pretty, "\n")

## build the TMB obj template (unfitted; doOptim=FALSE) against the real
## design/covariates, then overwrite its beta/betadisp/theta slots with
## the pretty parameters. The "b" (random effect) and "theta" slots in the
## template are irrelevant beyond providing the right length/positions --
## obj$simulate() redraws the random effects fresh from the RE prior
## implied by theta each call (verified: two calls with identical `par`
## give different output).
pre <- glmmTMB(form, data = epil2_red, family = family, doFit = FALSE)
obj <- glmmTMB:::fitTMB(pre, doOptim = FALSE)
pf <- obj$env$last.par.best
stopifnot(sum(names(pf) == "beta") == 6, sum(names(pf) == "betadisp") == 1,
          sum(names(pf) == "theta") == 3)

pf[names(pf) == "beta"] <- unname(beta_pretty)
pf[names(pf) == "betadisp"] <- -log(phi_pretty)
pf[names(pf) == "theta"] <- c(log(sd_pretty[1]), log(sd_pretty[2]), corr_to_glmmTMB_theta(rho_pretty))

args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args) >= 1) as.integer(args[[1]]) else 10
sim_seed <- 20260802
set.seed(sim_seed)
sim_y <- lapply(seq_len(B), function(b) obj$simulate(par = pf)$yobs)

cat("\n=== sanity check on simulated responses ===\n")
for (b in seq_len(B)) cat(sprintf("rep %2d: n=%d range=[%.3f, %.2f]\n",
                                    b, length(sim_y[[b]]), min(sim_y[[b]]), max(sim_y[[b]])))

sim_data <- lapply(sim_y, function(y) { d <- epil2_red; d$y <- y; d })

saveRDS(list(name = "epil2_complex", formula = form, family = family,
             base_data = epil2_red, sim_data = sim_data,
             pretty = list(beta = beta_pretty, sd = sd_pretty, corr = rho_pretty, phi = phi_pretty),
             fit_ref_summary = list(beta = beta_fit, sd = attr(vc_fit, "stddev"),
                                     corr = attr(vc_fit, "correlation")[1, 2], phi = phi_fit),
             B = B, seed = sim_seed),
        file.path(outdir, "epil2_complex_simdata.rds"))
cat("\nsaved to", file.path(outdir, "epil2_complex_simdata.rds"), "\n")
