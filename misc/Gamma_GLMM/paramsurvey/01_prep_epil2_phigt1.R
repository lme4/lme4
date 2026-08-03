## Prep step for the "epil2 (phi>1)" example: identical to epil2 (simple)
## (same real base dataset, same y ~ trt + (1|subject) formula, same
## beta_pretty/sd_pretty) except phi_pretty is set to 1.2 instead of 0.2
## -- a synthetic stress test of the phi>1 (CV>1) regime, which none of
## the real reference datasets in this survey reach (see discussion in
## misc/Gamma_GLMM/README_Gamma_GLMMs.md and the chat search for a real
## phi>1 example -- closest confirmed real case found was Gojanovic et al.
## 2026's reef-fish "max distance ventured from shelter" GLMM, phi ~= 2.95,
## but its random effects collapse to ~0, unlike this survey's other
## examples).

suppressMessages(library(glmmTMB))
source("/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/param_survey/toolkit.R")

outdir <- "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/param_survey"

data(epil2, package = "glmmTMB")
epil2$subject <- factor(epil2$subject)
epil2_red <- subset(epil2, y > 0)   ## Gamma requires y>0; same base data as epil2 (simple/complex)

form <- y ~ trt + (1 | subject)
family <- Gamma(link = "log")

fit_ref <- glmmTMB(form, data = epil2_red, family = family)
cat("=== reference fit (same as epil2 (simple); only phi_pretty changes below) ===\n")
print(summary(fit_ref))
stopifnot(!performance::check_singularity(fit_ref))

beta_fit <- fixef(fit_ref)$cond
cat("\nfitted beta:\n"); print(beta_fit)
beta_pretty <- c(`(Intercept)` = 1.9, trtprogabide = -0.2)
stopifnot(identical(names(beta_pretty), names(beta_fit)))

vc_fit <- VarCorr(fit_ref)$cond$subject
cat("\nfitted RE sd:\n"); print(attr(vc_fit, "stddev"))
sd_pretty <- c(`(Intercept)` = 0.8)

phi_fit <- sigma(fit_ref)^2
cat("\nfitted dispersion (phi):", phi_fit, "\n")
phi_pretty <- 1.2   ## artificially set >1 (CV = sqrt(1.2) ~= 1.10); real fitted
                     ## value here is ~0.2, same as epil2 (simple)/(complex)

cat("\n=== pretty (rounded) true parameters ===\n")
print(beta_pretty); cat("sd:", sd_pretty, " phi:", phi_pretty, "\n")

pre <- glmmTMB(form, data = epil2_red, family = family, doFit = FALSE)
obj <- glmmTMB:::fitTMB(pre, doOptim = FALSE)
pf <- obj$env$last.par.best
stopifnot(sum(names(pf) == "beta") == 2, sum(names(pf) == "betadisp") == 1,
          sum(names(pf) == "theta") == 1)

pf[names(pf) == "beta"] <- unname(beta_pretty)
pf[names(pf) == "betadisp"] <- -log(phi_pretty)
pf[names(pf) == "theta"] <- log(sd_pretty[1])

args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args) >= 1) as.integer(args[[1]]) else 10
sim_seed <- 20260803
set.seed(sim_seed)
sim_y <- lapply(seq_len(B), function(b) obj$simulate(par = pf)$yobs)

cat("\n=== sanity check on simulated responses ===\n")
for (b in seq_len(B)) cat(sprintf("rep %2d: n=%d range=[%.3f, %.2f]\n",
                                    b, length(sim_y[[b]]), min(sim_y[[b]]), max(sim_y[[b]])))

sim_data <- lapply(sim_y, function(y) { d <- epil2_red; d$y <- y; d })

saveRDS(list(name = "epil2_phigt1", formula = form, family = family,
             base_data = epil2_red, sim_data = sim_data,
             pretty = list(beta = beta_pretty, sd = sd_pretty, corr = NA_real_, phi = phi_pretty),
             fit_ref_summary = list(beta = beta_fit, sd = attr(vc_fit, "stddev"), phi = phi_fit),
             B = B, seed = sim_seed),
        file.path(outdir, "epil2_phigt1_simdata.rds"))
cat("\nsaved to", file.path(outdir, "epil2_phigt1_simdata.rds"), "\n")
