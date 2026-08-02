## Prep step for the "schizophrenia" example: bundled directly in lme4
## itself (data(schizophrenia), man/schizophrenia.Rd's own example model),
## fit as a Gamma GLMM. A genuinely different regime from the other three
## examples: much larger (1603 obs, 437 subjects, single random intercept)
## -- non-singular without any regularizing prior.

suppressMessages(library(glmmTMB))
suppressMessages(library(lme4))
source("/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/param_survey/toolkit.R")

outdir <- "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/param_survey"

data(schizophrenia)
schizophrenia$id <- factor(schizophrenia$id)

form <- imps79 ~ TxDrug * Week + (1 | id)
family <- Gamma(link = "log")

fit_ref <- glmmTMB(form, data = schizophrenia, family = family)
cat("=== reference fit (no prior needed -- already non-singular) ===\n")
print(summary(fit_ref))
stopifnot(!performance::check_singularity(fit_ref))

beta_fit <- fixef(fit_ref)$cond
cat("\nfitted beta:\n"); print(beta_fit)
beta_pretty <- c(`(Intercept)` = 1.65, TxDrug = 0, Week = -0.04, `TxDrug:Week` = -0.06)
stopifnot(identical(names(beta_pretty), names(beta_fit)))

vc_fit <- VarCorr(fit_ref)$cond$id
cat("\nfitted RE sd:\n"); print(attr(vc_fit, "stddev"))
sd_pretty <- c(id = 0.2)

phi_fit <- sigma(fit_ref)^2
cat("\nfitted dispersion (phi):", phi_fit, "\n")
phi_pretty <- 0.08

cat("\n=== pretty (rounded) true parameters ===\n")
print(beta_pretty); cat("sd:", sd_pretty, " phi:", phi_pretty, "\n")

pre <- glmmTMB(form, data = schizophrenia, family = family, doFit = FALSE)
obj <- glmmTMB:::fitTMB(pre, doOptim = FALSE)
pf <- obj$env$last.par.best
stopifnot(sum(names(pf) == "beta") == 4, sum(names(pf) == "betadisp") == 1,
          sum(names(pf) == "theta") == 1)

pf[names(pf) == "beta"] <- unname(beta_pretty)
pf[names(pf) == "betadisp"] <- -log(phi_pretty)
pf[names(pf) == "theta"] <- log(sd_pretty[1])

args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args) >= 1) as.integer(args[[1]]) else 10
sim_seed <- 20260802
set.seed(sim_seed)
sim_y <- lapply(seq_len(B), function(b) obj$simulate(par = pf)$yobs)

cat("\n=== sanity check on simulated responses ===\n")
for (b in seq_len(B)) cat(sprintf("rep %2d: n=%d range=[%.3f, %.2f]\n",
                                    b, length(sim_y[[b]]), min(sim_y[[b]]), max(sim_y[[b]])))

sim_data <- lapply(sim_y, function(y) { d <- schizophrenia; d$imps79 <- y; d })

saveRDS(list(name = "schizophrenia", formula = form, family = family,
             base_data = schizophrenia, sim_data = sim_data,
             pretty = list(beta = beta_pretty, sd = sd_pretty, corr = NA_real_, phi = phi_pretty),
             fit_ref_summary = list(beta = beta_fit, sd = attr(vc_fit, "stddev"), phi = phi_fit),
             B = B, seed = sim_seed),
        file.path(outdir, "schizophrenia_simdata.rds"))
cat("\nsaved to", file.path(outdir, "schizophrenia_simdata.rds"), "\n")
