## Prep step for the "epil2 (simple)" example: single random-intercept
## model y ~ trt + (1|subject) (the reduced form used earlier in this
## investigation, README #1-#2), as a companion to epil2 (complex) --
## same real dataset, much simpler RE structure (no correlation, no
## slope), and (unlike complex) non-singular without any regularizing
## prior, so the reference fit here is plain glmmTMB with no prior.

suppressMessages(library(glmmTMB))
source("/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/param_survey/toolkit.R")

outdir <- "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/param_survey"

data(epil2, package = "glmmTMB")
epil2$subject <- factor(epil2$subject)
epil2_red <- subset(epil2, y > 0)   ## Gamma requires y>0; same base data as epil2 (complex)

form <- y ~ trt + (1 | subject)
family <- Gamma(link = "log")

fit_ref <- glmmTMB(form, data = epil2_red, family = family)
cat("=== reference fit (no prior needed -- already non-singular) ===\n")
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
phi_pretty <- 0.2   ## same pretty value as epil2 (complex), and the disp=0.2 point in README #3

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
sim_seed <- 20260802
set.seed(sim_seed)
sim_y <- lapply(seq_len(B), function(b) obj$simulate(par = pf)$yobs)

cat("\n=== sanity check on simulated responses ===\n")
for (b in seq_len(B)) cat(sprintf("rep %2d: n=%d range=[%.3f, %.2f]\n",
                                    b, length(sim_y[[b]]), min(sim_y[[b]]), max(sim_y[[b]])))

sim_data <- lapply(sim_y, function(y) { d <- epil2_red; d$y <- y; d })

saveRDS(list(name = "epil2_simple", formula = form, family = family,
             base_data = epil2_red, sim_data = sim_data,
             pretty = list(beta = beta_pretty, sd = sd_pretty, corr = NA_real_, phi = phi_pretty),
             fit_ref_summary = list(beta = beta_fit, sd = attr(vc_fit, "stddev"), phi = phi_fit),
             B = B, seed = sim_seed),
        file.path(outdir, "epil2_simple_simdata.rds"))
cat("\nsaved to", file.path(outdir, "epil2_simple_simdata.rds"), "\n")
