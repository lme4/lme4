## Prep step: fit glmmTMB (gaussian(link="log")) to the real nlme::Rail
## data -- the exact example from GH #936 -- round the fitted parameters
## to nearby "pretty" values, then use glmmTMB::simulate_new() to draw B
## fresh datasets from those pretty parameters using Rail's real design
## (6 rails x 3 obs each, n=18).

suppressMessages(library(glmmTMB))

outdir <- here::here("misc", "Gamma_GLMM", "paramsurvey_loggaussian")

data("Rail", package = "nlme")

form <- travel ~ 1 + (1 | Rail)
family <- gaussian(link = "log")

fit_ref <- glmmTMB(form, data = Rail, family = family)
cat("=== reference fit (real Rail data, matches GH #936's glmmTMB numbers) ===\n")
print(summary(fit_ref))

beta_fit <- fixef(fit_ref)$cond[["(Intercept)"]]
sd_fit <- attr(VarCorr(fit_ref)$cond$Rail, "stddev")[["(Intercept)"]]
sigma_fit <- sigma(fit_ref)
cat(sprintf("\nfitted: beta=%.4f  sd=%.4f  sigma=%.4f\n", beta_fit, sd_fit, sigma_fit))

## pretty (rounded) true parameters -- sigma=4.0 not coincidentally close
## to GH #936's own "should be ~4.0" reference point
beta_pretty <- 4.1
sd_pretty <- 0.4
sigma_pretty <- 4.0
cat(sprintf("pretty : beta=%.4f  sd=%.4f  sigma=%.4f\n", beta_pretty, sd_pretty, sigma_pretty))

args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args) >= 1) as.integer(args[[1]]) else 10
sim_seed <- 20260806

newparams <- list(beta = beta_pretty, betadisp = log(sigma_pretty), theta = log(sd_pretty))
sim_y <- simulate_new(~ 1 + (1 | Rail), nsim = B, seed = sim_seed,
                       family = family, newdata = Rail, newparams = newparams)

cat("\n=== sanity check on simulated responses ===\n")
for (b in seq_len(B)) cat(sprintf("rep %2d: n=%d range=[%.3f, %.2f]\n",
                                    b, length(sim_y[[b]]), min(sim_y[[b]]), max(sim_y[[b]])))

sim_data <- lapply(sim_y, function(y) { d <- Rail; d$travel <- y; d })

saveRDS(list(name = "rail", formula = form, family = family,
             base_data = Rail, sim_data = sim_data,
             pretty = list(beta = beta_pretty, sd = sd_pretty, sigma = sigma_pretty),
             fit_ref_summary = list(beta = beta_fit, sd = sd_fit, sigma = sigma_fit),
             B = B, seed = sim_seed),
        file.path(outdir, "rail_simdata.rds"))
cat("\nsaved to", file.path(outdir, "rail_simdata.rds"), "\n")
