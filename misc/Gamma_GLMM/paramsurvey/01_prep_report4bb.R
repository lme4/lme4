## Prep step for the "Report4BB" example: the real dataset/model that
## kicked off this whole investigation (README #1), reconstructed from
## https://github.com/TiagoAMarques/Report4BB (cloned to
## ../report4bb_repo; see testing_lme4/Testing_lme4.R for the original
## data-prep + crglmer model code this mirrors). crate ~ (1|location) +
## (1|fyear), Gamma(log), n=103 -- TWO independent (crossed, not nested)
## single-term random-intercept grouping factors, with few levels each (7
## locations, 10 years): a genuinely different regime from both epil2
## examples (no correlation parameter at all, but two separate small-group
## RE terms instead of one).

suppressMessages({ library(tidyverse); library(glmmTMB) })
source("/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/param_survey/toolkit.R")

outdir <- "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/param_survey"
repo <- "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb_repo"

## --- reconstruct the `tags` dataset (verbatim from testing_lme4/Testing_lme4.R) ---
load(file.path(repo, "data_4_article_clickrates_deep_dive.rda"))
DDCs <- ddata1[ddata1$sonar != "sonar", ]
tags <- DDCs %>%
  group_by(tag) %>%
  summarise(location = unique(location), year = unique(year), sex = unique(sex),
            duration = sum(durations, na.rm = TRUE), nclicks = sum(nclick, na.rm = TRUE),
            crate = sum(nclick, na.rm = TRUE) / sum(durations, na.rm = TRUE),
            ddc = max(absdives + 1, na.rm = TRUE))
tags <- tags[tags$location != "Norway Andenes", ]
tags$year[tags$location == "Norway" & tags$year == 2009] <- 2010
tags$year[tags$location == "DOMINICA" & tags$year == 2017] <- 2016
tags$year[tags$location == "Mediterranean" & tags$year == 2001] <- 2003
tags$fyear <- as.factor(tags$year)
tags <- as.data.frame(tags)
tags$location <- factor(tags$location)
cat("nrow:", nrow(tags), " n location:", nlevels(tags$location), " n fyear:", nlevels(tags$fyear), "\n")

form <- crate ~ (1 | location) + (1 | fyear)
family <- Gamma(link = "log")

## non-singular without any prior (verified interactively) -- unlike the
## epil2 (complex) example, no regularizing prior is needed here.
fit_ref <- glmmTMB(form, data = tags, family = family)
cat("=== reference fit ===\n")
print(summary(fit_ref))
stopifnot(!performance::check_singularity(fit_ref))

beta_fit <- fixef(fit_ref)$cond
cat("\nfitted beta:\n"); print(beta_fit)
beta_pretty <- c(`(Intercept)` = -0.1)
stopifnot(identical(names(beta_pretty), names(beta_fit)))

vc_loc <- VarCorr(fit_ref)$cond$location
vc_yr <- VarCorr(fit_ref)$cond$fyear
cat("\nfitted sd(location):", attr(vc_loc, "stddev"), " sd(fyear):", attr(vc_yr, "stddev"), "\n")
sd_pretty <- c(location = 0.10, fyear = 0.06)

phi_fit <- sigma(fit_ref)^2
cat("fitted dispersion (phi):", phi_fit, "\n")
phi_pretty <- 0.10

cat("\n=== pretty (rounded) true parameters ===\n")
print(beta_pretty); cat("sd:", sd_pretty, " phi:", phi_pretty, "\n")

## build TMB obj template and overwrite beta/betadisp/theta with pretty
## values. theta here is TWO independent single-term grouping factors
## (location, fyear) -- both plain log(sd), no correlation transform
## needed (unlike epil2 complex's single 2-term block).
pre <- glmmTMB(form, data = tags, family = family, doFit = FALSE)
obj <- glmmTMB:::fitTMB(pre, doOptim = FALSE)
pf <- obj$env$last.par.best
stopifnot(sum(names(pf) == "beta") == 1, sum(names(pf) == "betadisp") == 1,
          sum(names(pf) == "theta") == 2)

pf[names(pf) == "beta"] <- unname(beta_pretty)
pf[names(pf) == "betadisp"] <- -log(phi_pretty)
pf[names(pf) == "theta"] <- log(sd_pretty)

args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args) >= 1) as.integer(args[[1]]) else 10
sim_seed <- 20260802
set.seed(sim_seed)
sim_y <- lapply(seq_len(B), function(b) obj$simulate(par = pf)$yobs)

cat("\n=== sanity check on simulated responses ===\n")
for (b in seq_len(B)) cat(sprintf("rep %2d: n=%d range=[%.3f, %.2f]\n",
                                    b, length(sim_y[[b]]), min(sim_y[[b]]), max(sim_y[[b]])))

sim_data <- lapply(sim_y, function(y) { d <- tags; d$crate <- y; d })

saveRDS(list(name = "report4bb", formula = form, family = family,
             base_data = tags, sim_data = sim_data,
             pretty = list(beta = beta_pretty, sd = sd_pretty, corr = NA_real_, phi = phi_pretty),
             fit_ref_summary = list(beta = beta_fit, sd = c(attr(vc_loc, "stddev"), attr(vc_yr, "stddev")),
                                     phi = phi_fit),
             B = B, seed = sim_seed),
        file.path(outdir, "report4bb_simdata.rds"))
cat("\nsaved to", file.path(outdir, "report4bb_simdata.rds"), "\n")
