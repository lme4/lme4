## Compare glmer() with glmmTMB and GLMMadaptive for binomial and Poisson
## GLMMs with non-integer responses (GH #180, #682, #1011; PR #1012).
##
## Laplace fits (glmer nAGQ = 1 vs glmmTMB) and adaptive GHQ fits
## (glmer nAGQ = 11 vs GLMMadaptive, whose default is 11 points) should
## agree, as should their log-likelihoods: all three use the continuous
## extension of the binomial/Poisson density (lgamma() for factorials).
##
## GLMMadaptive's Poisson log-density is already continuous, but its
## binomial one calls dbinom() (-Inf for non-integer counts), so we give
## it a binomial() family with a continuous log_dens. GLMMadaptive's
## weights are per cluster, so binomial data are passed as cbind(k, n - k).
##
## GLMMadaptive's default convergence settings stop short of glmer's AGQ
## optimum (differences ~1e-2 even for integer data), so we tighten them.
##
## glmmTMB's log-likelihood omits the binomial coefficient when the
## binomial size is 1 (fractional responses with weights 1), so for those
## data its logLik is lower than glmer's/GLMMadaptive's by
## sum(lchoose(1, y)) (continuous lchoose); parameters still agree.
##
## Data: integer examples (as a baseline) and the same examples with
## y + runif(n, -0.5, 0.5) (clipped to the support), plus the fractional
## binomial examples from GH #1011 / PR #1012.
##
## Usage: Rscript misc/noninteger_compare.R  [with R_LIBS pointing at the
## lme4 build to test]

suppressPackageStartupMessages({
    library(lme4)
    library(glmmTMB)
    library(GLMMadaptive)
})
cat("lme4", format(packageVersion("lme4")), "from", find.package("lme4"),
    "\nglmmTMB", format(packageVersion("glmmTMB")),
    " GLMMadaptive", format(packageVersion("GLMMadaptive")), "\n\n")

## continuous binomial log-density for GLMMadaptive (y is cbind(k, n - k))
lchoose_ext <- function(n, k) lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
binomial_cont <- function() {
    fam <- binomial()
    fam$log_dens <- function(y, eta, mu_fun, phis, eta_zi) {
        mu_y <- mu_fun(eta)
        k <- y[, 1L]; n <- y[, 1L] + y[, 2L]
        out <- lchoose_ext(n, k) + k * log(mu_y) + (n - k) * log1p(-mu_y)
        attr(out, "mu_y") <- mu_y
        out
    }
    fam
}

## tighter than GLMMadaptive's defaults (tol1 = 1e-4, tol2 = 1e-5,
## iter_qN_outer = 15, update_GH_every = 10, numeric_deriv = "fd")
ga_control <- list(tol1 = 1e-8, tol2 = 1e-8, iter_qN_outer = 50,
                   update_GH_every = 1, numeric_deriv = "cd")

quiet <- function(expr) suppressWarnings(suppressMessages(expr))

## summary: fixed effects, RE standard deviation, log-likelihood
summ_glmer <- function(f)
    c(fixef(f), sd = unname(getME(f, "theta")), logLik = as.numeric(logLik(f)))
summ_tmb <- function(f)
    c(fixef(f)$cond, sd = sqrt(unname(VarCorr(f)$cond[[1]][1, 1])),
      logLik = as.numeric(logLik(f)))
summ_ga <- function(f)
    c(fixef(f), sd = sqrt(unname(f$D[1, 1])), logLik = as.numeric(logLik(f)))

## fit one example with all methods; 'resp' is the response expression as a
## string (cbind(k, n - k) for binomial), 'rhs' the fixed-effect RHS
compare <- function(label, data, resp, rhs, grp, family) {
    ffix <- reformulate(rhs, response = str2lang(resp))
    fmix <- reformulate(c(rhs, sprintf("(1 | %s)", grp)), response = str2lang(resp))
    ga_fam <- if (family == "binomial") binomial_cont() else poisson()
    fam <- get(family)
    res <- rbind(
        glmer_L   = summ_glmer(quiet(glmer(fmix, data, family = fam))),
        glmmTMB_L = summ_tmb(quiet(glmmTMB(fmix, data, family = fam))),
        glmer_Q11 = summ_glmer(quiet(glmer(fmix, data, family = fam, nAGQ = 11))),
        GLMMad_Q11 = summ_ga(quiet(mixed_model(ffix,
                                               random = reformulate(sprintf("1 | %s", grp)),
                                               data = data, family = ga_fam,
                                               control = ga_control))))
    cat("====", label, "====\n")
    print(signif(res, 6))
    pars <- setdiff(colnames(res), "logLik")
    dpar <- function(a, b) max(abs(res[a, pars] - res[b, pars]))
    dll <- function(a, b) unname(res[a, "logLik"] - res[b, "logLik"])
    r <- c(par_L = dpar("glmer_L", "glmmTMB_L"), ll_L = dll("glmer_L", "glmmTMB_L"),
           par_Q = dpar("glmer_Q11", "GLMMad_Q11"), ll_Q = dll("glmer_Q11", "GLMMad_Q11"))
    cat(sprintf(paste("Laplace (glmer - glmmTMB): max |par diff| %.2e, logLik diff %.2e",
                      "AGQ (glmer - GLMMadaptive): max |par diff| %.2e, logLik diff %.2e\n\n",
                      sep = "\n"), r[1], r[2], r[3], r[4]))
    invisible(r)
}

jitter_counts <- function(y, lower = 0, upper = Inf, seed = 101) {
    set.seed(seed)
    pmin(pmax(y + runif(length(y), -0.5, 0.5), lower), upper)
}

## ---- data
## Poisson: grouseticks (single grouping factor for GLMMadaptive)
gt <- transform(grouseticks, cHEIGHT = HEIGHT - mean(HEIGHT))
gt$TICKS_j <- jitter_counts(gt$TICKS)
## Poisson: simulated example from GH #180
set.seed(101)
sim <- expand.grid(fac1 = LETTERS[1:2], grp = factor(1:10), obs = 1:10)
sim$y <- simulate(~ fac1 + (1 | grp), family = poisson, newdata = sim,
                  newparams = list(beta = c(2, 1), theta = 1))[[1]]
sim$y_j <- jitter_counts(sim$y)
## binomial: cbpp
cb <- cbpp
cb$inc_j <- jitter_counts(cb$incidence, upper = cb$size)
## binomial: fractional responses with weights 1 (GH #1011; theta = 0)
d1011 <- data.frame(
    y = c(0.2565, 0.3552, 0.344, 0.3393, 0.2216, 0.2628, 0.479, 0.2689, 0.2791,
          0.1885, 0.3799, 0.4617, 0.2209, 0.4506, 0.4922, 0.5739, 0.8292, 0.6354,
          0.776, 0.048, 0.5312, 0.0616, 0.1741, 0.9788, 0.917, 0.3709, 0.5611,
          0.111, 0.0479, 0.3284, 0.1255, 0.1129, 0.1322, 0.2209, 0.1343, 0.1795,
          0.0932, 0.1599, 0.3264),
    x = c(0.7458, 0.2496, 0.3323, 0.116, 2.0572, 1.0248, -0.6979, 1.0392, 1.1827,
          2.0272, 0.3178, -0.0954, -2.0399, -1.7392, -1.6437, -0.8398, 2.3697,
          0.0016, 0.1374, -0.769, 0.0934, -0.7415, -0.4593, 1.0565, 1.3049,
          0.1363, 0.9627, 0.0119, -0.3273, 0.6624, -2.0192, -1.7293, -0.9185,
          0.9509, -1.1285, -0.2256, 1.1649, -0.3001, -0.6718),
    f = factor(rep(c("a", "b", "c"), c(13, 13, 13))),
    g = factor(paste0("g", c(rep(1:6, each = 6), 1, 1, 1))))
## binomial: fractional responses with theta > 0 (PR #1012 discussion)
set.seed(11)
g <- factor(rep(1:10, each = 30)); x <- rnorm(300)
dth <- data.frame(y = plogis(-0.3 + 0.8 * x + rnorm(10, sd = 1.5)[g]), x = x, g = g)

## ---- comparisons
labels <- c("grouseticks, integer", "grouseticks, jittered",
            "GH #180 sim, integer", "GH #180 sim, jittered",
            "cbpp, integer", "cbpp, jittered",
            "GH #1011 fractional", "PR #1012 fractional, theta > 0")
out <- rbind(
    compare("grouseticks, integer", gt, "TICKS", c("YEAR", "cHEIGHT"), "BROOD", "poisson"),
    compare("grouseticks, jittered", gt, "TICKS_j", c("YEAR", "cHEIGHT"), "BROOD", "poisson"),
    compare("GH #180 sim, integer", sim, "y", "fac1", "grp", "poisson"),
    compare("GH #180 sim, jittered", sim, "y_j", "fac1", "grp", "poisson"),
    compare("cbpp, integer", cb, "cbind(incidence, size - incidence)", "period", "herd", "binomial"),
    compare("cbpp, jittered", cb, "cbind(inc_j, size - inc_j)", "period", "herd", "binomial"),
    compare("GH #1011 fractional", d1011, "cbind(y, 1 - y)", c("x", "f"), "g", "binomial"),
    compare("PR #1012 fractional, theta > 0", dth, "cbind(y, 1 - y)", "x", "g", "binomial"))
rownames(out) <- labels
cat("==== summary: max |diff| over fixed effects and RE sd; logLik differences ====\n")
print(signif(as.data.frame(out), 3))
cat("\nsize-1 binomial constant sum(lchoose(1, y)) (expected glmer - glmmTMB logLik gap):\n")
print(c(`GH #1011` = sum(lchoose_ext(1, d1011$y)),
        `PR #1012` = sum(lchoose_ext(1, dth$y))))
