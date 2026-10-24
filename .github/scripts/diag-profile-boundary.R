## Diagnostic for the "profiling correlations near or on the boundary" test
## (tests/testthat/test-methods.R), which fails on CRAN r-release-macos-x86_64
## because profile() does not warn "non-monotonic profile".
## Run with: Rscript .github/scripts/diag-profile-boundary.R
library(lme4)
options(digits = 17, width = 150)

cat("==== platform ====\n")
print(sessionInfo())
print(extSoftVersion())
print(La_library())
print(La_version())
cat("CC:  ", system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CC"), stdout = TRUE), "\n")
cat("CXX: ", system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CXX"), stdout = TRUE), "\n")
cat("lme4 version:", format(packageVersion("lme4")), "\n")

cat("\n==== fit (same code as the test) ====\n")
ctrl <- lmerControl(check.conv.singular = "ignore")
set.seed(101)
d2 <- data.frame(g = factor(rep(1:6, each = 5)), x = rnorm(30))
d2$y <- rnorm(30) + rep(rnorm(6), each = 5) * (1 + d2$x)
f2 <- lmer(y ~ x + (x | g), d2, control = ctrl)
print(summary(f2))
th <- getME(f2, "theta")
cat("theta (%a):", sprintf("%a", th), "\n")
cat("theta:     ", format(th, digits = 17), "\n")
cat("REMLcrit:  ", format(REMLcrit(f2), digits = 17), "\n")
cat("on boundary (th[3] == 0 && th[2] != 0):", th[3] == 0 && th[2] != 0, "\n")

cat("\n==== profile(f2, which = 3) ====\n")
warns <- character(0)
p2 <- withCallingHandlers(
    profile(f2, which = 3),
    warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
cat("warnings issued:", length(warns), "\n")
if (length(warns)) print(warns)

pp <- as.data.frame(p2)
pp <- pp[pp$.par == ".sig03", ]
pp <- pp[order(pp$.sig03), ]
cat("\nprofile slice for .sig03 (sorted):\n")
print(pp[, c(".zeta", ".sig01", ".sig02", ".sig03", ".sigma")], digits = 17)

d <- na.omit(pp[, c(".sig03", ".zeta")])
dz <- diff(d$.zeta)
cat("\ndiff(.zeta):\n")
print(dz, digits = 17)
cat("n points:", nrow(d), " duplicated .sig03:", anyDuplicated(d$.sig03),
    "\n#dz > 0:", sum(dz > 0), " #dz < 0:", sum(dz < 0), " #dz == 0:", sum(dz == 0),
    "\nstrictly monotone (=> monoHermite fallback, no warning):",
    all(dz > 0) || all(dz < 0), "\n")

cat("\nspline classes stored in profile object:\n")
cat("forward: ", class(attr(p2, "forward")[[".sig03"]]), "\n")
cat("backward:", class(attr(p2, "backward")[[".sig03"]]), "\n")

cat("\n==== confint ====\n")
print(suppressWarnings(confint(p2)), digits = 17)

saveRDS(list(theta = th, profile = p2, warnings = warns),
        file = "diag-profile-boundary.rds")
