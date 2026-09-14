library("testthat")
library("lme4")

test_that("vcov() falls back to RX when the Hessian-based SEs are implausible", {
    gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                 data = cbpp, family = binomial)
    expect_true(!is.null(gm1@optinfo$derivs$Hessian))
    ## an intact fit: Hessian-based SEs are within a factor of 2 of the RX ones
    ## and are used without a warning
    expect_silent(v1 <- vcov(gm1))
    expect_equal(as.matrix(v1), as.matrix(vcov(gm1, use.hessian = TRUE)))
    ## a corrupted Hessian (as observed on fits that also raise a max|grad|
    ## warning): SEs 100 times too small -> warning and RX-based matrix
    gm2 <- gm1
    gm2@optinfo$derivs$Hessian <- gm1@optinfo$derivs$Hessian * 1e4
    expect_warning(v2 <- vcov(gm2), "differ from the\nRX-based ones by a factor of up to")
    expect_equal(as.matrix(v2), as.matrix(suppressWarnings(vcov(gm1, use.hessian = FALSE))))
    ## the check can be switched off
    op <- options(lme4.vcov.hess.se.ratio = Inf); on.exit(options(op))
    expect_silent(v3 <- vcov(gm2))
    expect_equal(as.matrix(v3), as.matrix(v1) / 1e4)
})
