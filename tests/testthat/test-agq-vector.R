## Adaptive Gauss-Hermite quadrature (nAGQ > 1) for a single vector-valued
## random-effects term, e.g. (1 + x | g): glmerAGQ() in src/external.cpp,
## AGQrule() in R/GHrule.R.

test_that("AGQrule() reduces to GHrule() and integrates low-order moments", {
    expect_identical(lme4:::AGQrule(7), GHrule(7))
    sp1 <- lme4:::AGQrule(5, 1, "sparse")
    expect_equal(unname(sp1[order(sp1[, 1]), ]), unname(GHrule(5)))
    for (type in c("product", "sparse")) {
        r <- lme4:::AGQrule(5, 3, type)
        expect_equal(ncol(r), 5L)
        w <- r[, "w"]; z <- r[, 1:3]
        expect_equal(sum(w), 1)
        expect_equal(unname(colSums(w * z^2)), rep(1, 3))
        expect_equal(unname(colSums(w * z^4)), rep(3, 3))
        expect_equal(sum(w * z[, 1]^2 * z[, 2]^2), 1)
        expect_equal(r[, "ldnorm"], rowSums(dnorm(z, log = TRUE)))
    }
    expect_equal(nrow(lme4:::AGQrule(5, 2)), 25L)
    expect_lt(nrow(lme4:::AGQrule(5, 3, prune = 1e-4)), 125L)
})
