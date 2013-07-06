## data set and formula extracted from ?prLogisticDelta example
##   (Thailand, clustered-data) in prLogistic package
load(system.file("testdata","prLogistic.RData",package="lme4"))
library(lme4)
(testLevel <- if (nzchar(s <- Sys.getenv("LME4_TEST_LEVEL"))) as.numeric(s) else 1)

if (testLevel > 2) {
    lme4_est <- glmer(rgi ~ sex + pped + (1|schoolid),
                data = dataset, family=binomial)

    lme4_results <- list(sigma= sqrt(unname(unlist(VarCorr(lme4_est)))),
                         beta = fixef(lme4_est))


    ## stored results from other pkgs
    glmmML_est <- structure(list(sigma = 1.25365353546143,
                                 beta = structure(c(-2.19478801858317,
                                 0.548884468743364, -0.623835613907385), .Names = c("(Intercept)",
                                                                       "sex", "pped"))),
                            .Names = c("sigma", "beta"))
    lme4.0_est <- structure(list(sigma = 1.25369539060849, beta = structure(c(-2.19474529099587,
                      0.548900267825802, -0.623934772981894), .Names = c("(Intercept)",
                          "sex", "pped"))), .Names = c("sigma", "beta"))

    stopifnot(all.equal(lme4_results,glmmML_est,tol=3e-3))
    stopifnot(all.equal(lme4_results,lme4.0_est,tol=3e-3))
}
