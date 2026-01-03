testLevel <- if (nzchar(s <- Sys.getenv("LME4_TEST_LEVEL"))) as.numeric(s) else 1
if (testLevel>1) {

  L <- load(system.file("testdata", "lme-tst-fits.rda",
                        package="lme4", mustWork=TRUE))
  
  had_pars <- exists("pars", envir = globalenv(), inherits = FALSE)
  had_ctrl <- exists("ctrl", envir = globalenv(), inherits = FALSE)
  
  gm_all <- allFit(fit_cbpp_1, verbose=FALSE)
  gm_all_nostart <- allFit(fit_cbpp_1, verbose=FALSE, start_from_mle = FALSE)

  test_that("pars and ctrl did not leak into the environment", {
    expect_equal(exists("pars", envir = globalenv(), inherits = FALSE), 
                 had_pars)
    expect_equal(exists("ctrl", envir = globalenv(), inherits = FALSE), 
                 had_ctrl)
  })
  
  summary(gm_all)$times[,"elapsed"]
  summary(gm_all_nostart)$times[,"elapsed"]

  sapply(gm_all, function(x) x@optinfo$feval)
  sapply(gm_all_nostart, function(x) x@optinfo$feval)

  ## library(microbenchmark)
  ##  mb1 <- microbenchmark(
  ##   start = allFit(fit_cbpp_1, verbose=FALSE),
  ##   nostart = allFit(fit_cbpp_1, verbose=FALSE, start_from_mle = FALSE)
  ##  )

    test_that("allFit print/summary is fine", {
        expect_is(gm_all, "allFit")
        expect_is(summary(gm_all), "summary.allFit")
    })

    test_that("nloptwrap switches optimizer correctly", {
        expect_equal(attr(gm_all[["nloptwrap.NLOPT_LN_BOBYQA"]],"optCtrl"),
                     list(maxeval = 1e5, algorithm = "NLOPT_LN_BOBYQA"))
        expect_equal(attr(gm_all[["nloptwrap.NLOPT_LN_NELDERMEAD"]],"optCtrl"),
                     list(maxeval = 1e5, algorithm = "NLOPT_LN_NELDERMEAD"))

    })

    test_that("lmerControl() arg works too", {
        fm0 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
        fm  <- update(fm0,
                      control = lmerControl(optCtrl = list(xtol_rel = 1e-8,
                                                           ftol_rel = 1e-8),
                                            calc.derivs=FALSE))
        afm0 <- allFit(fm0,verbose=FALSE)
        afm  <- allFit(fm,verbose=FALSE) # used to fail
        drop_ <- function(x) {
            x[setdiff(names(x), c("times","feval"))]
        }
        ## should be approximately the same
        expect_equal(drop_(summary(afm0)),
                     drop_(summary(afm)), tolerance = 1e-2)
        ## should NOT be the same!
        expect_false(isTRUE(all.equal(drop_(summary(afm0)),
                                      drop_(summary(afm)), tolerance=1e-10)))

    })

    test_that("glmerControl() arg + optimizer", {
        ## GH #523?
        fit_cbpp_1u <- update(fit_cbpp_1,
                              control=glmerControl(optimizer="nloptwrap",
                                                   optCtrl=list(xtol_abs=1e-10, ftol_abs=1e-10)))
        af2 <- allFit(fit_cbpp_1u, verbose=FALSE)
        expect_equal(class(af2),"allFit")
    })

    test_that("i in model call is OK", {
        ## GH #538
        ## ugh, testthat scoping is insane ...
        ## if d and i are
        ## assigned normally with <- outside expect_true(), test fails
        ## BUT global assignment of 'd' breaks downstream tests in
        ##  'data= argument and formula evaluation' (test-formulaEval.R)
        ## ddd breaks similar test in 'fitting lmer models' (test-lmer.R)
        ##  (where 'd' is supposed to be nonexistent)
        ## if we do global assignment with <<-
        ##   can't figure out how to remove d (or ddd) after it's created to leave
        ##   the environment clean ...
        ## tried to encapsulate all the necessary assignments
        ## within expect_true({ ... }) but that fails in other ways

        nr <- nrow(sleepstudy)
        ..dd <<- list(sleepstudy[1:nr,], sleepstudy[-(1:nr)])
        i <<- 1
        fm0 <- lmer(Reaction ~ Days + (1 | Subject), data=..dd[[i]])
        aa <- allFit(fm0, verbose=FALSE)
        expect_true(
            all(summary(aa)$which.OK)
        )
    })

    test_that("allFit/update scoping", {
        ## GH #601
        fit_func <- function(dataset) {
            gm1 <- glmer(
                cbind(incidence, size - incidence) ~ period + (1 | herd),
                data = dataset, family = binomial
            )
            gm1@call$data <- dataset
            allFit(gm1, catch.errs=FALSE)
        }

        cc <- capture.output(ff <- fit_func(cbpp))
        expect_true(all(summary(ff)$which.OK))
    })

    test_that("maxfun works", {
        gm_it10 <- suppressWarnings(allFit(fit_cbpp_1, verbose=FALSE, maxfun = 10))
        v <- vapply(gm_it10, function(x) as.integer(x@optinfo$feval), FUN.VALUE=1L)
        ## function values are sometimes off a bit (due to initialization or Hessian calculation)
        ##  but close enough ...
        expect_true(all(is.na(v) | v < 12))
    })

}  ## testLevel
