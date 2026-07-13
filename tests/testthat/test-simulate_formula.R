quietly <- TRUE

## factory for making methods
mk_method <- function(class, print_dims=FALSE) {
    method <- sprintf("simulate.formula_lhs_%s",class)
    sim_generic <- function(object, nsim=1, seed=NULL, ...) {
        if (!quietly) message(sprintf("%s called",method))
        if (!quietly) cat(".Basis from attributes:\n")
        if (!quietly) print(attr(object,".Basis")) # NULL
        if (print_dims) {
            if (!quietly) print(dim(attr(object,".Basis")))
        }
        return(attr(object,".Basis"))
    }
    assign(method,sim_generic,.GlobalEnv)
    invisible(NULL)
}

## (**) these methods should (??) _mask_ package versions ...
## works in source(), not in devtools::test() ...
mk_method("NULL")
mk_method("numeric")
mk_method("array",print_dims=TRUE)
mk_method("")

test_that("simple numerics", {
    ## expect_equal(simulate(1~.),1)  ## FIXME re-enable if we resolve (**) above
    ## One-sided formula is not the same as an LHS that evaluates to NULL:
    expect_equal(simulate(NULL~.),NULL)
})

test_that("raw formulas", {
    expect_error(suppressWarnings(simulate(x~.)), "Error evaluating")
})

simulate.formula_lhs_character <- function(object, nsim=1, seed=NULL, ...) {
    if (!quietly) message("simulate.formula_lhs_character() called.")
    if (!quietly) print(ls(all.names=TRUE))
    NextMethod() # Calls simulate.formula(), resulting in an infinite recursion.
}

test_that("prevent recursion", {
     expect_error(simulate("a"~.), "No applicable method")
})

dd <- expand.grid(A=factor(1:3),B=factor(1:10),rep=1:10)
test_that("two-sided formula warning", {
    expect_error(suppressMessages(simulate(.~1 + (A|B),
                                             newdata=dd,
                                             newparams=list(beta=1,theta=rep(1,6),
                                                            sigma=1),
                                             family=gaussian,
                                             seed=101))[[1]],
                   "object '.' not found")
})

## cleanup: mk_method() assigns into .GlobalEnv, so remove explicitly from there
suppressWarnings(try(rm(list = c("simulate.formula_lhs_", "simulate.formula_lhs_numeric",
                                 "simulate.formula_lhs_array"),
                        envir = .GlobalEnv), silent = TRUE))

## a one-sided formula with neither 'family' nor 'newparams' used to fall
## through to the generic (and unhelpful) "No applicable method for LHS of
## type 'NULL'" error; it should instead report the missing arguments.
test_that("one-sided formula without family/newparams gives an informative error", {
    expect_error(simulate(~1 + (1|Days), newdata=sleepstudy),
                 "must specify all of.*formula.*newdata.*newparams")
})

## GH#948: simulate with re.form=NULL should work when a predictor is a 1-col
## matrix (e.g. from scale()) used inside poly() -- poly() with pre-specified
## coefs rejects matrix input, so mkNewReTrms must drop such columns first
test_that("simulate re.form=NULL works with poly(scale(x)) predictor (GH#948)", {
  ss_sc <- transform(sleepstudy, Days = scale(Days))
  fit1  <- lmer(Reaction ~ 1 + poly(Days, 2) + (1|Subject), data = sleepstudy)
  fit_sc <- update(fit1, data = ss_sc)
  s1 <- simulate(fit_sc, seed = 101, re.form = NULL)
  s2 <- simulate(fit_sc, seed = 101, re.form = NULL)
  expect_equal(s1, s2)
})

