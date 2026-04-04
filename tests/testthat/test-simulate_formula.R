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

## cleanup

## I can't figure out what environments these things actually live in so I'm going to
## give up and try() to remove them ...

## rmx <- function(s) if (exists(s, parent.frame())) rm(list=s, envir=parent.frame())
## rmx("simulate.formula_lhs_character")
## rmx("simulate.formula_lhs_")
## rmx("simulate.formula_lhs_numeric")

suppressWarnings(try(rm(list = c("simulate.formula_lhs_", "simulate.formula_lhs_numeric")),silent=TRUE))

## in general, we shouldn't have it such that weights and offsets leak into the
## global environment.

test_that("ensuring weights and offsets don't leak into global environment", {
  # weights have not yet been assigned
  form <- ~x + (1|f)
  # First simulate - with weights of 10
  d1 <- data.frame(x=rnorm(60), f=factor(rep(1:6,each=10)), w=rep(10,60))
  d1$y1 <- simulate(form, family=binomial, weights=d1$w, newdata=d1,
                   newparams=list(theta=0.01, beta=c(1,1)))[[1]]
  # ensures it was not leaked in the local environment!
  expect_false(exists("weights", envir = environment(form), 
                      inherits = FALSE))
  expect_false(exists("offset", envir = environment(form), 
                      inherits = FALSE))
  # weights and offsets have been assigned
  weights <- offset <- 1:10  
  offset <- 1:10
  d1$y2 <- simulate(form, family=binomial, weights=d1$w, newdata=d1,
                   newparams=list(theta=0.01, beta=c(1,1)))[[1]]
  expect_equal(weights, 1:10)
  expect_equal(offset, 1:10)
})
