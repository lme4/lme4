library("testthat")
library("lme4")

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

mk_method("NULL")
mk_method("numeric")
mk_method("array",print_dims=TRUE)
mk_method("")

test_that("simple numerics", {
    expect_equal(simulate(1~.),1)
    ## One-sided formula is not the same as an LHS that evaluates to NULL:
    expect_equal(simulate(NULL~.),NULL)
})

test_that("raw formulas", {
    expect_error(simulate(~.), '"newdata" is missing')
    expect_error(suppressWarnings(simulate(x~.)), "Error evaluating")
})

test_that("multielement classes", {
    expect_equal(simulate(diag(5)~.), diag(5))
    M <- diag(5)
    expect_equal(simulate(M~.), diag(5))
    A <- array(1,c(2,2,2))
    expect_equal(simulate(A~.), A)
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
