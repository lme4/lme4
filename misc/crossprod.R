library(Rcpp)
library(RcppEigen)
library(inline)
crossprodCpp <- readLines("crossprod.cpp") |> paste(collapse = "\n")
fcprd <- cxxfunction(signature(RV = "matrix", RVtV = "matrix"), crossprodCpp, plugin="RcppEigen",
                     verbose = TRUE)

testfun <- function(m, n) {
    M <- matrix(0.0, m, n)
    result <- matrix(NA_real_, n, n)
    fcprd(M, result)
    result
}

testfun(3,3)
testfun(0,0)
testfun(30,0)
## note that ASAN error is only printed the first time it occurs for any particular executable ...
## /usr/local/lib/R/site-library/RcppEigen/include/Eigen/src/Core/util/BlasUtil.h:193:18: runtime error: reference binding to misaligned address 0x000000000001 for type 'const double', which requires 8 byte alignment
## 0x000000000001: note: pointer points here
## <memory cannot be printed>

testfun(30,0)

library(lme4)
fm01ML <- lmer(Yield ~ 1|Batch, Dyestuff, REML = FALSE)
tpr  <- profile(fm01ML, optimizer="Nelder_Mead", which="beta_")
