library(Rcpp)
library(RcppEigen)
library(inline)
Sys.setenv("UBSAN_OPTIONS"='print_stacktrace=1')
crossprodCpp <- readLines("crossprod.cpp") |> paste(collapse = "\n")
fcprd <- cxxfunction(signature(RV = "matrix", RVtV = "matrix"), crossprodCpp, plugin="RcppEigen",
                     verbose = TRUE)


result <- matrix(NA_real_, 3, 3)
fcprd(matrix(9, 3, 3), result)
result
M <- matrix(nrow = 0, ncol = 0)
storage.mode(M) <- "numeric"
result <- matrix(0.0, 0, 0)
fcprd(M, result)

library(lme4)
fm01ML <- lmer(Yield ~ 1|Batch, Dyestuff, REML = FALSE)
tpr  <- profile(fm01ML, optimizer="Nelder_Mead", which="beta_")
