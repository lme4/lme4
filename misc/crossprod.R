library(Rcpp)
library(RcppEigen)
library(inline)
crossprodCpp <- readLines("crossprod.cpp") |> paste(collapse = "\n")
fcprd <- cxxfunction(signature(AA = "matrix"), crossprodCpp, plugin="RcppEigen",
                     verbose = TRUE)

fcprd(matrix(9, 3, 3))
M <- matrix(nrow = 0, ncol = 0)
storage.mode(M) <- "numeric"
fcprd(M)

## d_VtV.setZero().selfadjointView<Eigen::Upper>().rankUpdate(d_V.adjoint());

