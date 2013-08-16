lme4: Mixed-effects models in R. 
====

## Features

* Efficient for large data sets, using algorithms from the 
[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
linear algebra package via the [RcppEigen](http://cran.r-project.org/web/packages/RcppEigen/index.html)
interface layer.
* Allows arbitrarily many nested and crossed random effects.
* Fits generalized linear mixed models (GLMMs) and nonlinear mixed models (NLMMs) via Laplace approximation
or adaptive Gauss-Hermite quadrature; GLMMs allow user-defined families and link functions.
* Incorporates likelihood profiling and parametric bootstrapping.

## Installation

* From CRAN (note stable version 0.999999-2 will soon be superseded by stable release 1.0.+)
* Nearly up-to-date development binaries from `lme4` r-forge repository:
```
install.packages("lme4",
   repos=c("http://lme4.r-forge.r-project.org/repos",
          getOption("repos")["CRAN"]))
```
* Development version from github:
```
library("devtools"); install_github("lme4",user="lme4")
```
(The last approach requires that you build from source, i.e. `make` and compilers must be installed on your system -- see the R FAQ for your operating system; you may also need to install dependencies manually.)
