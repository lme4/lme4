lme4: Mixed-effects models in R.
====

[![Build Status](https://travis-ci.org/lme4/lme4.svg?branch=master)](https://travis-ci.org/lme4/lme4)
[![downloads](http://cranlogs.r-pkg.org/badges/lme4)](http://cranlogs.r-pkg.org/badges/lme4)
[![cran version](http://www.r-pkg.org/badges/version/lme4)](http://cran.rstudio.com/web/packages/lme4)

## Recent/release notes

* Release 1.1-9 is on CRAN now. There are no major user-visible changes.
   * Some improvements to \code{glmer.nb}, including the ability to simulate; note that \code{glmer.nb} is still somewhat unstable, e.g. for models with additional arguments such as \code{offset} or \code{contrasts}.
   * Bug fixes to \code{predict} for models with orthogonal polynomial/spline bases, and to \code{lmList} models with offsets and weights
* Otherwise, see the [NEWS file](https://github.com/lme4/lme4/blob/master/inst/NEWS.Rd) (or  `news(Version=="1.1.9",package="lme4")`).

## Features

* Efficient for large data sets, using algorithms from the
[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
linear algebra package via the [RcppEigen](http://cran.r-project.org/package=RcppEigen)
interface layer.
* Allows arbitrarily many nested and crossed random effects.
* Fits generalized linear mixed models (GLMMs) and nonlinear mixed models (NLMMs) via Laplace approximation
or adaptive Gauss-Hermite quadrature; GLMMs allow user-defined families and link functions.
* Incorporates likelihood profiling and parametric bootstrapping.

## Installation

### On current R (>= 3.0.0)

* From CRAN (stable release 1.0.+)
* Development version from Github:
```
library("devtools"); install_github("lme4/lme4",dependencies=TRUE)
```
(This requires `devtools` >= 1.6.1, and installs the "master" (development) branch.)
This approach builds the package from source, i.e. `make` and compilers must be installed on your system -- see the R FAQ for your operating system; you may also need to install dependencies manually. Specify `build_vignettes=FALSE` if you have trouble because your system is missing some of the `LaTeX/texi2dvi` tools.
* Usually up-to-date development binaries from `lme4` r-forge repository:
```
install.packages("lme4",
   repos=c("http://lme4.r-forge.r-project.org/repos",
          getOption("repos")[["CRAN"]]))
```
(these source and binary versions are updated manually, so may be out of date; if you believe they are, please contact the maintainers).

### On old R (pre-3.0.0)

It is possible to install (but not easily to check) `lme4` at least as recently as 1.1-7.

* make sure you have *exactly* these package versions: `Rcpp` 0.10.5, `RcppEigen` 3.2.0.2
* for installation, use `--no-inst`; this is necessary in order to prevent R from getting hung up by the `knitr`-based vignettes
* running `R CMD check` is difficult, but possible if you hand-copy the contents of the `inst` directory into the installed package directory ...

### Of `lme4.0`

* `lme4.0` is a maintained version of lme4 back compatible to CRAN versions of lme4 0.99xy,
  mainly for the purpose of  *reproducible research and data analysis* which was done with 0.99xy versions of lme4.
* there have been [some](http://stackoverflow.com/questions/23662589/r-reverting-to-lme4-0-and-still-getting-inconsistent-results) [reports](http://hlplab.wordpress.com/2014/06/24/more-on-old-and-new-lme4/) of problems with `lme4.0` on R version 3.1; if someone has a specific reproducible example they'd like to donate, please contact the maintainers.
* Notably, `lme4.0` features  `getME(<mod>, "..")` which is compatible (as much as sensibly possible) with the current `lme4`'s version of `getME()`.
* You can use the `convert_old_lme4()` function to take a fitted object created with `lme4` <1.0 and convert it for use with `lme4.0`.
* It currently resides on R-forge, and you should be able to install it with
```
install.packages("lme4.0",
                 repos=c("http://lme4.r-forge.r-project.org/repos",
                         getOption("repos")[["CRAN"]]))
```
(if the binary versions are out of date or unavailable for your system, please contact the maintainers).
