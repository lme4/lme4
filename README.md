lme4: Mixed-effects models in R. 
====

[![Build Status](https://travis-ci.org/lme4/lme4.svg?branch=master)](https://travis-ci.org/lme4/lme4)

## Recent/release notes

* Recent versions of `lme4` (e.g. 1.1-6) give false convergence warnings. There is a [summary post on r-sig-mixed-models](http://thread.gmane.org/gmane.comp.lang.r.lme4.devel/11893).  
   * If you get warnings about `max|grad|` but the  model passes this test:
```
dd <- fit@optinfo$derivs
with(dd,max(abs(solve(Hessian,gradient)))<2e-3)
```
then you are seeing a false-positive warning, and the problem will disappear in future versions (1.1-7 and up).  
   * For other warnings (e.g. about the Hessian being singular or having negative eigenvalues), you can try centering and/or scaling continuous predictor variables. 
   * You can also try (for `glmer` fits)  `control=glmerControl(optimizer="bobyqa")`, or use [this code](https://github.com/lme4/lme4/blob/master/misc/issues/allFit.R) to try your problem with a range of optimizers, to see if any of them work better.
* If your convergence warnings persist, the `lme4` maintainers would be happy to hear from you.


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

* From CRAN (stable release 1.0.+)
* Development version from Github:
```
library("devtools"); install_github("lme4",user="lme4")
```
(These commands install the "master" (development) branch; if you
want the release branch from Github add `ref="release"` to the
`install_github()` call.
The `install_github()` approach requires that you build from source, i.e. `make` and compilers must be installed on your system -- see the R FAQ for your operating system; you may also need to install dependencies manually. )
* Nearly up-to-date development binaries from `lme4` r-forge repository:
```
install.packages("lme4",
   repos=c("http://lme4.r-forge.r-project.org/repos",
          getOption("repos")[["CRAN"]]))
```

## Installation of `lme4.0`

* `lme4.0` is a maintained version of lme4 back compatible to CRAN versions of lme4 0.99xy,
  mainly for the purpose of  *reproducible research and data analysis* which was done with 0.99xy versions of lme4.
* Notably, `lme4.0` features  `getME(<mod>, "..")` which is compatible (as much as sensibly possible) with the current `lme4`'s version of `getME()`.
* You can use the `convert_old_lme4()` function to take a fitted object created with `lme4` <1.0 and convert it for use with `lme4.0`.
* It currently resides on R-forge, and you can install it with

```
install.packages("lme4.0", 
                 repos=c("http://lme4.r-forge.r-project.org/repos",
                         getOption("repos")[["CRAN"]]))
```
