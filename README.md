lme4: Mixed-effects models in R.
====

<!-- badges: start -->
[![R-CMD-check](https://github.com/lme4/lme4/workflows/R-CMD-check/badge.svg)](https://github.com/lme4/lme4/actions)
[![cran version](http://www.r-pkg.org/badges/version/lme4)](https://cran.r-project.org/package=lme4)
[![downloads](http://cranlogs.r-pkg.org/badges/lme4)](http://cranlogs.r-pkg.org/badges/lme4)
[![total downloads](http://cranlogs.r-pkg.org/badges/grand-total/lme4)](http://cranlogs.r-pkg.org/badges/grand-total/lme4)
<!-- badges: start -->

## Recent/release notes

* See the [NEWS file](https://github.com/lme4/lme4/blob/master/inst/NEWS.Rd)

## Where to get help

- [r-sig-mixed-models@r-project.org](https://stat.ethz.ch/mailman/listinfo/r-sig-mixed-models) for questions about `lme4` usage and more general mixed model questions; please read the info page, and subscribe, before posting ... (note that the mailing list does not support images or large/non-text attachments)
- https://github.com/lme4/lme4/issues for bug, infelicity, and wishlist reporting
- The [lme4 tag on StackOverflow](https://stackoverflow.com/questions/tagged/lme4) for programming-related or the [lme4-nlme tag on CrossValidated](https://stats.stackexchange.com/questions/tagged/lme4-nlme) for statistics-related questions
- maintainer e-mail only for urgent/private communications

## Support

If you choose to support `lme4` development financially, you can contribute to a fund at McMaster University (home institution of one of the developers) [here](https://secureca.imodules.com/s/1439/17/giving/form.aspx?sid=1439&gid=1&pgid=770&cid=1618&dids=2413&bledit=1&appealcode=18C9). The form will say that you are donating to the "Global Coding Fund"; this fund is available for use by the developers, under McMaster's research spending rules. We plan to use the funds, as available, to pay students to do maintenance and development work. There is no way to earmark funds or set up a bounty to direct funding toward particular features, but you can e-mail the maintainers and suggest priorities for your donation.

## Features

* Efficient for large data sets, using algorithms from the
[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
linear algebra package via the [RcppEigen](https://cran.r-project.org/package=RcppEigen)
interface layer.
* Allows arbitrarily many nested and crossed random effects.
* Fits generalized linear mixed models (GLMMs) and nonlinear mixed models (NLMMs) via Laplace approximation
or adaptive Gauss-Hermite quadrature; GLMMs allow user-defined families and link functions.
* Incorporates likelihood profiling and parametric bootstrapping.

## Installation

### On current R (>= 3.6.0)

* From CRAN (stable release 1.+)
* Development version from Github:
```r
library("devtools"); install_github("lme4/lme4",dependencies=TRUE)
```
(This requires `devtools` >= 1.6.1, and installs the "master" (development) branch.)
This approach builds the package from source, i.e. `make` and compilers must be installed on your system -- see the R FAQ for your operating system; you may also need to install dependencies manually. Specify `build_vignettes=FALSE` if you have trouble because your system is missing some of the `LaTeX/texi2dvi` tools.

* Development binaries from r-universe:
```r
install.packages('lme4', repos = c('https://lme4.r-universe.dev', getOption("repos")[["CRAN"]]))
```

## Development notes

`lme4` is developed in a mixture of

* traditional R package building tools, as documented in [Writing R Extensions](cran.r-project.org/doc/manuals/r-devel/R-exts.html#Documenting-functions)
   * NEWS in `inst/NEWS.Rd` (not a top-level `NEWS.md` file)
   * documentation as `.Rd` files (*not* `roxygen2`, although some functions have internal roxygen-style documentation [not used])
   * 'classic' tests in the `tests/` directory
   * some Sweave (`knitr`)/`Rnw`-format vignette, especially `vignettes/lmer.Rnw`
* 'tidyverse'-style tools, as documented in [R Packages](https://r-pkgs.org/) (Wickham and Bryan)
   * `testthat` tests, in `tests/testthat`
   * `pkgdown` web site (via [pkgdown.extras](https://github.com/HenrikBengtsson/pkgdown.extras), extensions to allow PDF vignettes); trigger manual builds [here](https://github.com/lme4/lme4/actions/workflows/pkgdown.yaml)
* GitHub 
   * primary development repository
   * [issues](https://github.com/lme4/lme4/issues)
   * testing on [GitHub actions](https://github.com/lme4/lme4/actions) (activated by specifying "[run ci]" at the end of a commit message)
   * [pull requests](https://github.com/lme4/lme4/pulls) are welcome, but please open a discussion as an issue first
