# Old README

Information for installing ancient versions of `lme4`, on ancient versions of R.

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
