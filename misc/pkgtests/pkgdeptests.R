## code for testing lme4 downstream packages
##
## include all downstream packages from CRAN, r-forge:
## packages to check, loaded from package-dependency scan


## these should eventually be default arguments
verbose <- TRUE
## need this because R can't handle '@CRAN@' magic default
## in non-interactive mode ...
options(repos=c(CRAN="http://probability.ca/cran"))
pkg <- "lme4"
do_parallel <- TRUE
testdir <- getwd()
tarballdir <- file.path(testdir,"tarballs")
libdir <-     file.path(testdir,"library")
checkdir <-     file.path(testdir,"check")
reinstall_pkg <- FALSE
locpkg <- "lme4_0.99999911-2.tar.gz"

## FIXME: should get these straight from DESCRIPTION file
pkgdep <- c("Rcpp","RcppEigen","minqa")
instPkgs <- installed.packages(lib.loc=libdir,noCache=TRUE)
pkgdepMiss <- setdiff(pkgdep,c("R",rownames(instPkgs)))
if (length(pkgdepMiss)>0)
    install.packages(pkgdepMiss,lib=libdir)

if (reinstall_pkg) {
    install.packages(locpkg,repos=NULL,lib=libdir)
}
                 
if (verbose) cat("retrieving dependency information\n")
source("http://developer.r-project.org/CRAN/Scripts/depends.R")
## FIXME: should download a local copy of this, for future-proofing

rr <- reverse_dependencies_with_maintainers(pkg)
source("lme4depfuns.R")  ## component/utility functions
## packages to skip (because dependencies don't install etc.
## skippkgs <- "polytomous"
skippkgs <- character(0)

## * must export R_LIBS_SITE=./library before running R CMD BATCH
##   and  make sure that .R/check.Renviron is set
##   (this is done in the 'runtests' script)

##  (repeated from lme4depfuns.R):
##   FIXME: check for/document local version of tarball more recent than R-forge/CRAN versions
##  currently just tries to find most recent version and check it, but could check all versions?
##  FIXME: consistent implementation of checkdir

## FIXME: set up an appropriate makefile structure for this ? (a little tricky if it also depends on
##   checking CRAN/R-forge versions?
##  might to be able to use update.packages() ...
  
## FIXME (?)/warning: need to make sure that latest/appropriate version of lme4 is installed locally ...

## FIXME: why are R2admb, RLRsim, sdtalt, Zelig not getting checked?




suppressWarnings(rm(list=c("availCRAN","availRforge"))) ## clean up

## require(tools)

## make directories ...
dir.create(tarballdir,showWarnings=FALSE)
dir.create(libdir,showWarnings=FALSE)
## want to install additional dependencies etc. out of the way
## to keep original installed base clean, but this may not be feasible
## it would be nice to use tools:::testInstalledPackages(), but I may simply
##  have to do R CMD check

pkgnames <- rr[,"Package"]

names(pkgnames) <- pkgnames ## so results are named
pkgnames <- pkgnames[!pkgnames %in% skippkgs]

if (verbose) {
    cat("packages to test:\n")
    print(unname(pkgnames),quote=FALSE)
}


if (do_parallel) {
    require(parallel)
    Apply <- mclapply
} else Apply <- lapply
## FIXME (maybe): mclapply doesn't work on Windows ?
##  and might hang Ubuntu VM?

## FIXME: not sure this is necessary/functional
testresults <- Apply(pkgnames,function(x) {
    if (verbose) cat("checking package",x,"\n")
    try(checkPkg(x,verbose=TRUE))
})
skipresults <- Apply(skippkgs,function(x) try(checkPkg(x,skip=TRUE,verbose=TRUE)))
testresults <- c(testresults,skipresults)

save("testresults",file="pkgtests_out.RData")

if (FALSE) {
    ## playing with results
    load("pkgtests_out.RData")
    checkPkg("HSAUR2")
}
## names(testresults) <- X$X  ## should be obsolete after next run
genReport(X,testresults)

