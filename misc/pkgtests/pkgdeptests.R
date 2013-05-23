## need this because R can't handle '@CRAN@' magic default
## in non-interactive mode ...
options(repos=c(CRAN="http://probability.ca/cran"))

source("pkgdepfuns.R")
## checkPkg("afex",verbose=TRUE,checkdir="check")
testresults <- doPkgDeptests("lme4",verbose=TRUE)
save("testresults",file="lme4tests_out.RData")
genReport(rr,testresults)

###

if (FALSE) {
    ## playing with results
    L <- load("pkgtests_out.RData")
    checkPkg("HSAUR2")
}



