## need this because R can't handle '@CRAN@' magic default
## in non-interactive mode ...
options(repos=c(CRAN="http://probability.ca/cran",
        rforge="http://r-forge.r-project.org",
        bioc="http://www.bioconductor.org/packages/release/bioc"))

source("pkgdepfuns.R")
rr <- getDepends("lme4")
## checkPkg("afex",verbose=TRUE,checkdir="check")
## checkPkg("cplm",verbose=TRUE,checkdir="check")
testresults <- doPkgDeptests("lme4",verbose=TRUE,do_parallel=FALSE)
save("testresults",file="lme4tests_out.RData")
genReport(rr,testresults)

###

if (FALSE) {
    ## playing with results
    L <- load("lme4tests_out.RData")
    rr <- getDepends("lme4")
    genReport(rr,testresults)
    checkPkg("HSAUR2")
    checkPkg("difR",checkdir="check",verbose=TRUE)
}



