## need this because R can't handle '@CRAN@' magic default
## in non-interactive mode ...
options(repos=c(CRAN="http://probability.ca/cran",
        rforge="http://r-forge.r-project.org",
        bioc="http://www.bioconductor.org/packages/release/bioc"))

source("pkgdepfuns.R")
rr <- getDepends("lme4")
pkgnotes <- read.csv("lme4_notes.csv")
testresults <- doPkgDeptests("lme4",verbose=TRUE,do_parallel=FALSE)
save("testresults",file="lme4tests_out.RData")
genReport(rr,testresults,extra.info=pkgnotes)

###

if (FALSE) {
    ## playing with results
    L <- load("lme4tests_out.RData")
    rr <- getDepends("lme4")
    xx <- read.csv("pkg_notes.csv")
    genReport(rr,testresults,extra=xx)
    checkPkg("HSAUR2")
    checkPkg("difR",checkdir="check",verbose=TRUE)
}



