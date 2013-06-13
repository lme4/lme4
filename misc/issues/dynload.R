## this is the more complex version of the code for testing/exercising
## https://github.com/lme4/lme4/issues/35
## see also ../../tests/dynload.R for the simpler version
pkg <- so_name <- "lme4"  ## package to test
doUnload <- FALSE         ## manual unloading?
doTest <- TRUE            ## run test between load/unload cycles?
use_devtools <- FALSE     ## use load_all/unload rather than library/detach; only works if wd=lme4 directory
otherpkg <- "digest"      ## package to load/unload between focal package loading/unloading
## for running e.g. from uninstalled binary of R with full valgrind:
##    .libPaths(c("/usr/local/lib/R/library","/usr/local/lib/R/site-library"))
## -- alternative options --
## pkg <- so_name <- "RcppEigen"; doUnload <- TRUE; doTest <- TRUE
## need to deal with the fact that DLL name != package name for lme4.0 ...
### pkg <- "lme4.0"; so_name <- "lme4"; doUnload <- TRUE
instPkgs <- as.data.frame(installed.packages(),stringsAsFactors=FALSE)
Load <- function() {
    if (use_devtools) {
        require("devtools")
        load_all(reset=TRUE)
    } else library(pkg,character.only=TRUE)
}
Unload <- function() {
    ld <- library.dynam()
    pnames <- sapply(ld,"[[","name")
    names(ld) <- pnames
    lp <- gsub("/libs/.*$","",ld[[so_name]][["path"]])
    cat("unloading from",lp,"\n")
    library.dynam.unload(so_name, lp)
}
Detach <- function() {
    if (use_devtools) {
        unload()
        ## don't really need this -- just load_all(reset=TRUE)
        ##  should work -- but this helps maintain the same sequence
    } else {
        detach(paste0("package:",pkg),character.only=TRUE,unload=TRUE)
        if (doUnload) Unload()
    }
}
tmpf <- function() {
    g <- getLoadedDLLs()
    lnames <- names(g)[is.na(instPkgs[names(g),"Priority"])]
    cat("loaded DLLs:",lnames,"\n")
    g <- g[na.omit(match(c(so_name,"nlme"),names(g)))]
    class(g) <- "DLLInfoList"
    g
}
test <- function() {
    if (doTest) {
        if (pkg %in% c("lme4","lme4.0")) {
            fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy,
                        devFunOnly=TRUE)
        }
        if (pkg=="RcppEigen") {
            data(trees, package="datasets")
            mm <- cbind(1, log(trees$Girth))   # model matrix
            y  <- log(trees$Volume)            # response
            ## bare-bones direct interface
            flm <- fastLmPure(mm, y)
        }
    }
}

## erase existing .onUnload hook
## nullfun <- function(libpath) { cat("fake unload\n") }
## assignInNamespace(".onUnload",nullfun,"lme4")

for (i in 1:6) {
    cat("Attempt #",i,"\n",sep="")
    cat("loading",pkg,"\n")
    Load()
    tmpf()
    test()
    cat("detaching",pkg,"\n")
    Detach()
    cat("loading",otherpkg,"\n")
    library(otherpkg,character.only=TRUE)
    tmpf()
    detach(paste0("package:",otherpkg),character.only=TRUE,unload=TRUE)
    cat("detaching",otherpkg,"\n")
}

## lme4/digest with unload=FALSE: OK
cycle1 <- function() {
    library(lme4); capture.output(example(lmer)); detach(package:lme4)
    library(digest); detach(package:digest,unload=TRUE)
}
invisible(replicate(10,cycle1()))
## RcppEigen/digest with unload=FALSE: OK
cycle2 <- function() {
    library(RcppEigen); capture.output(example(fastLm)); detach(package:RcppEigen)
    library(digest); detach(package:digest,unload=TRUE)
}
invisible(replicate(10,cycle2()))
## RcppEigen/digest with unload=TRUE: OK
cycle3 <- function() {
    library(RcppEigen); capture.output(example(fastLm)); detach(package:RcppEigen,unload=TRUE)
    library(digest); detach(package:digest,unload=TRUE)
}
invisible(replicate(10,cycle3()))
## lmer/digest with unload=TRUE: NOT OK
cycle4 <- function() {
    library(lme4); capture.output(example(lmer)); detach(package:lme4,unload=TRUE)
    library(digest); detach(package:digest,unload=TRUE)
}
invisible(replicate(10,cycle4()))
## lmer/digest with unload=TRUE but unload hook removed: OK
aa <- function() {
    fun0 <- function(libpath) { cat("fake unload\n") }
    assignInNamespace(".onUnload",fun0,"lme4")
}
cycle5 <- function() {
    library(lme4); capture.output(example(lmer)); aa(); detach(package:lme4,unload=TRUE)
    library(digest); detach(package:digest,unload=TRUE)
}
invisible(replicate(10,cycle5()))


### NOTES:
###   can get a crash with just merPredD$new(...)
##    if we don't initializePtr() everything is OK
##    the proximal problem is .Call(merPredDCreate,...)
##    now checking to see if the non-trivial stuff on ll. 52ff of predModule.cpp
##   (RX.compute, updateLamtUt, analyzePattern, etc.) are necessary -- they're not.
## so ... we should be able to make lme4_Lobotomized (with apologies to all mental
##  patients) which contains _only_ AllClass.R, predModule.cpp, external.cpp
##  i.e. just the definitions of merPredD class and its creation
## ... and this should cause the same problem ...

## related? old issue with finalization of ref class objects -- the
##  errors referred to herein don't cause segfaults, but do still happen
##  with current R ...
## https://stat.ethz.ch/pipermail/r-devel/2011-December/062917.html
