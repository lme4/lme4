pkg <- so_name <- "lme4"; doUnload <- FALSE; doTest <- TRUE
## for running uninstalled R binary:
## .libPaths(c("/usr/local/lib/R/library","/usr/local/lib/R/site-library"))
## pkg <- so_name <- "RcppEigen"; doUnload <- TRUE; doTest <- TRUE
## need to deal with the fact that DLL name != package name for lme4.0 ...
### pkg <- "lme4.0"; so_name <- "lme4"; doUnload <- TRUE
instPkgs <- as.data.frame(installed.packages(),stringsAsFactors=FALSE)
Unload <- function() {
    ld <- library.dynam()
    pnames <- sapply(ld,"[[","name")
    names(ld) <- pnames
    lp <- gsub("/libs/.*$","",ld[[so_name]][["path"]])
    cat("unloading from",lp,"\n")
    library.dynam.unload(so_name, lp)
}
Detach <- function() {
    detach(paste0("package:",pkg),character.only=TRUE,unload=TRUE)
    if (doUnload) Unload()
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
if (FALSE) {
  ## FIXME: disabled test for now
for (i in 1:6) {
    cat("Attempt #",i,"\n",sep="")
    cat("loading",pkg,"\n")
    library(pkg,character.only=TRUE)
    tmpf()
    test()
    cat("detaching",pkg,"\n")
    Detach()
    cat("loading nlme\n")
    library("nlme")
    tmpf()
    detach("package:nlme",unload=TRUE)
    cat("detaching nlme\n")
}
}
