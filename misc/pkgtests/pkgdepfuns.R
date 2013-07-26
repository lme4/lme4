## TO DO
## * store pkg version tested, timestamp
## * include suggests/depends/etc. in genReport, or allow subsetting
## * table of outcomes in genReport?

require("tools")
require("plyr") ## for rename()
## originally downloaded from http://developer.r-project.org/CRAN/Scripts/depends.R
## modified (BMB) to include package dependency type; return results as data frame with package rownames
reverse_dependencies_with_maintainers <-
function(packages, which = c("Depends", "Imports", "LinkingTo"),
         cols=c("Package", "Version", "Maintainer"),         
         recursive = FALSE, getDepType=TRUE)
{
    contrib.url(getOption("repos")["CRAN"],type="source") # trigger chooseCRANmirror() if required
    if (length(packages)>1 && getDepType) stop("can't do depType for >1 package")
    description <- sprintf("%s/web/packages/packages.rds",
                           getOption("repos")["CRAN"])
    con <- if(substring(description, 1L, 7L) == "file://")
        file(description, "rb")
    else
        url(description, "rb")
    on.exit(close(con))
    db <- readRDS(gzcon(con))
    rownames(db) <- NULL
    rdepends <- package_dependencies(packages, db, which,
                                     recursive = recursive,
                                     reverse = TRUE)
    rdepends <- sort(unique(unlist(rdepends)))
    pos <- match(rdepends, db[, "Package"], nomatch = 0L)
    d <- data.frame(db[pos,cols],stringsAsFactors=FALSE)
    rownames(d) <- d$Package
    if (getDepType) {
        getType <- function(r) {
            (names(r)[grep(pattern=paste0("(^|[ ,]|\\n)",packages,"([ ,]|\\n|$)"),r)])[1]
        }
        depType <- apply(db[pos, which],1,getType)
        d$depType <- depType
    }
    d
}

getDepends <- function(pkg="lme4",verbose=FALSE, getSuggests=TRUE) {    
    if (verbose) cat("retrieving dependency information\n")
    w <-  c("Depends", "Imports", "LinkingTo")
    if (getSuggests) w <- c(w,"Suggests")
    reverse_dependencies_with_maintainers(pkg,which=w)
}

pkgmax <- function(x,y) {
    if (missing(y)) {
        if (length(x)==1) return(x)
        return(Reduce(pkgmax,x))
    }
    if (package_version(x)>package_version(y)) x else y
}

checkPkg <- function(pn,verbose=FALSE,
                     tarballdir="./tarballs",libdir="./library",
                     checkdir=".",skip=FALSE,
                     check_time=TRUE,
                     upstreamPkg="lme4")
{

    ## expand paths to protect against setwd() for R CMD check
    tarballdir <- normalizePath(tarballdir)
    libdir <- normalizePath(libdir)

    reposURL <- c(CRAN=unname(getOption("repos")["CRAN"]),
                   rforge="http://r-forge.r-project.org",
                   bioconductor="http://www.bioconductor.org/packages/release/bioc")

    if (!exists("availList")) availList <<- list()
    for (i in names(reposURL)) {
        if (is.null(availList[[i]]) || nrow(availList[[i]])==0) {
            if (verbose) cat("getting list of available packages from ",i,"\n")
            availList[[i]] <<- available.packages(contriburl=contrib.url(reposURL[[i]],type="source"),
                                                  type="source")
        }
    }
    .libPaths(libdir)
    ## should include both system files and libdir
    instPkgs <- installed.packages()
    if (verbose) cat("checking package",pn,"\n")
    loc <- "none"  ## where is the package coming from?
    for (i in rev(names(availList))) {  ## check in BACKWARD order (bioc, rforge, CRAN)
        if (pn %in% rownames(availList[[i]])) {
            loc <- i
            pkginfo <- availList[[i]][pn,]
            break
        }
    }
    locpkg <- list.files(tarballdir,paste0("^",pn,"_[0-9]+"))
    if (length(locpkg)>0) {
        locver <- gsub(paste0(pn,"_([-0-9.]+).tar.gz"),"\\1",locpkg)
    } else locver <- NULL
    if (loc=="none" && is.null(locver)) stop("package seems to be unavailable")
    ver <- pkginfo["Version"]  ## FIXME: check that tarball matches latest available???
    if (!is.null(locver)) {
        if (length(locver)>1) {
            locver <- pkgmax(locver)
        }
        if (package_version(locver)>package_version(ver)) {
            ver <- locver
            loc <- "local"
            if (verbose) cat("local package is more recent than CRAN or R-forge\n")
        }
    }
    tn <- paste0(pn,"_",ver,".tar.gz")
    if (loc!="local" && !file.exists(tdn <- file.path(tarballdir,tn)))
    {
        if (verbose) cat("downloading tarball\n")
        basepath <- switch(loc,CRAN=contrib.url(reposURL["CRAN"],type="source"),
                           rforge=contrib.url(reposURL["rforge"],type="source"),
                           bioconductor=contrib.url(reposURL["bioconductor"],
                               type="source"),
                           loc=stop("tarball not available"))
        download.file(file.path(basepath,tn),
                      destfile=tdn)
    }
    ## install suggested packages that aren't already installed
    ## must have set R_LIBS, R_LIBS_SITE, R_LIBS_USER
    ##    in order to match R CMD check settings
    depList <- lapply(c("Suggests","Depends","Imports"),
                      tools:::package.dependencies,
                      x=pkginfo,
                      check=FALSE)
    depList <- unlist(lapply(depList,function(x) {
        if (!is.matrix(x[[1]])) character(0) else x[[1]][,1] }))
    depMiss <- setdiff(depList,c("R",rownames(instPkgs)))
    if (length(depMiss)>0) {
        if (verbose) cat("installing dependencies",depMiss,"\n")
        install.packages(depMiss,lib=libdir,dependencies=TRUE,type="source")
        ## FIXME: not used???
        rPath <- if (loc=="CRAN") reposURL["CRAN"] else reposURL
        instPkgs <- installed.packages(noCache=TRUE,lib.loc=libdir)  ## update installed package info
   }

    ## must have set check.Renviron here in order for R CMD check to respect libdir
    newer_check <- FALSE
    curCheckdir <- file.path(checkdir,paste0(pn,".Rcheck"))
    if (file.exists(curCheckdir)) {
        checktime <- file.info(curCheckdir)["mtime"]  ## check time
        tbinfo <- file.info(file.path(tarballdir,tn)) ## tarball time
        tbtime <- tbinfo["mtime"]
        upinfo <- file.info(file.path(libdir,upstreamPkg)) ## upstream pkg time
        uptime <- upinfo["mtime"]
        newer_check <- (checktime>tbtime && checktime>uptime)
        zero_tb <- tbinfo$size==0
        if (!check_time || !newer_check || zero_tb) unlink(curCheckdir)
    }
    if (!skip) {
        if (check_time && newer_check) {
            if (verbose) cat("check more recent than tarball, skipping\n")
            ss <- readLines(file.path(curCheckdir,"00check.log"))
            t0 <- NA
            stat <- NULL
        } else {
            if (verbose)
                cat("running R CMD check ...\n")
            setwd(checkdir)
            tt <- system.time(ss <- suppressWarnings(system(paste("R CMD check",
                                                                  file.path(tarballdir,tn)),
                                                            intern=TRUE)))
            if (verbose) print(ss)
            stat <- attr(ss,"status")
            ss <- paste0(seq(ss),": ",ss)
            t0 <- tt["elapsed"]
            setwd("..")
        }
    } else {
        stat <- "skipped"
        t0 <- NA
        msg <- ""
        ss <- ""
    }
    list(status=stat,msg=ss,time=t0,location=loc,version=ver)
}


dumbBrackets <- function(x) {
    gsub("<","&lt;",
         gsub(">","&gt",x))
}
dumbQuotes <- function(x) {
    gsub("[“”]","\"",
         gsub("[‘’]","'",x))
}
## old/obsolete
## colorCode <- function(x) {
##     fcol <- ifelse(grepl("skipped",x),"gray",
##                    ifelse(grepl("error_depfail",x),"purple",
##                           ifelse(grepl("error_[[:alpha:]]+",x),"red",
##                                  "blue")))
##     paste0("<font style=\"color:",fcol,"\">",x,"</font>")
## }
## test <- c("skipped","OK","error_depfail","error_other")
colorCode <- function(strvec,
                      colCodes=c(gray="skipped",purple="error_depfail",red="error_[[:alpha:]]",blue=NA)) {
    colvec <- names(colCodes)
    m <- sapply(colCodes,grepl,x=strvec)
    otherVal <- names(colCodes)[is.na(colCodes)]
    tmpf <- function(x) {
        if (sum(w <- which(na.omit(x)==1))==0) otherVal else names(x[w[1]])
    }
    fcol <- apply(m,1,tmpf)
    paste0("<font style=\"color:",fcol,"\">",strvec,"</font>")
}    

errstrings <- c(error_ex="checking examples \\.\\.\\. ERROR",
                error_depfail="Package (suggested|required) but not available",
                error_install="checking whether package.+can be installed.+ERROR",
                error_vignette="Errors in running code in vignettes")

procError <- function(z,pkgname=NULL,debug=FALSE,checkdir="check") {
    if (is.null(loc <- z$location)) loc <- ""
    if (is.null(ver <- z$version)) ver <- ""
    L <- list(pkgname=pkgname,location=loc,version=ver)
    if (!is.null(z$status) && z$status=="skipped") {
        m <- list(result="skipped",diag="")
    } else {
        if (any(grepl(errstrings["error_ex"],z$msg))) {
            m <- list(result="error_examples",diag=tail(z$msg,3))
        } else if (any(grepl(errstrings["error_depfail"],z$msg))) {
            m <- list(result="error_depfail",diag=grep(errstrings["error_depfail"],z$msg,value=TRUE))
        } else if (any(grepl(errstrings["error_install"],z$msg))) {
            m <- list(result="error_install",
                      diag=tail(readLines(file.path(checkdir,paste0(pkgname,".Rcheck"),"00install.out"))))
        } else if (any(grepl(errstrings["error_vignette"],z$msg))) {
            m <- list(result="error_vignette",diag=tail(z$msg,3))
        } else m <- list(result="OK",diag="")
    }
    c(L,m)
}

errLevels <- c(paste("error",
                     c("depfail","examples","install","vignette"),sep="_"),
               "OK")
genReport <- function(depmatrix,      ## results of reverse_dependencies_with_maintainer()
                      testresults,   ## list of packages with elements status, msg, time, location, version
                      contact="lme4-authors <at> r-forge.wu-wien.ac.at",
                      pkg="lme4",
                      outfn=paste0(pkg,"_compat_report"),
                      verbose=FALSE,
                      extra.info=NULL,
                      sortCols=c("result","pkgname")) {
    require(pkg,character.only=TRUE)  ## for package version  (FIXME: should be stored with test results!)
    ## FIXME: should store/pull date from test results too
    if (!require("R2HTML")) {
        ## auto-install because we may be missing it in the test environment ...
        install.packages("R2HTML"); library("R2HTML")
    }
    
    isOK <- !sapply(testresults,inherits,what="try-error")
    tOK <- testresults[isOK]

    if (FALSE) {
        ## testing
        for (i in seq_along(tOK)) {
            print(procError(tOK[[i]],names(tOK)[[i]]))
            scan()
        }
    }
    rpt <- mapply(procError,tOK,names(tOK),SIMPLIFY=FALSE)
    rpt <- lapply(rpt,
                  function(x) {
                      x$diag <- paste(dumbQuotes(x$diag),collapse="<br>")
                      data.frame(x,stringsAsFactors=FALSE)
                  })
    rpt <- do.call(rbind,rpt)
    ## add info from notes, rr
    rpt <- merge(rpt,as.data.frame(depmatrix),by.x="pkgname",by.y="Package")
    ## table of results by package status
    sumtab <- with(rpt,table(result,depType))
    rpt <- rpt[,c("pkgname",
                  "depType","location","version",
                  "Maintainer","result","diag")] ## drop e-mail, reorder
    rpt <- rename(rpt,c(Maintainer="maintainer"))
    if (!is.null(extra.info))
        rpt <- merge(rpt,extra.info,by="pkgname",all.x=TRUE)
    ## mess with ordering by result *before* altering result!
    rpt <- rpt[do.call(order,rpt[sortCols]),]
    ## HTML table formatting
    rpt$maintainer <- dumbBrackets(rpt$maintainer)
    rpt$result <- colorCode(as.character(rpt$result))
    rpt$depType <- colorCode(rpt$depType,
       colCodes=c(blue="Depends",green="Suggests",purple="Imports",red=NA))
    ############# now write file
    title <- paste0(pkg,": downstream package report")
    HTMLInitFile(filename=outfn,outdir=".",
                 Title=title)
    HTML("<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\">")  ## for special chars in names etc.
    HTML.title(title,HR=1)
    HTML(paste(pkg,"version:",dumbQuotes(packageVersion(sessionInfo()$otherPkgs$lme4))))
    HTML(paste("Test date: ",date()))
    HTML.title("Notes",HR=2)
    HTML("<ul>")
    HTMLli(paste("contact:",dumbBrackets(contact)))
    HTMLli("error_depfail indicates a dependency problem")
    ## HTMLli("'error_install' results due to missing dependencies are probably spurious (packages that are installed elsewhere on my machine but not seen during testing")
    HTML("</ul>")
    HTML("<hr>")
    HTML(sumtab)
    HTML("<hr>")
    HTML(rpt,innerBorder=1,sortableDF=TRUE)
    HTMLEndFile()
    outfn
}

doPkgDeptests <- function(pkg="lme4",
                          do_parallel=TRUE,
                          testdir=getwd(),
                          tarballdir=file.path(testdir,"tarballs"),
                          libdir=file.path(testdir,"library"),
                          checkdir=file.path(testdir,"check"),
                          pkg_tarball=NULL,
                          skippkgs=character(0),
                          verbose=TRUE) {

    if (!file.exists(libdir)) dir.create(libdir,showWarnings=FALSE)
    if (!file.exists(checkdir)) dir.create(checkdir,showWarnings=FALSE)
    if (!file.exists(tarballdir)) dir.create(tarballdir,showWarnings=FALSE)

    ## install focal package
    ##   expects recent source tarball in working directory

    ## FIXME: package dependencies
    ##   lme4-specific; should get these straight from DESCRIPTION file
    pkgdep <- c("Rcpp","RcppEigen","minqa")
    if (missing(pkg_tarball) && is.null(pkg_tarball)) {
         pkg_tarball <- list.files(pattern=paste0(pkg,".*.tar.gz"))
         if (length(pkg_tarball)==0) warning("can't find package tarball: not re-installing focal package")
    }
    instPkgs <- installed.packages(lib.loc=libdir,noCache=TRUE)
    pkgdepMiss <- setdiff(pkgdep,c("R",rownames(instPkgs)))
    if (length(pkgdepMiss)>0)
        install.packages(pkgdepMiss,lib=libdir, type="source")
    if (!is.null(pkg_tarball)) {
        tb0times <- file.info(pkg_tarball)$mtime
        pkg_tarball <- pkg_tarball[which.max(tb0times)]
        tb0time <- max(tb0times)
        pkg_inst <- file.exists(file.path(libdir,pkg))
        pkgtime <- if (!pkg_inst) -Inf else {
            file.info(file.path(libdir,pkg))$mtime
        }
        if (tb0time>pkgtime)
            install.packages(pkg_tarball,repos=NULL,lib=libdir,type="source")
    }
    ## * must export R_LIBS_SITE=./library before running R CMD BATCH
    ##   and  make sure that .R/check.Renviron is set
    ##   (this is done by setTestEnv, called from 'runCheck')

    ##  FIXME: consistent implementation of checkdir

    ## FIXME: set up an appropriate makefile structure for this ? (a little tricky if it also depends on
    ##   checking CRAN/R-forge versions?
    ##  might to be able to use update.packages() ...

    suppressWarnings(rm("availList"))
    ## suppressWarnings(rm(list=c("availCRAN","availRforge"))) ## clean up

    ## want to install additional dependencies etc. out of the way
    ## to keep original installed base clean, but this may not be feasible
    ## it would be nice to use tools:::testInstalledPackages(), but I may simply
    ##  have to do R CMD check

    rr <- getDepends(pkg,verbose)
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
    ## FIXME (maybe): mclapply doesn't work on Windows??

    testresults <- Apply(pkgnames,function(x) {
        ## if (verbose) cat("checking package",x,"\n")  ## redundant
        try(checkPkg(x,verbose=TRUE,checkdir=checkdir))
    })
    skipresults <- Apply(skippkgs,function(x) try(checkPkg(x,skip=TRUE,verbose=TRUE)))
    testresults <- c(testresults,skipresults)
    testresults
}
