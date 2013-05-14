pkgmax <- function(x,y) {
    if (missing(y)) {
        if (length(x)==1) return(x)
        return(Reduce(pkgmax,x))
    }
    if (package_version(x)>package_version(y)) x else y
}


checkPkg <- function(pn,verbose=FALSE,tarballdir="./tarballs",libdir="./library",
                     checkdir=".",skip=FALSE) {
    ## FIXME: check for/document local version of tarball more recent than R-forge/CRAN versions
    ## FIXME: consistent implementation of checkdir
    rforge <- "http://r-forge.r-project.org"  
    if (!exists("availCRAN")) {
        if (verbose) cat("getting list of available packages from CRAN\n")
        availCRAN <<- available.packages()
    }
    if (!exists("availRforge")) {
        if (verbose) cat("getting list of available packages from R-forge\n")
        availRforge <<- available.packages(contriburl=contrib.url(rforge))
    }
    ## FIXME: maybe this is not effective/not what we want to do?
    ##  *don't* want to look at any already-installed packages
    .libPaths(libdir) 
    ## we definitely want this to check for packages in the local library directory;
    ## not sure if we want to check in the rest of the standard library paths
    instPkgs <- installed.packages(lib.loc=libdir)
    if (verbose) cat("checking package",pn,"\n")
    loc <- "none"  ## where is the package coming from?
    if (pn %in% rownames(availRforge)) {
        ## if available at R-forge ... take the R-forge version
        loc <- "Rforge"
        pkginfo <- availRforge[pn,]
    } else if (pn %in% rownames(availCRAN)) {
        loc <- "CRAN"
        pkginfo <- availCRAN[pn,]
    }
    ## FIXME: look for local tarballs???
    locpkg <- list.files(tarballdir,paste0("^",pn))
    if (length(locpkg)>0) {
        locver <- gsub(paste0(pn,"_([-0-9.]+).tar.gz"),"\\1",locpkg)
    } else locver <- NULL
    if (loc=="none" && is.null(locver)) stop("package seems to be unavailable")
    ver <- pkginfo["Version"]  ## FIXME: check that tarball matches latest available???
    if (!is.null(locver)) {
        if (length(locver)>1) {
            locver <- pkgmax(locver)
        }
        ## FIXME: could figure out most recent version?
        if (package_version(locver)>package_version(ver)) {
            ver <- locver
            loc <- "local"
            if (verbose) cat("local package is more recent than CRAN or R-forge\n")
        }
    }
    ## FIXME: can we safely assume if the file is in 'available.packages(Rforge)' that
    ##   the tarball is really there??
    tn <- paste0(pn,"_",ver,".tar.gz")
    if (loc!="local" && !file.exists(tdn <- file.path(tarballdir,tn)))
    {
        if (verbose) cat("downloading tarball\n")
        basepath <- switch(loc,CRAN=contrib.url(getOption("repos")),
                           Rforge=contrib.url(rforge),
                           loc=stop("tarball not available"))
        download.file(file.path(basepath,tn),
                      destfile=tdn)
        ## install suggested packages that aren't already installed
        depList <- lapply(c("Suggests","Depends"),
                         tools:::package.dependencies,
                         x=pkginfo,
                         check=FALSE)
        depList <- unlist(lapply(depList,function(x) {
            if (!is.matrix(x[[1]])) character(0) else x[[1]][,1] }))
        depMiss <- setdiff(depList,c("R",rownames(instPkgs)))
        if (length(depMiss)>0) {
            if (verbose) cat("installing dependencies",depMiss,"\n")
            install.packages(depMiss,lib=libdir)
            rPath <- if (loc=="CRAN") getOption("repos") else c(rforge,getOption("repos"))
            instPkgs <- installed.packages(noCache=TRUE,lib.loc=libdir)  ## update installed package info
        }
    }
    ## must have set check.Renviron here in order for R CMD check to respect libdir
    if (verbose) cat("running R CMD check ...\n")
    unlink(file.path(checkdir,paste0(pn,".Rcheck")))  ## erase existing check directory
    ## FIXME: run R CMD check in checkdir ...
    if (!skip) {
        tt <- system.time(ss <- suppressWarnings(system(paste("R CMD check",
                                                              file.path(tarballdir,tn)),
                                                        intern=TRUE)))
        if (verbose) print(ss)
        stat <- attr(ss,"status")
        ss <- paste0(seq(ss),": ",ss)
        t0 <- tt["elapsed"]
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
colorCode <- function(x) {
    fcol <- ifelse(grepl("skipped",x),"gray",
                   ifelse(grepl("error_depfail",x),"purple",
                          ifelse(grepl("error_[[:alpha:]]+",x),"red",
                                 "blue")))
    paste0("<font style=\"color:",fcol,"\">",x,"</font>")
}

errstrings <- c(error_ex="checking examples \\.\\.\\. ERROR",
                error_depfail="Package (suggested|required) but not available",
                error_install="checking whether package.+can be installed.+ERROR",
                error_vignette="Errors in running code in vignettes")

procError <- function(z,pkgname=NULL,debug=FALSE) {
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
                      diag=tail(readLines(file.path(paste0(pkgname,".Rcheck"),"00install.out"))))
        } else if (any(grepl(errstrings["error_vignette"],z$msg))) {
            m <- list(result="error_vignette",diag=tail(z$msg,3))
        } else m <- list(result="OK",diag="")
    }
    c(L,m)
}

genReport <- function(X,testresults,
                      contact="lme4-authors <at> r-forge.wu-wien.ac.at",
                      verbose=FALSE) {
    require(lme4)  ## for lme4 package version  (FIXME: should be stored with test results!)
    ## FIXME: should store/pull date from test results too
    require(R2HTML)
    
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
    rpt$result <- colorCode(rpt$result)
    rpt <- merge(rpt,X,by.x="pkgname",by.y="X")
    rpt <- rpt[,c("pkgname","location","version","name","result","diag")] ## drop e-mail, reorder
    names(rpt)[4] <- "maintainer"
    rpt$maintainer <- dumbBrackets(rpt$maintainer)

    ## write file
    HTMLInitFile(filename="lme4_compat_report",outdir=".",
                 Title="lme4 downstream package report")
    HTML("<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\">")  ## for special chars in names etc.
    HTML.title("downstream lme4 package test report",HR=1)
    HTML(paste("lme4 version: ",dumbQuotes(packageVersion(sessionInfo()$otherPkgs$lme4))))
    HTML(paste("Test date: ",date()))
    HTML.title("Notes",HR=2)
    HTML("<ul>")
    HTMLli(paste("contact:",dumbBrackets(contact)))
    HTMLli("'error_depfail' results are due to packages I haven't managed to get installed yet")
    ## HTMLli("'error_install' results due to missing dependencies are probably spurious (packages that are installed elsewhere on my machine but not seen during testing")
    HTMLli("'PIRLS step failed' errors are due to known, current GLMM issues (not your fault)")
    HTML("</ul>")
    HTML(rpt,innerBorder=1,sortableDF=TRUE)
    HTMLEndFile()
}
