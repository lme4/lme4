## .... CLASSES ........................................................

##
## Objects inheriting from virtual class "Covariance" represent
## 'nc'-by-'nc' covariance matrices 'S' using 'par', a numeric vector
## of length at most nc*(nc+1)/2 storing parameters.
##
## Subclasses must define the mapping 'F' from 'par' to 'theta', the
## numeric vector of length nc*(nc+1)/2 storing in column-major order
## the lower triangular entries of the lower triangular Cholesky factor
## of 'S'.
##
##     F(par) = theta := L[lower.tri(L, diag = TRUE)]
##                  L := t(chol(S))
##

setClass("Covariance",
         contains = "VIRTUAL",
         slots = c(nc = "integer", par = "numeric"),
         prototype = list(nc = 0L),
         validity =
         function (object) {
             if (length(nc <- object@nc) != 1L || is.na(nc) || nc < 0L)
                 gettextf("'%s' is not a non-negative integer",
                          "nc")
             else if (!is.double(object@par))
                 gettextf("type of '%s' is not \"%s\"",
                          "par", "double")
             else TRUE
         })

setClass("Covariance.us",
         contains = "Covariance",
         validity =
         function (object) {
             nc <- object@nc
             if (length(object@par) != (nc * (nc - 1L)) %/% 2L + nc)
                 gettextf("length of '%s' is not %s",
                          "par", "nc*(nc+1)/2")
             else TRUE
         })

setClass("Covariance.diag",
         contains = "Covariance",
         slots = c(hom = "logical"),
         prototype = list(hom = FALSE),
         validity =
         function (object) {
             nc <- object@nc
             if (length(hom <- object@hom) != 1L || is.na(hom))
                 gettextf("'%s' is not %s or %s",
                          "hom", "TRUE", "FALSE")
             else if (length(par <- object@par) !=
                      if (hom) nc > 0L else nc) {
                 if (hom)
                     gettextf("length of '%s' is not %d",
                              "par", nc > 0L)
                 else gettextf("length of '%s' is not '%s'",
                               "par", "nc")
             }
             else TRUE
         })

setClass("Covariance.cs",
         contains = "Covariance",
         slots = c(hom = "logical"),
         prototype = list(hom = FALSE),
         validity =
         .fn <-
         function (object) {
             nc <- object@nc
             if (length(hom <- object@hom) != 1L || is.na(hom))
                 gettextf("'%s' is not %s or %s",
                          "hom", "TRUE", "FALSE")
             else if (length(par <- object@par) !=
                      (if (hom) nc > 0L else nc) + (nc > 1L))
                 gettextf("length of '%s' is not %s",
                          "par", if (hom) "(nc>0)+(nc>1)" else "nc+(nc>1)")
             else TRUE
         })

setClass("Covariance.ar1",
         contains = "Covariance",
         slots = c(hom = "logical"),
         prototype = list(hom = FALSE),
         validity = .fn)

rm(.fn)


## .... GENERIC FUNCTIONS ..............................................

## Unused and perhaps not needed
getLambda <- getLambdat <-
function (object) {
    ## stopifnot(is(object, "Covariance"))
    nc <- object@nc
    ans <- matrix(0, nc, nc)
    if (nc > 0L) {
        i <- sequence.default(
            from = seq.int(from = 1L, by = nc + 1L, length.out = nc),
            by = .BY.,
            nvec = nc:1L)
        ans[i] <- getTheta(object)
    }
    ans
}
body(getLambda ) <-
    do.call(substitute, list(body(getLambda ), list(.BY. = 1L)))
body(getLambdat) <-
    do.call(substitute, list(body(getLambdat), list(.BY. = quote(nc))))

getPar <-
function (object) {
    ## stopifnot(is(object, "Covariance"))
    object@par
}

setPar <-
function (object, value) {
    ## stopifnot(is(object, "Covariance"))
    if (!is.double(value))
        stop(gettextf("type of '%s' is not \"%s\"",
                      "value", "double"),
             domain = NA)
    if (length(value) != length(object@par))
        stop(gettextf("length of '%s' is not equal to length of '%s'",
                      "value", "par"),
             domain = NA)
    object@par <- value
    object
}

setGeneric("getTheta",
           function (object, ...)
               standardGeneric("getTheta"))

setGeneric("setTheta",
           function (object, value, ...)
               standardGeneric("setTheta"))

setGeneric("getLower",
           function (object, ...)
               standardGeneric("getLower"))


## .... METHODS ........................................................

setMethod("initialize",
          c(.Object = "Covariance.us"),
          function (.Object, nc, par, ...) {
              if (missing(par) && !missing(nc) &&
                  is.integer(nc) && length(nc) == 1L && !is.na(nc) &&
                  nc >= 0L) {
                  .Object@nc <- nc
                  .Object@par <- `[<-`(double((nc * (nc - 1L)) %/% 2L + nc),
                                       cumsum(c(if (nc > 0L) 1L, if (nc > 1L) nc:2L)),
                                       1)
                  .Object
              }
              else callNextMethod()
          })

setMethod("initialize",
          c(.Object = "Covariance.diag"),
          function (.Object, nc, par, hom = FALSE, ...) {
              if (missing(par) && !missing(nc) &&
                  is.integer(nc) && length(nc) == 1L && !is.na(nc) &&
                  nc >= 0L &&
                  is.logical(hom) && length(hom) == 1L && !is.na(hom)) {
                  .Object@nc <- nc
                  .Object@hom <- hom
                  .Object@par <- rep(1, if (hom) nc > 0L else nc)
                  .Object
              }
              else callNextMethod()
          })

setMethod("initialize",
          c(.Object = "Covariance.cs"),
          .fn <-
          function (.Object, nc, par, hom = FALSE, ...) {
              if (!missing(par) && !missing(nc) &&
                  is.integer(nc) && length(nc) == 1L && !is.na(nc) &&
                  nc >= 0L &&
                  is.logical(hom) && length(hom) == 1L && !is.na(hom)) {
                  .Object@nc <- nc
                  .Object@hom <- hom
                  .Object@par <- c(rep(1, if (hom) nc > 0L else nc),
                                   if (nc > 1L) 0)
                  .Object
              }
              else callNextMethod()
          })

setMethod("initialize",
          c(.Object = "Covariance.ar1"),
          .fn)

rm(.fn)

setMethod("getTheta",
          c(object = "Covariance.us"),
          function (object, ...)
              object@par)

setMethod("getTheta",
          c(object = "Covariance.diag"),
          function (object, ...) {
              nc <- object@nc
              `[<-`(double((nc * (nc - 1L)) %/% 2L + nc),
                    cumsum(c(if (nc > 0L) 1L, if (nc > 1L) nc:2L)),
                    object@par)
          })

setMethod("getTheta",
          c(object = "Covariance.cs"),
          function (object, ...) {
              .NotYetImplemented()
          })

setMethod("getTheta",
          c(object = "Covariance.ar1"),
          function (object, ...) {
              .NotYetImplemented()
          })

setMethod("setTheta",
          c(object = "Covariance.us", value = "numeric"),
          function (object, value, pos = 0L, ...) {
              nc <- object@nc
              nt <- (nc * (nc - 1L)) %/% 2L + nc
              if (!is.double(value))
                  stop(gettextf("type of '%s' is not \"%s\"",
                                "value", "double"),
                       domain = NA)
              if (nt > length(value) - pos)
                  stop(gettextf("attempt to read past end of '%s'",
                                "value"),
                       domain = NA)
              object@par <-
              if (nt == length(value))
                  value
              else {
                  i <- seq.int(from = pos + 1L, length.out = nt)
                  value[i]
              }
              object
          })

setMethod("setTheta",
          c(object = "Covariance.diag", value = "numeric"),
          function (object, value, pos = 0L, ...) {
              nc <- object@nc
              nt <- (nc * (nc - 1L)) %/% 2L + nc
              if (!is.double(value))
                  stop(gettextf("type of '%s' is not \"%s\"",
                                "value", "double"),
                       domain = NA)
              if (nt > length(value) - pos)
                  stop(gettextf("attempt to read past end of '%s'",
                                "value"),
                       domain = NA)
              hom <- object@hom
              i <- if (hom) { if (nc > 0L) pos + 1L } else cumsum(c(if (nc > 0L) pos + 1L, if (nc > 1L) nc:2L))
              object@par <- value[i]
              object
          })

setMethod("setTheta",
          c(object = "Covariance.cs", value = "numeric"),
          function (object, value, pos = 0L, ...) {
              .NotYetImplemented()
              object
          })

setMethod("setTheta",
          c(object = "Covariance.ar1", value = "numeric"),
          function (object, value, pos = 0L, ...) {
              .NotYetImplemented()
              object
          })

setMethod("getLower",
          c(object = "Covariance.us"),
          function (object, ...) {
              nc <- object@nc
              `[<-`(rep(-Inf, (nc * (nc - 1L)) %/% 2L + nc),
                    cumsum(c(if (nc > 0L) 1L, if (nc > 1L) nc:2L)),
                    0)
          })

setMethod("getLower",
          c(object = "Covariance.diag"),
          function (object, ...)
              double(length(object@par)))

setMethod("getLower",
          c(object = "Covariance.cs"),
          .fn <-
          function (object, ...) {
              nc <- object@nc
              hom <- object@hom
              c(double(if (hom) nc > 0L else nc),
                if (nc > 1L) -Inf)
          })

setMethod("getLower",
          c(object = "Covariance.ar1"),
          .fn)

rm(.fn)


## .... HELPERS ........................................................

## return a function that maps
##     c(theta_1, ..., theta_k) -> c(par_1, ..., par_k)
mkMkPar <-
function (reCovs) {
    if (all(vapply(reCovs, is, FALSE, "Covariance.us")))
        return(function (theta) theta)
    nc <- vapply(reCovs, slot, 0L, "nc")
    nt <- (nc * (nc - 1L)) %/% 2L + nc
    np <- lengths(lapply(reCovs, getPar), use.names = FALSE)
    snt <- sum(nt)
    snp <- sum(np)
    jt <- split(seq_len(nt), rep(seq_along(nt), nt))
    jp <- split(seq_len(np), rep(seq_along(np), np))
    par <- double(snp)
    ii <- seq_along(reCovs)
    function (theta) {
        for (i in ii)
            par[jp[[i]]] <<-
                getPar(setTheta(reCovs[[i]], theta[jt[[i]]]))
        par
    }
}

## return a function that maps
##     c(par_1, ..., par_k) -> c(theta_1, ..., theta_k)
mkMkTheta <-
function (reCovs) {
    if (all(vapply(reCovs, is, FALSE, "Covariance.us")))
        return(function (par) par)
    nc <- vapply(reCovs, slot, 0L, "nc")
    nt <- (nc * (nc - 1L)) %/% 2L + nc
    np <- lengths(lapply(reCovs, getPar))
    snt <- sum(nt)
    snp <- sum(np)
    jt <- split(seq_len(nt), rep(seq_along(nt), nt))
    jp <- split(seq_len(np), rep(seq_along(np), np))
    theta <- double(snt)
    ii <- seq_along(reCovs)
    function (par) {
        for (i in ii)
            theta[jt[[i]]] <<-
                getTheta(setPar(reCovs[[i]], par[jp[[i]]]))
        theta
    }
}

## update mkReTrms(..., calc.lambdat = FALSE) so that it has components
## 'Lambdat', 'Lind', 'theta', 'lower', 'par', 'reCovs' where
##
##     length(lower) == length(par) <= length(theta)
##
## note that previously there was no 'par' and no 'reCovs' as 'par' and
## 'theta' were identical by construction
upReTrms <-
function (reTrms, spCalls) {
    spNames <- as.character(lapply(spCalls, `[[`, 1L))
    hom1 <- function (spCall) {
        ## FIXME? not evaluating 'hom' in environment of formula
        hom <- spCall$hom
        if (is.null(hom))
            FALSE
        else if (is.logical(hom) && length(hom) == 1L && !is.na(hom))
            hom
        else stop(gettextf("'%s' is not %s or %s in special call",
                           "hom", "TRUE", "FALSE"),
                  domain = NA)
    }
    hom <- vapply(spCalls, hom1, FALSE, USE.NAMES = FALSE)
    nc <- lengths(reTrms$cnms, use.names = FALSE)
    nl <- reTrms$nl
    ncnl <- rep(nc, nl)
    nt <- (nc * (nc - 1L)) %/% 2L + nc
    ## Inhale ...
    Lt.dp <- sequence.default(from = 1L,
                              by = 1L,
                              nvec = ncnl)
    Lt.p <- cumsum(c(0L, Lt.dp))
    Lt.i <- sequence.default(from = rep(cumsum(c(0L, ncnl)[seq_along(ncnl)]), ncnl),
                             by = 1L,
                             nvec = Lt.dp)
    Lind1 <- function (from, nc)
        rep(seq.int(from = from - nc, by = 1L, length.out = nc),
            seq_len(nc)) +
        cumsum(seq.int(from = nc, by = -1L, length.out = nc))[
            sequence.default(from = 1L, by = 1L, nvec = seq_len(nc))]
    Lind <- unlist(rep(.mapply(Lind1,
                               list(from = cumsum(c(1L, nt)[seq_along(nt)]),
                                    nc = nc),
                               NULL),
                       nl),
                   FALSE, FALSE)
    ## Exhale ...
    reCovs <- .mapply(new,
                      list(Class = paste0("Covariance.", spNames),
                           nc = nc,
                           hom = hom),
                      NULL)
    reTrms$Lambdat <- new("dgCMatrix", Dim = rep(length(Lt.dp), 2L),
                          p = Lt.p, i = Lt.i, x = as.double(Lind))
    reTrms$Lind <- Lind
    reTrms$theta <- unlist(lapply(reCovs, getTheta), FALSE, FALSE)
    reTrms$lower <- unlist(lapply(reCovs, getLower), FALSE, FALSE)
    reTrms$par   <- unlist(lapply(reCovs, getPar  ), FALSE, FALSE)
    reTrms$reCovs <- reCovs
    reTrms
}

## use c(theta_1, ..., theta_k) to update par_1, ..., par_k
upReCovs <-
function (reCovs, theta) {
    ## FIXME?  function (reCovs, par) instead?
    pos <- 0L
    for (i in seq_along(reCovs)) {
        elt <- reCovs[[i]]
        nc <- elt@nc
        nt <- (nc * (nc - 1L)) %/% 2L + nc
        reCovs[[i]] <- setTheta(elt, theta, pos)
        pos <- pos + nt
    }
    reCovs
}
