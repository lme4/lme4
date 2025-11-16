## .... CLASSES ........................................................

##
## Objects inheriting from virtual class "Covariance" represent
## 'nc'-by-'nc' covariance matrices 'S' using 'par', a vector storing
## nc*(nc+1)/2 or fewer parameters.
##
## Subclasses must define a mapping
##
##     par [double]  ---->  ( theta [double] ,  thetaIndex [integer] )
##
## such that theta[thetaIndex] stores in column-major order the
## structurally nonzero entries of Lambdat := chol(S).  By convention,
## subclasses will define 'theta' as storing in *row-major* order only
## the *independent* structurally nonzero entries of 'Lambdat', keeping
## the length of 'theta' to a minimum.
##
## MJ:
## Why do we change the storage order, when doing so greatly complicates
## determination of 'thetaIndex'?
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

setGeneric("getPar",
           function (object)
               standardGeneric("getPar"))

setGeneric("getParLength",
           function (object)
               standardGeneric("getParLength"))

setGeneric("getParNames",
           function (object, cnm, gnm, prf = NULL, old = TRUE)
               standardGeneric("getParNames"),
           signature = c("object", "cnm", "gnm"))

setGeneric("setPar",
           function (object, value, pos = 0L)
               standardGeneric("setPar"))

setGeneric("getTheta",
           function (object)
               standardGeneric("getTheta"))

setGeneric("getThetaLength",
           function (object)
               standardGeneric("getThetaLength"))

setGeneric("getThetaNames",
           function (object, cnm, gnm, prf = NULL, old = TRUE)
               standardGeneric("getThetaNames"),
           signature = c("object", "cnm", "gnm"))

setGeneric("getThetaIndex",
           function (object)
               standardGeneric("getThetaIndex"))

setGeneric("getThetaIndexLength",
           function (object)
               standardGeneric("getThetaIndexLength"))

setGeneric("setTheta",
           function (object, value, pos = 0L)
               standardGeneric("setTheta"),
           signature = c("object", "value"))

setGeneric("getLower",
           function (object)
               standardGeneric("getLower"))

setGeneric("getUpper",
           function (object)
               standardGeneric("getUpper"))

setGeneric("getLambda",
           function (object)
               standardGeneric("getLambda"))

setGeneric("getLambdat.dp",
           function (object)
               standardGeneric("getLambdat.dp"))

setGeneric("getLambdat.i",
           function (object)
               standardGeneric("getLambdat.i"))

setGeneric("getVC",
           function (object)
               standardGeneric("getVC"))

setGeneric("getVCNames",
           function (object, cnm, gnm, prf = NULL, old = TRUE)
               standardGeneric("getVCNames"),
           signature = c("object", "cnm", "gnm"))

setGeneric("setVC",
           function (object, vcomp, ccomp)
               standardGeneric("setVC"))

setGeneric("getProfPar",
           function (object, profscale, sc = NULL)
               standardGeneric("getProfPar"),
           signature = "object")

setGeneric("setProfPar",
           function (object, profscale, sc = NULL, value, pos = 0L)
               standardGeneric("setProfPar"),
           signature = c("object", "value"))

setGeneric("getProfLower",
           function (object, profscale, sc = NULL)
               standardGeneric("getProfLower"),
           signature = "object")

setGeneric("getProfUpper",
           function (object, profscale, sc = NULL)
               standardGeneric("getProfUpper"),
           signature = "object")


## .... METHODS ........................................................

setMethod("initialize",
          c(.Object = "Covariance.us"),
          function (.Object, nc, par, simulate = FALSE, ...) {
              if (missing(par) && !missing(nc) &&
                  is.integer(nc) && length(nc) == 1L && !is.na(nc) &&
                  nc >= 0L) {
                  .Object@nc <- nc
                  .Object@par <-
                  `[<-`((if (simulate) rnorm else double)(
                            (nc * (nc - 1L)) %/% 2L + nc),
                        cumsum(c(if (nc > 0L) 1L, if (nc > 1L) nc:2L)),
                        if (simulate) rlnorm(nc) else 1)
                  .Object
              }
              else callNextMethod()
          })

setMethod("initialize",
          c(.Object = "Covariance.diag"),
          function (.Object, nc, par, hom = FALSE, simulate = FALSE, ...) {
              if (missing(par) && !missing(nc) &&
                  is.integer(nc) && length(nc) == 1L && !is.na(nc) &&
                  nc >= 0L &&
                  is.logical(hom) && length(hom) == 1L && !is.na(hom)) {
                  .Object@nc <- nc
                  .Object@hom <- hom
                  .Object@par <-
                  (if (simulate) rlnorm else function (n) rep(1, n))(
                      if (hom) 1L * (nc > 0L) else nc)
                  .Object
              }
              else callNextMethod()
          })

setMethod("initialize",
          c(.Object = "Covariance.cs"),
          function (.Object, nc, par, hom = FALSE, simulate = FALSE, ...) {
              if (missing(par) && !missing(nc) &&
                  is.integer(nc) && length(nc) == 1L && !is.na(nc) &&
                  nc >= 0L &&
                  is.logical(hom) && length(hom) == 1L && !is.na(hom)) {
                  .Object@nc <- nc
                  .Object@hom <- hom
                  .Object@par <-
                  c((if (simulate) rlnorm else function (n) rep(1, n))(
                        if (hom) 1L * (nc > 0L) else nc),
                    if (nc > 1L) { if (simulate) runif(1L, -1/(nc - 1), 1) else 0 })
                  .Object
              }
              else callNextMethod()
          })

setMethod("initialize",
          c(.Object = "Covariance.ar1"),
          function (.Object, nc, par, hom = FALSE, simulate = FALSE, ...) {
              if (missing(par) && !missing(nc) &&
                  is.integer(nc) && length(nc) == 1L && !is.na(nc) &&
                  nc >= 0L &&
                  is.logical(hom) && length(hom) == 1L && !is.na(hom)) {
                  .Object@nc <- nc
                  .Object@hom <- hom
                  .Object@par <-
                  c((if (simulate) rlnorm else function (n) rep(1, n))(
                        if (hom) 1L * (nc > 0L) else nc),
                    if (nc > 1L) { if (simulate) runif(1L, -1, 1) else 0 })
                  .Object
              }
              else callNextMethod()
          })

setMethod("getPar",
          c(object = "Covariance"),
          function (object)
              object@par)

setMethod("getParLength",
          c(object = "Covariance"),
          function (object)
              length(object@par))

setMethod("getParNames",
          c(object = "Covariance.us", cnm = "character", gnm = "character"),
          getParNames.us <-
          function (object, cnm, gnm, prf = NULL, old = TRUE) {
              nc <- object@nc
              stopifnot(length(cnm) == nc, length(gnm) == 1L)
              inc <- seq_len(nc)
              dec <- seq.int(from = nc, by = -1L, length.out = nc)
              i <- sequence.default(from = inc, nvec = dec)
              j <- rep(inc, dec)
              ans <-
                  if (old)
                      paste0(gnm, ".", cnm[i], ".", cnm[j])
                  else if (is.null(prf))
                      paste0(cnm[i], ".", cnm[j], "|", gnm)
                  else paste0(prf[[2L]], "_", cnm[i], ".", cnm[j], "|", gnm)
              ans[cumsum(c(if (nc > 0L) 1L, if (nc > 1L) nc:2L))] <-
                  if (old)
                      paste0(gnm, ".", cnm)
                  else if (is.null(prf))
                      paste0(cnm, "|", gnm)
                  else paste0(prf[[1L]], "_", cnm, "|", gnm)
              ans
          })

setMethod("getParNames",
          c(object = "Covariance.diag", cnm = "character", gnm = "character"),
          getParNames.diag <-
          function (object, cnm, gnm, prf = NULL, old = TRUE) {
              nc <- object@nc
              stopifnot(length(cnm) == nc, length(gnm) == 1L)
              if (object@hom && nc > 0L)
                  cnm <- "*"
              if (old)
                  paste0(gnm, ".", cnm)
              else if (is.null(prf))
                  paste0(cnm, "|", gnm)
              else paste0(prf[[1L]], "_", cnm, "|", gnm)
          })

setMethod("getParNames",
          c(object = "Covariance.cs", cnm = "character", gnm = "character"),
          .fn <-
          function (object, cnm, gnm, prf = NULL, old = TRUE) {
              nc <- object@nc
              stopifnot(length(cnm) == nc, length(gnm) == 1L)
              if (object@hom && nc > 0L)
                  cnm <- "*"
              if (old)
                  paste0(gnm, ".", c(cnm, if (nc > 1L) "rho"))
              else if (is.null(prf))
                  paste0(c(cnm, if (nc > 1L) "rho"), "|", gnm)
              else c(paste0(prf[[1L]], "_", cnm, "|", gnm),
                     paste0(prf[[2L]], "_", "rho", "|", gnm))
          })

setMethod("getParNames",
          c(object = "Covariance.ar1", cnm = "character", gnm = "character"),
          .fn)

rm(.fn)

setMethod("setPar",
          c(object = "Covariance", value = "numeric"),
          function (object, value, pos = 0L) {
              np <- length(object@par)
              validValuePos(np, value, pos)
              object@par <-
              if (np == length(value))
                  value
              else
                  value[seq.int(from = pos + 1L, length.out = np)]
              object
          })

setMethod("getTheta",
          c(object = "Covariance.us"),
          .fn <-
          function (object)
              object@par)

## L[i, j]/sigma[i] =
##     if (i == j)
##         1
##     else 0
setMethod("getTheta",
          c(object = "Covariance.diag"),
          .fn)

rm(.fn)

## L[i, j]/sigma[i] =
##     if (i == j)
##         sqrt(1 - rho^2 * a[j])
##     else (rho - rho^2 * a[j])/sqrt(1 - rho^2 * a[j])
## where
##     a[    1] := 0
##     a[j + 1] := a[j] + (1 - rho * a[j])^2/(1 - rho^2 * a[j])
setMethod("getTheta",
          c(object = "Covariance.cs"),
          function (object) {
              nc <- object@nc
              par <- object@par
              if (nc <= 1L)
                  return(par)
              rho2 <- (rho <- par[length(par)])^2
              a <- double(nc); a. <- 0
              for (j in 1L:nc) {
                  a[j] <- a.
                  a. <- a. + (1 - rho * a.)^2/(1 - rho2 * a.)
              }
              v. <- sqrt(1 - rho2 * a)
              v <- rbind(v., (rho - rho2 * a)/v.)
              if (object@hom)
                  par[1L] * v[1L:(2L * nc - 1L)]
              else
                  par[sequence.default(from = 1L:nc, nvec = nc:1L)] *
                      rep(v, rbind(1L, (nc - 1L):0L))
          })

## L[i, j]/sigma[i] =
##     if (j == 1)
##         rho^(i - 1)
##     else rho^(i - j) * sqrt(1 - rho^2)
setMethod("getTheta",
          c(object = "Covariance.ar1"),
          function (object) {
              nc <- object@nc
              par <- object@par
              if (nc <= 1L)
                  return(par)
              rho2 <- (rho <- par[length(par)])^2
              v1 <- rho^(0L:(nc - 1L))
              v2 <- rho^(0L:(nc - 2L)) * sqrt(1 - rho2)
              if (object@hom)
                  par[1L] * c(v1, v2)
              else
                  par[sequence.default(from = 1L:nc, nvec = nc:1L)] *
                      c(v1, v2[sequence.default(from = 1L, nvec = (nc - 1L):1L)])
          })

setMethod("getThetaLength",
          c(object = "Covariance.us"),
          .fn <-
          function (object)
              length(object@par))

setMethod("getThetaLength",
          c(object = "Covariance.diag"),
          .fn)

rm(.fn)

setMethod("getThetaLength",
          c(object = "Covariance.cs"),
          .fn <-
          function (object) {
              nc <- object@nc
              if (object@hom)
                  2L * nc - (nc > 0L)
              else (nc * (nc - 1L)) %/% 2L + nc
          })

setMethod("getThetaLength",
          c(object = "Covariance.ar1"),
          .fn)

rm(.fn)

setMethod("getThetaNames",
          c(object = "Covariance.us", cnm = "character", gnm = "character"),
          function (object, cnm, gnm, prf = NULL, old = TRUE)
              getParNames.us(object, cnm, gnm, prf, old))

setMethod("getThetaNames",
          c(object = "Covariance.diag", cnm = "character", gnm = "character"),
          function (object, cnm, gnm, prf = NULL, old = TRUE)
              getParNames.diag(object, cnm, gnm, prf, old))

setMethod("getThetaNames",
          c(object = "Covariance.cs", cnm = "character", gnm = "character"),
          function (object, cnm, gnm, prf = NULL, old = TRUE) {
              if (object@hom) {
                  nc <- object@nc
                  rbind(if (old)
                            paste0(gnm, ".", cnm)
                        else if (is.null(prf))
                            paste0(cnm, "|", gnm)
                        else paste0(prf[[1L]], "_", cnm, "|", gnm),
                        if (old)
                            paste0(gnm, ".", "*", ".", cnm)
                        else if (is.null(prf))
                            paste0("*", ".", cnm, "|", gnm)
                        else paste0(prf[[2L]], "_", "*", ".", cnm, "|", gnm))[seq_len(2L * nc - 1L)]
              }
              else getParNames.us(object, cnm, gnm, prf, old)
          })

setMethod("getThetaNames",
          c(object = "Covariance.ar1", cnm = "character", gnm = "character"),
          function (object, cnm, gnm, prf = NULL, old = TRUE) {
              if (object@hom) {
                  nc <- object@nc
                  i <- if (nc > 1L) c(1L:nc, 2L:nc) else nc
                  j <- if (nc > 1L) rep(1L:2L, nc:(nc - 1L)) else nc
                  k <- seq_len(min(2L, nc))
                  ans <-
                      if (old)
                          paste0(gnm, ".", cnm[i], ".", cnm[j])
                      else if (is.null(prf))
                          paste0(cnm[i], ".", cnm[j], "|", gnm)
                      else paste0(prf[[2L]], "_", cnm[i], ".", cnm[j], "|", gnm)
                  ans[c(1L, nc + 1L)[k]] <-
                      if (old)
                          paste0(gnm, ".", cnm[k])
                      else if (is.null(prf))
                          paste0(cnm[k], "|", gnm)
                      else paste0(prf[[2L]], "_", cnm[k], "|", gnm)
                  ans
              }
              else getParNames.us(object, cnm, gnm, prf, old)
          })

setMethod("getThetaIndex",
          c(object = "Covariance.us"),
          getThetaIndex.us <-
          function (object) {
              snc <- seq_len(nc <- object@nc)
              rep(seq.int(from = 1L - nc, length.out = nc),
                  snc) +
              cumsum(seq.int(from = nc, by = -1L, length.out = nc))[
                  sequence.default(from = 1L, nvec = snc)]
          })

setMethod("getThetaIndex",
          c(object = "Covariance.diag"),
          function (object) {
              nc <- object@nc
              if (object@hom) rep(1L, nc) else seq_len(nc)
          })

setMethod("getThetaIndex",
          c(object = "Covariance.cs"),
          function (object) {
              if (object@hom) {
                  snc <- seq_len(nc <- object@nc)
                  `[<-`(sequence.default(from = 2L, by = 2L, nvec = snc),
                        cumsum(snc),
                        seq.int(from = 1L, by = 2L, length.out = nc))
              }
              else getThetaIndex.us(object)
          })

setMethod("getThetaIndex",
          c(object = "Covariance.ar1"),
          function (object) {
              if (object@hom) {
                  snc <- seq_len(nc <- object@nc)
                  `[<-`(sequence.default(from = seq.int(from = nc + 1L, length.out = nc),
                                         by = -1L,
                                         nvec = snc),
                        1L + cumsum(seq.int(from = 0L, length.out = nc)),
                        snc)
              }
              else getThetaIndex.us(object)
          })

setMethod("getThetaIndexLength",
          c(object = "Covariance.us"),
          function (object) {
              nc <- object@nc
              (nc * (nc - 1L)) %/% 2L + nc
          })

setMethod("getThetaIndexLength",
          c(object = "Covariance.diag"),
          function (object)
              object@nc)

setMethod("getThetaIndexLength",
          c(object = "Covariance.cs"),
          function (object) {
              nc <- object@nc
              (nc * (nc - 1L)) %/% 2L + nc
          })

setMethod("getThetaIndexLength",
          c(object = "Covariance.ar1"),
          function (object) {
              nc <- object@nc
              (nc * (nc - 1L)) %/% 2L + nc
          })

setMethod("setTheta",
          c(object = "Covariance.us", value = "numeric"),
          .fn <-
          function (object, value, pos = 0L) {
              nt <- length(object@par)
              validValuePos(nt, value, pos)
              object@par <-
              if (nt == length(value))
                  value
              else
                  value[seq.int(from = pos + 1L, length.out = nt)]
              object
          })

setMethod("setTheta",
          c(object = "Covariance.diag", value = "numeric"),
          .fn)

rm(.fn)

setMethod("setTheta",
          c(object = "Covariance.cs", value = "numeric"),
          function (object, value, pos = 0L) {
              nc <- object@nc
              hom <- object@hom
              nt <- if (hom) 2L * nc - (nc > 0L) else (nc * (nc - 1L)) %/% 2L + nc
              validValuePos(nt, value, pos)
              object@par <-
              if (nc <= 1L)
                  value[if (nc > 0L) 1L]
              else {
                  ## L[2, 1] = sigma[2] * rho
                  l21 <- value[pos + 1L + 1L]
                  ## L[2, 2] = sigma[2] * sqrt(1 - rho^2)
                  l22 <- value[pos + 1L + if (hom) 2L else nc]
                  if (l21 == 0 && l22 == 0) # FIXME
                      stop(gettextf("boundary case not yet supported by method for '%s' with signature (%s=\"%s\", %s=\"%s\")",
                                    "setTheta", "object", "Covariance.cs", "value", "numeric"),
                           domain = NA)
                  rho <- sign(l21) * 1/sqrt(1 + (l22/l21)^2)
                  sigma <-
                  if (hom)
                      value[pos + 1L]
                  else if (rho == 0)
                      value[cumsum(c(pos + 1L, nc:2L))]
                  else
                      value[seq.int(from = pos + 1L, length.out = nc)]/
                          rep(c(1, rho), c(1L, nc - 1L))
                  c(sigma, rho)
              }
              object
          })

setMethod("setTheta",
          c(object = "Covariance.ar1", value = "numeric"),
          function (object, value, pos = 0L) {
              nc <- object@nc
              hom <- object@hom
              nt <- if (hom) 2L * nc - (nc > 0L) else (nc * (nc - 1L)) %/% 2L + nc
              validValuePos(nt, value, pos)
              object@par <-
              if (nc <= 1L)
                  value[if (nc > 0L) 1L]
              else {
                  ## L[2, 1] = sigma[2] * rho
                  l21 <- value[pos + 1L + 1L]
                  ## L[2, 2] = sigma[2] * sqrt(1 - rho^2)
                  l22 <- value[pos + 1L + nc]
                  if (l21 == 0 && l22 == 0) # FIXME
                      stop(gettextf("boundary case not yet supported by method for '%s' with signature (%s=\"%s\", %s=\"%s\")",
                                    "setTheta", "object", "Covariance.ar1", "value", "numeric"),
                           domain = NA)
                  rho <- sign(l21) * 1/sqrt(1 + (l22/l21)^2)
                  sigma <-
                  if (hom)
                      value[pos + 1L]
                  else if (rho == 0)
                      value[cumsum(c(pos + 1L, nc:2L))]
                  else # avoid underflow of powers of abs(rho)
                      exp(log(abs(value[seq.int(from = pos + 1L, length.out = nc)])) - (0L:(nc - 1L)) * log(abs(rho)))
                  c(sigma, rho)
              }
              object
          })

setMethod("getLower",
          c(object = "Covariance.us"),
          function (object) {
              nc <- object@nc
              ## FIXME: make more scrutable!
              `[<-`(rep(-Inf, length(object@par)),
                    cumsum(c(if (nc > 0L) 1L, if (nc > 1L) nc:2L)),
                    0)
          })

setMethod("getLower",
          c(object = "Covariance.diag"),
          function (object)
              double(length(object@par)))

setMethod("getLower",
          c(object = "Covariance.cs"),
          function (object) {
              nc <- object@nc
              c(double(if (object@hom) 1L * (nc > 0L) else nc),
                if (nc > 1L) -1/(nc - 1L))
          })

setMethod("getLower",
          c(object = "Covariance.ar1"),
          function (object) {
              nc <- object@nc
              c(double(if (object@hom) 1L * (nc > 0L) else nc),
                if (nc > 1L) -1)
          })

setMethod("getUpper",
          c(object = "Covariance.us"),
          function (object)
              rep(Inf, length(object@par)))

setMethod("getUpper",
          c(object = "Covariance.diag"),
          function (object)
              rep(Inf, length(object@par)))

setMethod("getUpper",
          c(object = "Covariance.cs"),
          function (object) {
              nc <- object@nc
              c(rep(Inf, if (object@hom) 1L * (nc > 0L) else nc),
                if (nc > 1L) 1)
          })

setMethod("getUpper",
          c(object = "Covariance.ar1"),
          function (object) {
              nc <- object@nc
              c(rep(Inf, if (object@hom) 1L * (nc > 0L) else nc),
                if (nc > 1L) 1)
          })

setMethod("getLambda",
          c(object = "Covariance.us"),
          function (object) {
              nc <- object@nc
              ans <- matrix(0, nc, nc)
              if (nc > 0L)
              ans[lower.tri(ans, diag = TRUE)] <- object@par
              ans
          })

setMethod("getLambda",
          c(object = "Covariance.diag"),
          function (object) {
              nc <- object@nc
              diag(object@par, nc, nc)
          })

setMethod("getLambda",
          c(object = "Covariance.cs"),
          function (object) {
              nc <- object@nc
              ans <- matrix(0, nc, nc)
              if (nc > 0L)
              ans[lower.tri(ans, diag = TRUE)] <-
                  if (object@hom)
                      rep(getTheta(object), rbind(1L, (nc - 1L):0L)[seq_len(2L * nc - 1L)])
                  else getTheta(object)
              ans
          })

setMethod("getLambda",
          c(object = "Covariance.ar1"),
          function (object) {
              nc <- object@nc
              ans <- matrix(0, nc, nc)
              if (nc > 0L)
              ans[lower.tri(ans, diag = TRUE)] <-
                  if (object@hom)
                      getTheta(object)[sequence.default(from = rep(c(1L, nc + 1L), c(1L, nc - 1L)), nvec = nc:1L)]
                  else getTheta(object)
              ans
          })

setMethod("getLambdat.dp",
          c(object = "Covariance.us"),
          function (object)
              seq_len(object@nc))

setMethod("getLambdat.dp",
          c(object = "Covariance.diag"),
          function (object)
              rep(1L, object@nc))

setMethod("getLambdat.dp",
          c(object = "Covariance.cs"),
          function (object)
              seq_len(object@nc))

setMethod("getLambdat.dp",
          c(object = "Covariance.ar1"),
          function (object)
              seq_len(object@nc))

setMethod("getLambdat.i",
          c(object = "Covariance.us"),
          function (object)
              sequence.default(from = 0L, nvec = seq_len(object@nc)))

setMethod("getLambdat.i",
          c(object = "Covariance.diag"),
          function (object)
              seq.int(from = 0L, length.out = object@nc))

setMethod("getLambdat.i",
          c(object = "Covariance.cs"),
          function (object)
              sequence.default(from = 0L, nvec = seq_len(object@nc)))

setMethod("getLambdat.i",
          c(object = "Covariance.ar1"),
          function (object)
              sequence.default(from = 0L, nvec = seq_len(object@nc)))

setMethod("getVC",
          c(object = "Covariance.us"),
          function (object) {
              nc <- object@nc
              if (nc <= 1L) {
                  vcomp <- object@par
                  ccomp <- double(0L)
              } else {
                  ii <- seq.int(from = 1L, by = nc + 1L, length.out = nc)
                  i0 <- sequence.default(from = seq.int(from = 2L, by = nc + 1L, length.out = nc - 1L),
                                         nvec = (nc - 1L):1L)
                  i1 <- sequence.default(from = ii,
                                         nvec = nc:1L)
                  L <- matrix(0, nc, nc)
                  L[i1] <- object@par
                  S <- tcrossprod(L)
                  vcomp <- sqrt(S[ii])
                  ccomp <- (S/vcomp/rep(vcomp, each = nc))[i0]
              }
              list(vcomp = vcomp, ccomp = ccomp)
          })

setMethod("getVC",
          c(object = "Covariance.diag"),
          function (object)
              list(vcomp = object@par, ccomp = double(0L)))

setMethod("getVC",
          c(object = "Covariance.cs"),
          .fn <-
          function (object) {
              nc <- object@nc
              par <- object@par
              if (nc <= 1L)
                  list(vcomp = par, ccomp = double(0L))
              else list(vcomp = par[seq_len(length(par) - 1L)], ccomp = par[length(par)])
          })

setMethod("getVC",
          c(object = "Covariance.ar1"),
          .fn)

rm(.fn)

setMethod("getVCNames",
          c(object = "Covariance.us", cnm = "character", gnm = "character"),
          function (object, cnm, gnm, prf = NULL, old = TRUE) {
              nms <- getParNames(object, cnm, gnm, prf, old)
              nc <- object@nc
              ii <- cumsum(c(if (nc > 0L) 1L, if (nc > 1L) nc:2L))
              list(vcomp = nms[ii], ccomp = nms[-ii])
          })

setMethod("getVCNames",
          c(object = "Covariance.diag", cnm = "character", gnm = "character"),
          function (object, cnm, gnm, prf = NULL, old = TRUE) {
              nms <- getParNames(object, cnm, gnm, prf, old)
              list(vcomp = nms, ccomp = character(0L))
          })

setMethod("getVCNames",
          c(object = "Covariance.cs", cnm = "character", gnm = "character"),
          .fn <-
          function (object, cnm, gnm, prf = NULL, old = TRUE) {
              nms <- getParNames(object, cnm, gnm, prf, old)
              nc <- object@nc
              if (nc <= 1L)
                  list(vcomp = nms, ccomp = character(0L))
              else list(vcomp = nms[seq_len(nc)], ccomp = nms[nc + 1L])
          })

setMethod("getVCNames",
          c(object = "Covariance.ar1", cnm = "character", gnm = "character"),
          .fn)

rm(.fn)

setMethod("setVC",
          c(object = "Covariance.us", vcomp = "numeric", ccomp = "numeric"),
          function (object, vcomp, ccomp) {
              nc <- object@nc
              stopifnot(is.double(vcomp),
                        length(vcomp) == nc,
                        is.double(ccomp),
                        length(ccomp) == length(object@par) - nc)
              object@par <-
              if (nc <= 1L)
                  vcomp
              else {
                  i0 <- sequence.default(from = seq.int(from = nc + 1L, by = nc + 1L, length.out = nc - 1L),
                                         by = nc,
                                         nvec = (nc - 1L):1L)
                  i1 <- sequence.default(from = seq.int(from = 1L, by = nc + 1L, length.out = nc),
                                         by = nc,
                                         nvec = nc:1L)
                  S <- diag(1, nc, nc)
                  S[i0] <- ccomp
                  (safe_chol(S) * rep(vcomp, each = nc))[i1]
              }
              object
          })

setMethod("setVC",
          c(object = "Covariance.diag", vcomp = "numeric", ccomp = "numeric"),
          function (object, vcomp, ccomp) {
              stopifnot(is.double(vcomp),
                        length(vcomp) == length(object@par),
                        is.double(ccomp),
                        length(ccomp) == 0L)
              object@par <- vcomp
              object
          })

setMethod("setVC",
          c(object = "Covariance.cs", vcomp = "numeric", ccomp = "numeric"),
          .fn <-
          function (object, vcomp, ccomp) {
              nc <- object@nc
              stopifnot(is.double(vcomp),
                        length(vcomp) == length(object@par) - (nc > 1L),
                        is.double(ccomp),
                        length(ccomp) == (nc > 1L))
              object@par <- c(vcomp, ccomp)
              object
          })

setMethod("setVC",
          c(object = "Covariance.ar1", vcomp = "numeric", ccomp = "numeric"),
          .fn)

rm(.fn)

setMethod("getProfPar",
          c(object = "Covariance.us"),
          function (object, profscale, sc = NULL) {
              vc <- getVC(object)
              vcomp <- vc$vcomp
              ccomp <- vc$ccomp
              if (!is.null(sc))
                  vcomp <- vcomp * sc
              if (profscale == "varcov") {
                  nc <- object@nc
                  if (nc > 1L)
                      ccomp <- ccomp *
                          vcomp[sequence.default(from = 2L:nc, nvec = (nc - 1L):1L)] *
                          vcomp[rep(1L:(nc - 1L), (nc - 1L):1L)]
                  vcomp <- vcomp * vcomp
              }
              c(vcomp, ccomp)
          })

setMethod("getProfPar",
          c(object = "Covariance.diag"),
          function (object, profscale, sc = NULL) {
              vcomp <- getPar(object)
              if (!is.null(sc))
                  vcomp <- vcomp * sc
              if (profscale == "varcov")
                  vcomp <- vcomp * vcomp
              vcomp
          })

setMethod("setProfPar",
          c(object = "Covariance.us", value = "numeric"),
          function (object, profscale, sc = NULL, value, pos = 0L) {
              np <- length(object@par)
              validValuePos(np, value, pos)
              nc <- object@nc
              vcomp <- value[seq.int(from = pos + 1L     , length.out =      nc)]
              ccomp <- value[seq.int(from = pos + 1L + nc, length.out = np - nc)]
              if (profscale == "varcov") {
                  vcomp <- sqrt(vcomp)
                  if (nc > 1L)
                      ccomp <- ccomp/
                          vcomp[sequence.default(from = 2L:nc, nvec = (nc - 1L):1L)]/
                          vcomp[rep(1L:(nc - 1L), (nc - 1L):1L)]
              }
              if (!is.null(sc))
                  vcomp <- vcomp/sc
              setVC(object, vcomp, ccomp)
          })

setMethod("setProfPar",
          c(object = "Covariance.diag", value = "numeric"),
          function (object, profscale, sc = NULL, value, pos = 0L) {
              np <- length(object@par)
              validValuePos(np, value, pos)
              vcomp <-
              if (np == length(value))
                  value
              else
                  value[seq.int(from = pos + 1L, length.out = np)]
              if (profscale == "varcov")
                  vcomp <- sqrt(vcomp)
              if (!is.null(sc))
                  vcomp <- vcomp/sc
              setPar(object, vcomp)
          })

setMethod("getProfLower",
          c(object = "Covariance.us"),
          function (object, profscale, sc = NULL) {
              nc <- object@nc
              rep(c(0, if (profscale == "varcov") -Inf else -1),
                  c(nc, (nc * (nc - 1L)) %/% 2L))
          })

setMethod("getProfLower",
          c(object = "Covariance.diag"),
          function (object, profscale, sc = NULL) {
              nc <- object@nc
              rep(0, if (object@hom) 1L * (nc > 0L) else nc)
          })

setMethod("getProfUpper",
          c(object = "Covariance.us"),
          function (object, profscale, sc = NULL) {
              nc <- object@nc
              rep(c(Inf, if (profscale == "varcov") Inf else 1),
                  c(nc, (nc * (nc - 1L)) %/% 2L))
          })

setMethod("getProfUpper",
          c(object = "Covariance.diag"),
          function (object, profscale, sc = NULL) {
              nc <- object@nc
              rep(Inf, if (object@hom) 1L * (nc > 0L) else nc)
          })


## .... HELPERS ........................................................

## used in methods for set* which have formal arguments 'value', 'pos'
validValuePos <-
function (n, value, pos) {
    if (!is.double(value))
        stop(gettextf("type of '%s' is not \"%s\"",
                      "value", "double"),
             domain = NA)
    if (n > length(value) - pos)
        stop(gettextf("attempt to read past end of '%s'",
                      "value"),
             domain = NA)
    TRUE
}

## return a function that maps
##     c(theta_1, ..., theta_k) -> c(par_1, ..., par_k)
mkMkPar <-
function (reCovs) {
    if (all(vapply(reCovs, is, FALSE, "Covariance.us")))
        return(function (theta) theta)
    nt <- vapply(reCovs, getThetaLength, 0L, USE.NAMES = FALSE)
    np <- vapply(reCovs, getParLength, 0L, USE.NAMES = FALSE)
    snt <- sum(nt)
    snp <- sum(np)
    jt <- split(seq_len(snt), rep(seq_along(nt), nt))
    jp <- split(seq_len(snp), rep(seq_along(np), np))
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
    nt <- vapply(reCovs, getThetaLength, 0L, USE.NAMES = FALSE)
    np <- vapply(reCovs, getParLength, 0L, USE.NAMES = FALSE)
    snt <- sum(nt)
    snp <- sum(np)
    jt <- split(seq_len(snt), rep(seq_along(nt), nt))
    jp <- split(seq_len(snp), rep(seq_along(np), np))
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
## 'Lambdat', 'Lind', 'theta', 'par', 'lower', 'upper', 'reCovs' where
##
##     length(theta) >= length(par) == length(lower) == length(upper)
##
## note that previously there were no 'par', 'upper', and 'reCovs' as
## previously 'par' and 'theta' were identical by construction
upReTrms <-
function (reTrms, spCalls) {
    getHom <- function (spCall) {
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
    spNames <- as.character(lapply(spCalls, `[[`, 1L))
    nc <- lengths(reTrms$cnms, use.names = FALSE)
    hom <- vapply(spCalls, getHom, FALSE, USE.NAMES = FALSE)
    reCovs <- .mapply(new,
                      list(Class = paste0("Covariance.", spNames),
                           nc = nc,
                           hom = hom),
                      NULL)
    theta <- lapply(reCovs, getTheta)
    thetaIndex <- lapply(reCovs, getThetaIndex)
    nt <- lengths(theta)
    nti <- lengths(thetaIndex)
    nl <- reTrms$nl
    nc.nl <- rep(nc, nl)
    nti.nl <- rep(nti, nl)
    R.dp <- unlist(rep(lapply(reCovs, getLambdat.dp), nl), FALSE, FALSE)
    R.p <- cumsum(c(0L, R.dp))
    R.i <- rep(cumsum(c(0L, nc.nl)[seq_along(nc.nl)]), nti.nl) +
        unlist(rep(lapply(reCovs, getLambdat.i), nl), FALSE, FALSE)
    R.x <- rep(cumsum(c(0L, nt)[seq_along(nt)]), nti * nl) +
        unlist(rep(thetaIndex, nl), FALSE, FALSE)
    reTrms$Lambdat <- new("dgCMatrix", Dim = rep(length(R.dp), 2L),
                          p = R.p, i = R.i, x = as.double(R.x))
    reTrms$Lind    <- R.x
    reTrms$theta   <- unlist(                   theta, FALSE, FALSE)
    reTrms$par     <- unlist(lapply(reCovs, getPar  ), FALSE, FALSE)
    reTrms$lower   <- unlist(lapply(reCovs, getLower), FALSE, FALSE)
    reTrms$upper   <- unlist(lapply(reCovs, getUpper), FALSE, FALSE)
    reTrms$reCovs  <- reCovs
    reTrms
}

## use c(theta_1, ..., theta_k) to update par_1, ..., par_k
upReCovs <-
function (reCovs, theta) {
    ## FIXME?  function (reCovs, par) instead?
    pos <- 0L
    for (i in seq_along(reCovs)) {
        elt <- reCovs[[i]]
        reCovs[[i]] <- setTheta(elt, theta, pos)
        pos <- pos + getThetaLength(elt)
    }
    reCovs
}

## return 'reCovs' given 'object' of class "merMod", w/ or w/o attribute
getReCovs <-
function (object) {
    ## stopifnot(is(object, "merMod"))
    if (!is.null(ans <- attr(object, "reCovs")))
        return(ans)
    nc <- lengths(object@cnms, use.names = FALSE)
    upReCovs(lapply(nc, function (nc) new("Covariance.us", nc = nc)),
             object@theta)
}

convParToProfPar <-
function (value, object, profscale, sc = NULL) {
    reCovs <- getReCovs(object)
    np <- vapply(reCovs, getParLength, 0L, USE.NAMES = FALSE)
    pos <- cumsum(c(0L, np))[seq_along(np)]
    L <- .mapply(function (object, pos)
                     getProfPar(setPar(object, value, pos), profscale, sc),
                 list(object = reCovs, pos = pos),
                 NULL)
    c(L, sc, recursive = TRUE, use.names = FALSE)
}

convProfParToPar <-
function (value, object, profscale, sc = NULL) {
    reCovs <- getReCovs(object)
    np <- vapply(reCovs, getParLength, 0L, USE.NAMES = FALSE)
    pos <- cumsum(c(0L, np))[seq_along(np)]
    ## FIXME?  it could make sense to here do something like:
    ## if (is.null(sc) && length(value) > (snp <- sum(np)))
    ##     sc <- value[snp + 1L]
    L <- .mapply(function (object, pos)
                     getPar(setProfPar(object, profscale, sc, value, pos)),
                 list(object = reCovs, pos = pos),
                 NULL)
    c(L, recursive = TRUE, use.names = FALSE)
}


## .... METHODS [for class "merMod"] ...................................

setMethod("getPar",
          c(object = "merMod"),
          function (object)
              `length<-`(object@optinfo[["val"]],
                         length(object@lower)))

setMethod("getParLength",
          c(object = "merMod"),
          function (object)
              length(object@lower))

setMethod("getParNames",
          c(object = "merMod", cnm = "missing", gnm = "missing"),
          function (object, cnm, gnm, prf = NULL, old = TRUE) {
              reCovs <- getReCovs(object)
              gnm <- names(cnm <- object@cnms)
              unlist(.mapply(getParNames, list(object = reCovs, cnm = cnm, gnm = gnm), list(prf = prf, old = old)), FALSE, FALSE)
          })

setMethod("getTheta",
          c(object = "merMod"),
          function (object)
              object@theta)

setMethod("getThetaLength",
          c(object = "merMod"),
          function (object)
              length(object@theta))

setMethod("getThetaNames",
          c(object = "merMod", cnm = "missing", gnm = "missing"),
          function (object, cnm, gnm, prf = NULL, old = TRUE) {
              reCovs <- getReCovs(object)
              gnm <- names(cnm <- object@cnms)
              unlist(.mapply(getThetaNames, list(object = reCovs, cnm = cnm, gnm = gnm), list(prf = prf, old = old)), FALSE, FALSE)
          })

setMethod("getLower",
          c(object = "merMod"),
          function (object)
              object@lower)

setMethod("getUpper",
          c(object = "merMod"),
          function (object)
              attr(object, "upper") %||% rep(Inf, length(object@lower)))

setMethod("getVCNames",
          c(object = "merMod", cnm = "missing", gnm = "missing"),
          function (object, cnm, gnm, prf = NULL, old = TRUE) {
              reCovs <- getReCovs(object)
              gnm <- names(cnm <- object@cnms)
              L <- .mapply(getVCNames, list(object = reCovs, cnm = cnm, gnm = gnm), list(prf = prf, old = old))
              list(vcomp = lapply(L, `[[`, 1L),
                   ccomp = lapply(L, `[[`, 2L))
          })

setMethod("getProfPar",
          c(object = "merMod"),
          function (object, profscale, sc = NULL) {
              reCovs <- getReCovs(object)
              L <- lapply(reCovs, getProfPar, profscale = profscale, sc = sc)
              c(L, sc,
                recursive = TRUE, use.names = FALSE)
          })

setMethod("getProfLower",
          c(object = "merMod"),
          function (object, profscale, sc = NULL) {
              reCovs <- getReCovs(object)
              L <- lapply(reCovs, getProfLower, profscale = profscale)
              c(L, if (!is.null(sc)) 0,
                recursive = TRUE, use.names = FALSE)
          })

setMethod("getProfUpper",
          c(object = "merMod"),
          function (object, profscale, sc = NULL) {
              reCovs <- getReCovs(object)
              L <- lapply(reCovs, getProfUpper, profscale = profscale)
              c(L, if (!is.null(sc)) Inf,
                recursive = TRUE, use.names = FALSE)
          })

anyStructuredCov <- function(object) {
  re <- getReCovs(object)
  !all(vapply(re, function(x) inherits(x, "Covariance.us"),
              FUN.VALUE = logical(1)))
}
