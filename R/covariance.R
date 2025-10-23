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
         validity = function (object) {
             if (length(nc <- object@nc) != 1L || is.na(nc) || nc < 0L)
                 "'nc' is not a non-negative integer"
             else if (!is.double(object@par))
                 "type of 'par' is not \"double\""
             else TRUE
         })

setClass("UnstructuredCovariance",
         contains = "Covariance",
         validity = function (object) {
             nc <- object@nc
             if (length(object@par) != (nc * (nc - 1L)) %/% 2L + nc)
                 "length of 'par' is not nc*(nc+1)/2"
             else TRUE
         })

setClass("DiagonalCovariance",
         contains = "Covariance",
         validity = function (object) {
             nc <- object@nc
             if (length(par <- object@par) != 1L && length(par) != nc)
                 "length of 'par' is not 1 and not 'nc'"
             else TRUE
         })


## .... GENERIC FUNCTIONS ..............................................

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
function (object, ...) {
    ## stopifnot(is(object, "Covariance"))
    object@par
}

setGeneric("setPar",
           function (object, value, ...)
               standardGeneric("setPar"))

setGeneric("getTheta",
           function (object, ...)
               standardGeneric("getTheta"))

setGeneric("setTheta",
           function (object, value, ...)
               standardGeneric("setTheta"))


## .... METHODS ........................................................

## For class "UnstructuredCovariance"

setMethod("initialize",
          c(.Object = "UnstructuredCovariance"),
          function (.Object, nc, ...) {
              if (nargs() == 2L && !missing(nc) && is.integer(nc) &&
                  length(nc) == 1L && !is.na(nc) && nc >= 0L) {
                  .Object@par <- `[<-`(double((nc * (nc - 1L)) %/% 2L + nc),
                                       cumsum(c(if (nc > 0L) 1L, if (nc > 1L) nc:2L)),
                                       1)
                  .Object
              }
              else callNextMethod()
          })

setMethod("setPar",
          c(object = "UnstructuredCovariance", value = "numeric"),
          function (object, value, ...) {
              nc <- object@nc
              if (length(value) != (nc * (nc - 1L)) %/% 2L + nc)
                  stop("length of 'value' is not nc*(nc+1)/2")
              if (!is.double(value))
                  storage.mode(value) <- "double"
              object@par <- value
              object
          })

setMethod("getTheta",
          c(object = "UnstructuredCovariance"),
          function (object, ...)
              object@par)

setMethod("setTheta",
          c(object = "UnstructuredCovariance", value = "numeric"),
          function (object, value, ...) {
              nc <- object@nc
              if (length(value) != (nc * (nc - 1L)) %/% 2L + nc)
                  stop("length of 'value' is not nc*(nc+1)/2")
              if (!is.double(value))
                  storage.mode(value) <- "double"
              object@par <- value
              object
          })


## For class "DiagonalCovariance"

setMethod("initialize",
          c(.Object = "DiagonalCovariance"),
          function (.Object, nc, ...) {
              if (nargs() == 2L && !missing(nc) && is.integer(nc) &&
                  length(nc) == 1L && !is.na(nc) && nc >= 0L) {
                  .Object@par <- rep(1, nc)
                  .Object
              }
              else callNextMethod()
          })

setMethod("setPar",
          c(object = "DiagonalCovariance", value = "numeric")
          function (object, value, ...) {
              nc <- object@nc
              if (length(value) != 1L && length(value) != nc)
                  stop("length of 'value' is not 1 and not 'nc'")
              if (!is.double(value))
                  storage.mode(value) <- "double"
              object@par <- value
              object
          })

setMethod("getTheta",
          c(object = "DiagonalCovariance"),
          function (object, ...) {
              nc <- object@nc
              `[<-`(double((nc * (nc - 1L)) %/% 2L + nc),
                    cumsum(c(if (nc > 0L) 1L, if (nc > 1L) nc:2L)),
                    object@par)
          })

setMethod("setTheta",
          c(object = "DiagonalCovariance", value = "numeric"),
          function (object, value, ...) {
              nc <- object@nc
              if (length(value) != (nc * (nc - 1L)) %/% 2L + nc)
                  stop("length of 'value' is not nc*(nc+1)/2")
              if (!is.double(value))
                  storage.mode(value) <- "double"
              i <-
              if (length(object@par) == 1L)
                  1L
              else cumsum(c(if (nc > 0L) 1L, if (nc > 1L) nc:2L))
              object@par <- value[i]
              object
          })
