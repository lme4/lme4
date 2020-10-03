## NOTE: Unlike the rest of the package, the functions in this file
## are licensed under the MIT License to encourage incorporation.
## 
## Copyright (c) 2020 Pavel N. Krivitsky and Benjamin Bolker
## 
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
## 
## The above copyright notice and this permission notice shall be included in all
## copies or substantial portions of the Software.
## 
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.
##

#' A `simulate` Method for `formula` objects that dispatches based on the Left-Hand Side
#'
#' This method evaluates the left-hand side (LHS) of the given formula and
#' dispatches it to an appropriate method based on the result by
#' setting an nonce class name on the formula.
#'
#' @param object a one- or two-sided [`formula`].
#' @param nsim,seed number of realisations to simulate and the random
#'   seed to use; see [simulate()].
#' @param ... additional arguments to methods.
#' @param basis if given, overrides the LHS of the formula for the
#'   purposes of dispatching.
#' @param newdata,data if passed, the `object`'s LHS is evaluated in
#'   this environment; at most one of the two may be passed.
#'
#' The dispatching works as follows:
#'
#' 1. If `basis` is not passed, and the formula has an LHS the
#'    expression on the LHS of the formula in the `object` is
#'    evaluated in the environment `newdata` or `data` (if given), in
#'    any case enclosed by the environment of `object`. Otherwise,
#'    `basis` is used.
#' 
#' 1. The result is set as an attribute `".Basis"` on `object`. If
#'    there is no `basis` or LHS, it is not set.
#' 
#' 1. The class vector of `object` has `c("formula_lhs_\var{CLASS}",
#'    "formula_lhs")` prepended to it, where \var{CLASS} is the class
#'    of the LHS value or `basis`. If LHS or `basis` has multiple
#'    classes, they are all prepended; if there is no LHS or `basis`,
#'    `c("formula_lhs_", "formula_lhs")` is.
#' 
#' 1. [simulate()] generic is evaluated on the new `object`, with all
#'    arguments passed on, excluding `basis`; if `newdata` or `data`
#'    are missing, they too are not passed on. The evaluation takes
#'    place in the parent's environment.
#'
#' A "method" to receive a formula whose LHS evaluates to \var{CLASS}
#' can therefore be implemented by a function
#' `simulate.formula_lhs_\var{CLASS}()`. This function can expect a
#' [`formula`] object, with additional attribute `.Basis` giving the
#' evaluated LHS (so that it does not need to be evaluated again).
#'
#' @export
## See https://github.com/lme4/lme4/issues/566 for further discussion
simulate.formula <- function(object, nsim=1, seed=NULL, ..., basis, newdata, data) {
  ## utility fun for generating new class
  cfun <- function(cc) {
    c(paste0("formula_lhs_", cc), "formula_lhs", class(object))
  }

  ## grab the arguments and the call and replace the function to be called with stats::simulate
  cl <- match.call()
  cl[[1L]] <- quote(stats::simulate)
  
  if (length(object)==3 || !missing(basis)) {  ## two-sided formula or basis given
    if (missing(basis)) { # If basis is not passed, evaluate LHS.
      if (!missing(data) && !missing(newdata)) stop("At most one of ", sQuote("data"), " or ", sQuote("newdata"), " can be specified.")

      evaldata <- if (!missing(data)) data
                  else if(!missing(newdata)) newdata
                  else environment(object)

      lhs <- object[[2L]]
      .Basis <-  try(eval(lhs, envir=evaldata,
                          enclos=environment(object)), silent=TRUE)

      if (inherits(.Basis,"try-error")) {
        ## can't evaluate LHS
        stop(simpleError(paste("Error evaluating the left-hand side of the formula:", .Basis)))
      }

    } else { # Otherwise, override LHS.
      .Basis <- basis
    }

    ## Set the basis object and class.
    attr(object,".Basis") <- .Basis
    class(object) <- cfun(class(.Basis))

  } else {  ## one-sided
    class(object) <- cfun("")
  }

  ## Replace the dispatched-on argument (object) with the updated formula.
  cl[["object"]] <- object
  ## If data argument was not actually passed, remove it from the call.
  if(missing(data)) cl <- cl[names(cl)!="data"]
  ## If newdata argument was not actually passed, remove it from the call.
  if(missing(newdata)) cl <- cl[names(cl)!="newdata"]
  ## Always remove basis from the call (since it's an attribute of object now).
  cl <- cl[names(cl)!="basis"]

  # Evaluate the modified call as if in the environment from which simulate.formula() was called. (A poor man's NextMethod().)
  eval(cl, parent.frame())
}

#' @describeIn simulate.formula A function to catch the situation when there is no method implemented for the class to which the LHS evaluates.
#' 
#' @export
simulate.formula_lhs <- function(object, nsim=1, seed=NULL, ...){
  stop("No applicable method for LHS of type ", paste0(sQuote(class(attr(object, ".Basis"))), collapse=", "), ".")
}
