#' Catch errors and warnings and store them for subsequent evaluation
#'
#' Factory modified from a version written by Martin Morgan on Stack Overflow (see below).  
#' Factory generates a function which is appropriately wrapped by error handlers.  
#' If there are no errors and no warnings, the result is provided.  
#' If there are warnings but no errors, the result is provided with a warn attribute set.
#' If there are errors, the result retutrns is a list with the elements of warn and err.
#' This is a nice way to recover from a problems that may have occurred during loop evaluation or during cluster usage.
#' Check the references for additional related functions.
#' I have not included the other factory functions included in the original Stack Overflow answer because they did not play well with the return item as an S4 object.
#' @export
#' @param fun The function to be turned into a factory
#' @param debug print debugging statements?
#' @param errval the value to be returned from the function if an error is thrown
#' @return The result of the function given to turn into a factory.  If this function was in error "An error as occurred" as a character element.  factory-error and factory-warning attributes may also be set as appropriate.
#' @references
#' \url{http://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function}
#' @author Martin Morgan; Modified by Russell S. Pierce and Ben Bolker
#' @examples 
#' f.log <- factory(log)
#' f.log("a")
#' f.log.NA <- factory(log,errval=NA)
#' f.log.NA("a")
#' f.as.numeric <- factory(as.numeric)
#' f.as.numeric(c("a","b",1))
factory <- function (fun, debug=FALSE, errval="An error occurred in the factory function") {
    function(...) {
    errorOccurred <- FALSE
    warn <- err <- msg <- NULL
    res <- withCallingHandlers(tryCatch(fun(...),
                                        error = function(e) {
        if (debug) cat("error: ",conditionMessage(e),"\n")
        err <<- conditionMessage(e)
        errorOccurred <<- TRUE
        NULL
    }), warning = function(w) {
      warn <<- append(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      if (debug) cat("message: ",conditionMessage(m),"\n")
      msg <<- append(msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
    if (errorOccurred) {
      res <- errval
    }

    if (is.character(msg)) {
        res <- setattr(res,"factory-message",msg)
    } else {
        res <- setattr(res,"factory-message",NULL) 
    }

    if (is.character(warn)) {
        res <- setattr(res,"factory-warning",warn)
    } else {
        res <- setattr(res,"factory-warning",NULL) 
    }

    if (is.character(err)) {
      res <- setattr(res,"factory-error",err)
    } else {
      res <- setattr(res, "factory-error", NULL)
    }  
    return(res)
  }
}

setattr <- function(x, attrib, value) {
    attr(x,attrib) <- value
    x
}
    
