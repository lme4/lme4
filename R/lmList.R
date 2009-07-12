## List of linear models according to a grouping factor

## Extract the model formula
modelFormula <- function(form)
{
    if (class(form) != "formula" || length(form) != 3)
        stop("formula must be a two-sided formula object")
    rhs <- form[[3]]
    if (class(rhs) != "call" || rhs[[1]] != as.symbol('|'))
        stop("rhs of formula must be a conditioning expression")
    form[[3]] <- rhs[[2]]
    list(model = form, groups = rhs[[3]])
}

setMethod("lmList", signature(formula = "formula", data = "data.frame"),
          function(formula, data, family, subset, weights,
                   na.action, offset, pool, ...)
      {
          mCall <- mf <- match.call()           
          m <- match(c("family", "data", "subset", "weights",
                       "na.action", "offset"), names(mf), 0)
          mf <- mf[c(1, m)]
          ## substitute `+' for `|' in the formula
          mf$formula <- subbars(formula) 
          mf$x <- mf$model <- mf$y <- mf$family <- NULL
          mf$drop.unused.levels <- TRUE
          mf[[1]] <- as.name("model.frame")
          frm <- eval(mf, parent.frame())
          mform <- modelFormula(formula)
          if (missing(family)) {
              val <- lapply(split(frm, eval(mform$groups, frm)),
                            function(dat, formula)
                        {
                            ans <- try({
                                data <- as.data.frame(dat)
                                lm(formula, data)
                            })
                            if (inherits(ans, "try-error"))
                                NULL
                            else ans
                        }, formula = mform$model)
          } else {
              val <- lapply(split(frm, eval(mform$groups, frm)),
                            function(dat, formula, family)
                        {
                            ans <- try({
                                data <- as.data.frame(dat)
                                glm(formula, family, data)
                            })
                            if (inherits(ans, "try-error"))
                                NULL
                            else ans
                        }, formula = mform$model, family = family)
          }
          if (missing(pool)) pool <- TRUE
          new("lmList", val, call = mCall, pool = pool)
      })


setMethod("coef", signature(object = "lmList"),
          ## Extract the coefficients and form a  data.frame if possible
          function(object, augFrame = FALSE, data = NULL,
                   which = NULL, FUN = mean, omitGroupingFactor = TRUE, ...)
      {
          coefs <- lapply(object, coef)
          non.null <- !unlist(lapply(coefs, is.null))
          if (sum(non.null) > 0) {
              template <- coefs[non.null][[1]]
              if (is.numeric(template)) {
                  co <- matrix(template,
                               ncol = length(template),
                               nrow = length(coefs),
                               byrow = TRUE,
                               dimnames = list(names(object), names(template)))
                  for (i in names(object)) {
                      co[i,] <- if (is.null(coefs[[i]])) { NA } else coefs[[i]]
                  }
                  coefs <- as.data.frame(co)
                  effectNames <- names(coefs)
                  if(augFrame) {
                      if (is.null(data)) {
                          data <- getData(object)
                      }
                      data <- as.data.frame(data)
                      if (is.null(which)) {
                          which <- 1:ncol(data)
                      }
                      data <- data[, which, drop = FALSE]
                      ## eliminating columns with same names as effects
                      data <- data[, is.na(match(names(data), effectNames)), drop = FALSE]
                      data <- gsummary(data, FUN = FUN, groups = getGroups(object))
                      if (omitGroupingFactor) {
                          data <- data[, is.na(match(names(data),
                                                     names(getGroupsFormula(object,
                                                                            asList = TRUE)))),
                                       drop = FALSE]
                      }
                      if (length(data) > 0) {
                          coefs <- cbind(coefs, data[row.names(coefs),,drop = FALSE])
                      }
                  }
                  attr(coefs, "level") <- attr(object, "level")
                  attr(coefs, "label") <- "Coefficients"
                  attr(coefs, "effectNames") <- effectNames
                  attr(coefs, "standardized") <- FALSE
                  #attr(coefs, "grpNames") <- deparse(getGroupsFormula(object)[[2]])
                  #class(coefs) <- c("coef.lmList", "ranef.lmList", class(coefs))
              }
          }
          coefs
      })

pooledSD <- function(x, ...)
{
    stopifnot(is(x, "lmList"))
    sumsqr <- apply(sapply(x,
                           function(el) {
                               if (is.null(el)) {
                                   c(0,0)
                               } else {
                                   res <- resid(el)
                                   c(sum(res^2), length(res) - length(coef(el)))
                               }
                           }), 1, sum)
    if (sumsqr[2] == 0) {
        stop("No degrees of freedom for estimating std. dev.")
    }
    val <- sqrt(sumsqr[1]/sumsqr[2])
    attr(val, "df") <- sumsqr[2]
    val
}

setMethod("show", signature(object = "lmList"),
          function(object)
      {
          mCall <- object@call
          cat("Call:", deparse(mCall), "\n")
          cat("Coefficients:\n")
          invisible(print(coef(object)))
          if (object@pool) {
              cat("\n")
              poolSD <- pooledSD(object)
              dfRes <- attr(poolSD, "df")
              RSE <- c(poolSD)
              cat("Degrees of freedom: ", length(unlist(lapply(object, fitted))),
                  " total; ", dfRes, " residual\n", sep = "")
              cat("Residual standard error:", format(RSE))
              cat("\n")
          }
      })

setMethod("confint", signature(object = "lmList"),
          function(object, parm, level = 0.95, ...)
      {
          mCall <- match.call()
          if (length(object) < 1)
              return(new("lmList.confint", array(numeric(0), c(0,0,0))))
          mCall$object <- object[[1]]
          template <- eval(mCall)
          val <- array(template, c(dim(template), length(object)),
                       c(dimnames(template), list(names(object))))
          pool <- list(...)$pool
          if (is.null(pool)) pool <- object$pool
          if (length(pool) > 0 && pool[1]) {
              sd <- pooledSD(object)
              a <- (1 - level)/2
              fac <- sd * qt(c(a, 1 - a)/2, attr(sd, "df"))
              parm <- dimnames(template)[[1]]
              for (i in seq_along(object))
                  val[ , , i] <-
                      coef(object[[i]])[parm] +
                          sqrt(diag(summary(object[[i]],
                                            corr = FALSE)$cov.unscaled
                                    )[parm]) %o% fac
          } else {
              for (i in seq_along(object)) {
                  mCall$object <- object[[i]]
                  val[ , , i] <- eval(mCall)
              }
          }
          new("lmList.confint", aperm(val, c(3, 2, 1)))
      }, valueClass = "lmList.confint")

setMethod("plot", signature(x = "lmList.confint"),
          function(x, y, ...)
      {
          stopifnot(require("lattice"))
          arr <- as(x, "array")
          dd <- dim(arr)
          dn <- dimnames(arr)
          levs <- dn[[1]]
          dots <- list(...)
          if (length(dots$order) > 0 &&
                     (ord <- round(dots$order[1])) %in% seq(dd[3]))
              levs <- levs[order(rowSums(arr[ , , ord]))]
          ll <- length(arr)
          df <- data.frame(group =
                           ordered(rep(dn[[1]], dd[2] * dd[3]),
                                   levels = levs),
                           intervals = as.vector(arr),
                           what = gl(dd[3], dd[1] * dd[2], len = ll, lab = dn[[3]]),
                           end = gl(dd[2], dd[1], len = ll))
          strip <- dots[["strip"]]
          if (is.null(strip)) {
              strip <- function(...) strip.default(..., style = 1)
          }
          xlab <- dots[["xlab"]]
          if (is.null(xlab)) xlab <- ""
          ylab <- dots[["ylab"]]
          if (is.null(ylab)) ylab <- ""
          dotplot(group ~ intervals | what,
                  data = df,
                  scales = list(x="free"),
                  strip = strip,
                  xlab = xlab, ylab = ylab,
                  panel = function(x, y, pch = dot.symbol$pch,
                  col = dot.symbol$col, cex = dot.symbol$cex,
                  font = dot.symbol$font, ...)
              {
                  x <- as.numeric(x)
                  y <- as.numeric(y)
                  ok <- !is.na(x) & !is.na(y)
                  yy <- y[ok]
                  xx <- x[ok]
                  dot.symbol <- trellis.par.get("dot.symbol")
                  dot.line <- trellis.par.get("dot.line")
                  panel.abline(h = yy, lwd = dot.line$lwd, lty = dot.line$lty, col =
                               dot.line$col)
                  lpoints(xx, yy, pch = "|", col = col, cex = cex, font = font, ...)
                  lower <- tapply(xx, yy, min)
                  upper <- tapply(xx, yy, max)
                  nams <- as.numeric(names(lower))
                  lsegments(lower, nams, upper, nams, col = 1, lty = 1, lwd =
                            if (dot.line$lwd) {
                                dot.line$lwd
                            } else {
                                2
                            })
              }, ...)
      })

setMethod("update", signature(object = "lmList"),
          function(object, formula., ..., evaluate = TRUE)
      {
          call <- object@call
          if (is.null(call))
              stop("need an object with call slot")
          extras <- match.call(expand.dots = FALSE)$...
          if (!missing(formula.))
              call$formula <- update.formula(formula(object), formula.)
          if (length(extras) > 0) {
              existing <- !is.na(match(names(extras), names(call)))
              for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
              if (any(!existing)) {
                  call <- c(as.list(call), extras[!existing])
                  call <- as.call(call)
              }
          }
          if (evaluate)
              eval(call, parent.frame())
          else call
      })

setMethod("formula", signature(x = "lmList"),
          function(x, ...) x@call[["formula"]])
