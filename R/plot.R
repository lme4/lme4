## copied/modified from nlme
splitFormula <-
  ## split, on the nm call, the rhs of a formula into a list of subformulas
  function(form, sep = "/")
{
  if (inherits(form, "formula") ||
      mode(form) == "call" && form[[1]] == as.name("~"))
    return(splitFormula(form[[length(form)]], sep = sep))
  if (mode(form) == "call" && form[[1]] == as.name(sep))
    return(do.call("c", lapply(as.list(form[-1]), splitFormula, sep = sep)))
  if (mode(form) == "(") return(splitFormula(form[[2]], sep = sep))
  if (length(form) < 1) return(NULL)
  list(asOneSidedFormula(form))
}

## Recursive version of all.vars
allVarsRec <- function(object)
{
  if (is.list(object)) {
    unlist(lapply(object, allVarsRec))
  } else {
    all.vars(object)
  }
}

## crippled version of getData.gnls from nlme
getData <-  function(object)
{
  mCall <- object@call
  data <- eval(mCall$data)
  if (is.null(data)) return(data)
  ## FIXME: deal with NAs, subset appropriately
  ## naPat <- eval(mCall$naPattern)
  ## if (!is.null(naPat)) {
  ##   data <- data[eval(naPat[[2]], data), , drop = FALSE]
  ## }
  ## naAct <- eval(mCall$na.action)
  ## if (!is.null(naAct)) {
  ##   data <- naAct(data)
  ## }
  ## subset <- mCall@subset
  ## if (!is.null(subset)) {
  ##   subset <- eval(asOneSidedFormula(subset)[[2]], data)
  ##   data <- data[subset, ]
  ## }
  return(data)
}

asOneFormula <-
  ## Constructs a linear formula with all the variables used in a
  ## list of formulas, except for the names in omit
  function(..., omit = c(".", "pi"))
{
    names <- unique(allVarsRec(list(...)))
    names <- names[is.na(match(names, omit))]
    if (length(names))
        as.formula(paste("~", paste(names, collapse = "+"))) # else NULL
}

## Return the formula(s) for the groups associated with object.
## The result is a one-sided formula unless asList is TRUE in which case
## it is a list of formulas, one for each level.
##
## @title
## @param object
## @param asList
## @param sep
## @return
getGroupsFormula <- function(object, asList = FALSE, sep = "+")
    UseMethod("getGroupsFormula")



getGroupsFormula.default <-
  ## Return the formula(s) for the groups associated with object.
  ## The result is a one-sided formula unless asList is TRUE in which case
  ## it is a list of formulas, one for each level.
  function(object, asList = FALSE, sep = "/")
{
  form <- formula(object)
  if (!inherits(form, "formula")){
    stop("\"Form\" argument must be a formula")
  }
  form <- form[[length(form)]]
  if (!((length(form) == 3) && (form[[1]] == as.name("|")))) {
    ## no conditioning expression
    return(NULL)
  }
  ## val <- list( asOneSidedFormula( form[[ 3 ]] ) )
  val <- splitFormula(asOneSidedFormula(form[[3]]), sep = sep)
  names(val) <- unlist(lapply(val, function(el) deparse(el[[2]])))
#  if (!missing(level)) {
#    if (length(level) == 1) {
#      return(val[[level]])
#    } else {
#      val <- val[level]
#    }
#  }
  if (asList) as.list(val)
  else as.formula(paste("~", paste(names(val), collapse = sep)))
}

getGroupsFormula.merMod <- function(object,asList=FALSE, sep="+") {
   if (asList) {
     lapply(names(object@flist),asOneSidedFormula)
   } else {
    asOneSidedFormula(paste(names(object@flist),collapse=sep))
   }
 }

getCovariateFormula <- function (object)
{
    form <- formula(object)
    if (!(inherits(form, "formula"))) {
        stop("formula(object) must return a formula")
    }
    form <- form[[length(form)]]
    if (length(form) == 3 && form[[1]] == as.name("|")) {
        form <- form[[2]]
    }
    eval(substitute(~form))
  }

getResponseFormula <-
  function(object)
{
  ## Return the response formula as a one sided formula
  form <- formula(object)
  if (!(inherits(form, "formula") && (length(form) == 3))) {
    stop("\"Form\" must be a two sided formula")
  }
  as.formula(paste("~", deparse(form[[2]])))
}

##--- needs Trellis/Lattice :
plot.merMod <-
  function(x, form = resid(., type = "pearson") ~ fitted(.), abline,
	   id = NULL, idLabels = NULL, idResType = c("pearson", "normalized"),
           grid, ...)
  ## Diagnostic plots based on residuals and/or fitted values
{
  object <- x
  if (!inherits(form, "formula"))
    stop("\"form\" must be a formula")
  ## constructing data
  ## can I get away with using object@frame???
  allV <- all.vars(asOneFormula(form, id, idLabels))
  allV <- allV[is.na(match(allV,c("T","F","TRUE","FALSE")))]
  if (length(allV) > 0) {
    data <- getData(object)
    if (is.null(data)) {		# try to construct data
       alist <- lapply(as.list(allV), as.name)
       names(alist) <- allV
       alist <- c(list(as.name("data.frame")), alist)
       mode(alist) <- "call"
       data <- eval(alist, sys.parent(1))
   } else if (any(naV <- is.na(match(allV, names(data)))))
         stop(allV[naV], " not found in data")
   } else data <- NULL

  ## this won't do because there may well be variables we want
  ##  that were not in the model call

  ## data <- object@frame

  ## argument list
  dots <- list(...)
  args <- if (length(dots) > 0) dots else list()
  ## appending object to data
  data <- as.list(c(as.list(data), . = list(object)))
  ## covariate - must always be present
  covF <- getCovariateFormula(form)
  .x <- eval(covF[[2]], data)
  if (!is.numeric(.x)) {
    stop("Covariate must be numeric")
  }
  argForm <- ~ .x
  argData <- data.frame(.x = .x, check.names = FALSE)
  if (is.null(args$xlab)) {
      if (is.null(xlab <- attr(.x, "label")))
	  xlab <- deparse(covF[[2]])
      args$xlab <- xlab
  }

  ## response - need not be present
  respF <- getResponseFormula(form)
  if (!is.null(respF)) {
    .y <- eval(respF[[2]], data)
    if (is.null(args$ylab)) {
	if (is.null(ylab <- attr(.y, "label")))
	    ylab <- deparse(respF[[2]])
	args$ylab <- ylab
    }
    argForm <- .y ~ .x
    argData[, ".y"] <- .y
  }

  ## groups - need not be present
  grpsF <- getGroupsFormula(form)
  if (!is.null(grpsF)) {
      ## ?? FIXME ???
      gr <- splitFormula(grpsF, sep = "*")
      for(i in seq_along(gr)) {
          auxGr <- all.vars(gr[[i]])
          for(j in auxGr)
              argData[[j]] <- eval(as.name(j), data)
      }
      argForm <-
	  as.formula(paste(if (length(argForm) == 2)
			   "~ .x |" else ".y ~ .x |",
			   deparse(grpsF[[2]])))
  }
  ## adding to args list
   args <- c(list(argForm, data = argData), args)
   if (is.null(args$strip)) {
      args$strip <- function(...) strip.default(..., style = 1)
    }
  if (is.null(args$cex)) args$cex <- par("cex")
  if (is.null(args$adj)) args$adj <- par("adj")

  if (!is.null(id)) {			# identify points in plot
    idResType <- match.arg(idResType)
    id <-
      switch(mode(id),
	     numeric = {
	       if (id <= 0 || id >= 1)
                   stop("Id must be between 0 and 1")
	       as.logical(abs(resid(object, type = idResType)) > -qnorm(id / 2))
	     },
	     call = eval(asOneSidedFormula(id)[[2]], data),
	     stop("\"Id\" can only be a formula or numeric.")
	     )
    if (is.null(idLabels)) {
      ## FIXME: Have no example yet, where this is triggered...
      message("**** getGroups() case in plot.merMod() ****")
      idLabels <- getGroupsFormula(object)
      if (length(idLabels) == 0) idLabels <- 1:object$dims$N
      idLabels <- as.character(idLabels)
    } else {
      if (mode(idLabels) == "call") {
	idLabels <-
	  as.character(eval(asOneSidedFormula(idLabels)[[2]], data))
      } else if (is.vector(idLabels)) {
	if (length(idLabels <- unlist(idLabels)) != length(id)) {
	  stop("\"IdLabels\" of incorrect length")
	}
	idLabels <- as.character(idLabels)
      } else {
	stop("\"IdLabels\" can only be a formula or a vector")
      }
    }
  }

  ## defining abline, if needed
  if (missing(abline)) {
      abline <- if (missing(form)) # r ~ f
          c(0, 0) else NULL
  }

  #assign("id", id , where = 1)
  #assign("idLabels", idLabels, where = 1)
  #assign("abl", abline, where = 1)
  assign("abl", abline)

  ## defining the type of plot
  if (length(argForm) == 3) {
    if (is.numeric(.y)) {		# xyplot
      plotFun <- "xyplot"
      if (is.null(args$panel)) {
        args <- c(args,
                  panel = list(function(x, y, subscripts, ...)
		    {
                      x <- as.numeric(x)
                      y <- as.numeric(y)
                      dots <- list(...)
		      if (grid) panel.grid()
		      panel.xyplot(x, y, ...)
                      if (any(ids <- id[subscripts])){
                          ltext(x[ids], y[ids], idLabels[subscripts][ids],
                                cex = dots$cex, adj = dots$adj)
                      }
		      if (!is.null(abl)) {
			if (length(abl) == 2) panel.abline(a = abl, ...) else panel.abline(h = abl, ...)
		      }
		    }))
      }
    } else {				# assume factor or character
      plotFun <- "bwplot"
      if (is.null(args$panel)) {
        args <- c(args,
                  panel = list(function(x, y, ...)
		    {
		      if (grid) panel.grid()
		      panel.bwplot(x, y, ...)
		      if (!is.null(abl)) {
			panel.abline(v = abl[1], ...)
		      }
		    }))
      }
    }
  } else {
    plotFun <- "histogram"
    if (is.null(args$panel)) {
      args <- c(args,
                panel = list(function(x, ...)
		  {
		    if (grid) panel.grid()
		    panel.histogram(x, ...)
		    if (!is.null(abl)) {
		      panel.abline(v = abl[1], ...)
		    }
		  }))
    }
  }

  ## defining grid
  if (missing(grid)) {
    grid <- (plotFun == "xyplot")
  }
  # assign("grid", grid, where = 1)
  do.call(plotFun, as.list(args))
}


fortify.merMod <- function(model, data=getData(model), ...) {
  ## FIXME: get influence measures via influence.ME?
  ##   (expensive, induces dependency ...)
  ## FIXME: different kinds of residuals?
  ## FIXME: deal with na.omit/predict etc.
  data$.fitted <- predict(model)
  data$.resid <- resid(model)
  data
}
## e.g. qplot(.fitted,.resid,data=fm1,colour=Sex)

## autoplot???

##  plot method for plot.summary.mer ... coefplot-style
##  horizontal, vertical? other options???
##  scale?
plot.summary.mer <- function(object, type="fixef", ...) {
  if(any(!type %in% c("fixef","vcov")))
      stop("'type' not yet implemented: ", type)
  stop("FIXME -- not yet implemented")

}

