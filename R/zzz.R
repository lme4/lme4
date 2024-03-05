.onLoad <- function(libname, pkgname) {
  ## don't do this in production; also flags problems in downstream packages
  ## options(Matrix.warnDeprecatedCoerce = 3)
  options(lme4.summary.cor.max = 12)
  if((Rv <- getRversion()) < "4.1.0") {
    ## https://stackoverflow.com/questions/49056642/how-to-make-variable-available-to-namespace-at-loading-time/67664852#67664852
    ## not quite equivalent; this *forces* ... entries whereas true ...length()  doesn't
    assign('...names', envir = topenv(),
           function() eval(quote(names(list(...))), sys.frame(-1L)))
    if(Rv < "4.0.0") {
      ## NB: R >= 4.0.0's deparse1() is a generalization of our previous safeDeparse()
      assign('deparse1', envir = topenv(),
             function (expr, collapse = " ", width.cutoff = 500L, ...)
               paste(deparse(expr, width.cutoff, ...), collapse = collapse)
             )
      ## not equivalent ...
      assign('...length', envir = topenv(),
             function() eval(quote(length(list(...))), sys.frame(-1L))
             )
      if (Rv < "3.6.0") {
        assign('reformulate', envir = topenv(),
               function(..., env = parent.env) {
                   f <- stats::reformulate(...)
                   environment(f) <- env
                   return(f)
               })
        if (Rv < "3.2.1") {
        assign('lengths', envir = topenv(),
               function (x, use.names = TRUE) vapply(x, length, 1L, USE.NAMES = use.names)
               )
        if(Rv < "3.1.0") {
          assign('anyNA', envir = topenv(),
                 function(x) any(is.na(x))
                 )
          if(Rv < "3.0.0") {
            assign('rep_len', envir = topenv(),
                   function(x, length.out) rep(x, length.out=length.out)
                   )
            if(Rv < "2.15") {
              assign('paste0', envir = topenv(),
                     function(...) paste(..., sep = '')
                     )
            } ## R < 2.15
          } ## R < 3.0.0
        } ## R < 3.1.0
        } ## R < 3.2.1
      } ## R < 3.6.0
    } ## R < 4.0.0
  } ## R < 4.1.0
  rm(Rv)
  ## check Matrix ABI version
  check_dep_version()  
}

## https://github.com/lme4/lme4/issues/768
## https://github.com/kaskr/adcomp/issues/387
get_abi_version <- function() {
    if (utils::packageVersion("Matrix") < "1.6-2") return(numeric_version("0"))
    Matrix::Matrix.Version()[["abi"]]
}

.Matrix.abi.build.version <- get_abi_version()

## simplified version of glmmTMB package checking
##' @param this_pkg downstream package being tested
##' @param dep_pkg upstream package on which \code{this_pkg} depends
##' @param dep_type "ABI" or "package"
##' @param built_version a \code{numeric_version} object indicating what version of \code{dep_pkg} was used to  build \code{this_pkg}
##' @param warn (logical) warn if condition not met?
##' @noRd
check_dep_version <- function(this_pkg = "lme4",  dep_pkg = "Matrix", dep_type = "ABI",
                              built_version = .Matrix.abi.build.version,
                              warn = TRUE) {
    cur_version <- get_abi_version()
    result_ok <- cur_version == built_version
    if (!result_ok) {
        warning(
            sprintf("%s version mismatch: \n", dep_type),
            sprintf("%s was built with %s %s version %s\n",
                    this_pkg, dep_pkg, dep_type, built_version),
            sprintf("Current %s %s version is %s\n",
                    dep_pkg, dep_type, cur_version),
            sprintf("Please re-install %s from source ", this_pkg),
            "or restore original ",
            sQuote(dep_pkg), " package"
        )
    }
    return(result_ok)
}

.onUnload <- function(libpath) {
  gc()
  if (is.loaded("lmer_Deviance", PACKAGE="lme4")) {
    library.dynam.unload("lme4", libpath)
  }
}
