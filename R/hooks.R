.onLoad <- function(libname, pkgname) {
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
}

.onUnload <- function(libpath) {
  gc()
  if (is.loaded("lmer_Deviance", PACKAGE="lme4")) {
    library.dynam.unload("lme4", libpath)
  }
}
