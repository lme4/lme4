library(tinyplot)
## https://grantmcdermott.com/tinyplot/vignettes/introduction.html
## FIXME: allow abbreviation of optimizer names?
plot.allFit <- function(object, which = "fixed",
                        conf.level = 0.95,
                        facet.args = list(free = TRUE),
                        ...) {
  
  do_ci <- !is.null(conf.level) && !is.na(conf.level)
  get_fe <- function(x, opt) {
    cc <- coef(summary(x))
    dd <- data.frame(optimizer = opt, term = rownames(cc), estimate = cc[,"Estimate"], sd = cc[,"Std. Error"])
    if (do_ci) {
      qq <- qnorm((1+conf.level)/2)
      dd <- transform(dd,
                      lwr = estimate - qq*sd,
                      upr = estimate + qq*sd)
    }
    rownames(dd) <- NULL
    dd
  }
  mm <- Map(get_fe, object, names(object))
  dd <- do.call(rbind, mm)
  tinyplot_args <- list(
    x = estimate ~ optimizer,
    facet = ~ term,
    facet.args = facet.args,
    data = dd,
    ...)
  if (do_ci) {
    tinyplot_args <- c(tinyplot_args,
                       list(ymin = dd$lwr,
                            ymax = dd$upr,
                            type = "pointrange"))
  } else tinyplot_args <- c(tinyplot_args, list(type = "p"))
  do.call(tinyplot, tinyplot_args)
}

plot(gm_all)
par(las=1); tpar(fmar = c(1,12,1,1))
plot(gm_all, flip = TRUE)
plot(gm_all, conf.level = NA)
## not working yet ...
## plot(gm_all, conf.level = NA, by = "optimizer")
