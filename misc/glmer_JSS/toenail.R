data("toenail", package = "lme4")
toemod <- glmer(outcome ~ time*treatment + (1 | patientID), 
                     data = toenail, family = binomial(link = "logit"))

zeta <- function(m, zmin=-3, zmax=3, npts=301L, rank = FALSE) {
    stopifnot (is(m, "glmerMod"),
               length(m@flist) == 1L,    # single grouping factor
              length(m@cnms[[1]]) == 1L) # single column for that grouping factor
    pp <- m
    rr <- m@resp
    u0 <- getME(pp,"u")
    sd <- 1/getME(pp,"L")@x
    ff <- as.integer(getME(pp,"flist")[[1]])
    fc <- getME(pp,"X") %*% getME(pp,"beta") # fixed-effects contribution to linear predictor
    ZL <- t(getME(pp,"Lambdat") %*% getME(pp,"Zt"))
    dc <- function(z) { # evaluate the unscaled conditional density on the deviance scale
        uu <- u0 + z * sd
        rr$updateMu(fc + ZL %*% uu)
        unname(as.vector(tapply(rr$devResid(), ff, sum))) + uu * uu
    }
    zvals <- seq(zmin, zmax, length.out = npts)
    d0 <- dc(0) # because this is the last evaluation, the model is restored to its incoming state
    # signed square root
    sqrtmat <- t(sqrt(vapply(zvals, dc, d0, USE.NAMES=FALSE) - d0)) *
      array(ifelse(zvals < 0, -1, 1), c(npts, length(u0)))
    if (rank) {
      ## order by 'badness' (variance of normalized profile)
      nvals <- exp(-0.5*sqrtmat^2)/sqrt(2*pi)/dnorm(zvals)
      colnames(sqrtmat) <- seq(ncol(sqrtmat))
      sqrtmat <- sqrtmat[,order(apply(sqrtmat, 2, var), decreasing = TRUE)]
    }
    list(zvals=zvals,
         sqrtmat= sqrtmat 
         )
}

zm0 <- zeta(toemod, rank = FALSE)
zm <- zeta(toemod, rank = TRUE)


plot_zeta <- function(z, norm = FALSE,
                      ylab = if (!norm) "density" else "t(z)",
                      which = NULL, layout = c(5,3)) {
  dmat <- exp(-0.5*z$sqrtmat^2)/sqrt(2*pi)
  if (!is.null(which)) {
    dmat <- dmat[, which]
  }
  if (norm) dmat <- dmat/dnorm(z$zvals)
  xyplot(as.vector(dmat) ~ rep.int(z$zvals, ncol(dmat))|gl(ncol(dmat), nrow(dmat)),
         type=c("g","l"), aspect=0.6, layout=layout,
         xlab="z", ylab=ylab,
         panel=function(...){
           if (!norm) panel.lines(zm$zvals, dnorm(zm$zvals), lty=2)
           panel.xyplot(...)}
         )
}


plot_zeta(zm, which = 1:10, norm = TRUE)
set.seed(101); s <- sample(ncol(zm$sqrtmat), size = 9)
library(gridExtra)
grid.arrange(plot_zeta(zm, which = s, layout = c(3,3)),
             plot_zeta(zm, which = s, norm = TRUE, layout = c(3,3)), nrow = 1)



x11()
plot_zeta(zm0, n = 10, norm = TRUE)

plot_zeta(zm, norm = TRUE)
dmat <- exp(-0.5*zm$sqrtmat^2)/sqrt(2*pi)
dmat2 <- dmat/dnorm(zm$zvals)
vals <- apply(dmat2, 2, function(x) diff(range(x)))

zm_bad <- ...

worst <- which(vals > quantile(vals, 0.975) | vals < quantile(vals, )
library("lattice")
dmat_bad <- dmat[,worst]
## find worst examples?
xyplot(as.vector(dmat_bad) ~ rep.int(zm$zvals, ncol(dmat_bad))|gl(ncol(dmat_bad), nrow(dmat_bad)),
       type=c("g","l"), aspect=0.6, layout=c(5,3),
       xlab="z", ylab="density",
       panel=function(...){
           panel.lines(zm$zvals, dnorm(zm$zvals), lty=2)
           panel.xyplot(...)}
       )

xyplot(as.vector(dmat_bad/dnorm(zm$zvals)) ~ rep.int(zm$zvals, ncol(dmat_bad))|gl(ncol(dmat_bad), nrow(dmat_bad)),
       type=c("g","l"), aspect=0.6, layout=c(5,3),
       xlab="z", ylab="t(z)")
