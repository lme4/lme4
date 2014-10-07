## Laplace diagnostics

Copied/modified from DMB's code in glmer paper draft ...


```r
library(lme4)
library(lattice)
library(ggplot2); theme_set(theme_bw())
library(plyr)      ## for reshaping
library(abind)     ## ditto
library(reshape2)  ## for melt() generic
```

For now, I'm going to leave this function echoed in all its glory ...

```r
zetaDevfun <- function(m,grps=NULL) {    
    stopifnot (is(m, "glmerMod"),
               length(m@flist) == 1L)    # single grouping factor
    nvar <- length(m@cnms[[1]])
    ff <- getME(m,"flist")[[1]]
    if (is.null(grps)) grps <- seq(length(levels(ff)))  ## DRY ...
    ngrps <- length(grps)
    rr <- m@resp                         ## extract response module
    u0 <- getME(m,"u")                  ## conditional modes
    L <- getME(m,"L")
    ## sd <- 1/getME(pp,"L")@x
    ## filled elements of L matrix==diag for simple case
    ## for more general case need the following -- still efficient
    sd <- sqrt(diag(chol2inv(L)))
    ## fixed-effects contribution to linear predictor
    fc <- getME(m,"X") %*% getME(m,"beta") 
    ZL <- t(getME(m,"Lambdat") %*% getME(m,"Zt"))
    ## evaluate the unscaled conditional density on the deviance scale
    dc <- function(z) {
        uu <- u0 + z * sd    ## displace conditional modes
        ##  should still work if z is a vector (by recycling, because u values
        ##  applying to each group are stored adjacent to each other)
        rr$updateMu(fc + ZL %*% uu)     ## update linear predictor
        drc <- unname(as.vector(tapply(rr$devResid(), ff, sum)))
        uuc <- colSums(matrix(uu * uu,nrow=nvar))
        (drc + uuc)[grps]
    }
    return(dc)
}
zeta <- function(m, zmin=-3, zmax=3, npts=NULL,
                 grps=NULL, center=TRUE,
                 zvals = seq(zmin, zmax, length.out = npts)) {
    ff <- getME(m,"flist")[[1]]
    if (is.null(grps)) grps <- seq(length(levels(ff)))  ## DRY ...
    ngrps <- length(grps)
    nvar <- length(m@cnms[[1]])
    if (nvar>2) stop("can't handle vector RE with length >2")
    if (is.null(npts)) npts <- if (nvar>1) 31L else 301L
    dc <- zetaDevfun(m,grps)
    if (nvar==1) { # scalar-valued random effects
        vv <- vapply(zvals,dc,numeric(ngrps), USE.NAMES=FALSE)
        vv <- t(vv)  ## n.z * n.id
    } else { # vector-valued random effects
        nz <- length(zvals)
        vv <- mapply(function(x,y) { dc(c(x,y)) },
                         rep(zvals,nz),rep(zvals,each=nz))
        ## result: nu*(nz^2) matrix; want nz*nz*nu array
        ## *with* each nu slice being a nz^2 matrix for one group
        ## I'm sure there's a clever way to do this with array/aperm,
        ## but I just couldn't figure it out.  Instead,
        ## (1) take each row of vv and make into a matrix, return as list
        ##     of matrices
        ## (2) bind matrices into an array
        vv <- do.call(abind,c(alply(vv,1,matrix,nrow=nz),list(along=3)))
    }
    d0 <- dc(0) # restore the model to its incoming state
    devarr <- vv
    if (center) {
        sweep.margin <- if (nvar==1) 2 else 3 
        devarr <- sweep(devarr,sweep.margin,d0,"-")
    }
    ## computing deviance rather than signed sqrt, since we're not using it
    ## anyway and it's harder to generalize to >1 dimension ...
    rr <- list(zvals=zvals,
               devarr=devarr)
    ## signed square root
    ## array(ifelse(zvals < 0, -1, 1), c(npts, length(u0))))
    class(rr) <- "laplaceDiag"
    rr
}
```

Not shown: `melt` function for these objects that converts them from a list ($z$ value vector plus array of deviances) to a data frame ...


Replicate glmer paper Figs 2/3:


```r
m1 <- glmer(cbind(incidence, size-incidence) ~ period + (1|herd),
                  cbpp, binomial)
m1.z <- zeta(m1)
```


```r
plot(m1.z,layout=c(5,3))
```

![plot of chunk fig2](figure/fig2.png) 


```r
plot(m1.z,scaled=TRUE,layout=c(5,3))
```

![plot of chunk fig3](figure/fig3.png) 

### Example with vector-valued RE

I tried simulating a Poisson random-slopes model,
but at least at first glance it looked too good.
Try the toenail data 


```r
toenail <- read.csv("toenail.csv")
m2 <- glmer(outcome~treatment+visit+(visit|patient),toenail,
            family=binomial)
```


```r
m2.z <- zeta(m2,grps=1:25)
```

This is a nice collection of mussels ...

```r
plot(m2.z)
```

![plot of chunk plotmussels](figure/plotmussels.png) 

Still trying to work out the analogue of Figure 3 (i.e.,
a version where we scale by the bivariate normal).
I think this actually *should* work, but this example
is very badly behaved in the corner ...

Just look at patient #1 to try to sort out what's going on here ...

```r
zz <- m2.z$zvals
mm <- 2*pi*dnorm2d(sqrt(m2.z$devarr[,,1]))
m0 <- 2*pi*dnorm2d(sqrt(outer(zz^2,zz^2,"+")))
par(mfrow=c(2,2))
persp(zz,zz,mm,col="gray",main="conditional density")
persp(zz,zz,m0,col="lightblue",main="bivariate normal")
persp(zz,zz,mm/m0,col="pink",main="ratio")
persp(zz,zz,log(mm/m0),col="lightgreen",main="log ratio")
```

![plot of chunk ratios](figure/ratios.png) 

Does this really matter, or are we only in trouble if
we put quadrature points there?

### Further thoughts

How do we actually put GHQ into practice based on this information?
Can we come up with a sort of a score statistic for GHQ (i.e., what
would the difference be in log-likehood *conditional on parameters*
for different numbers of quadrature points)?

Check that `GHrule` and `GQdk` are equivalent (they're based on
the same underlying code, so they should be!)

```r
(gh5 <- GHrule(5))
```

```
##              z          w     ldnorm
## [1,] -2.856970 0.01125741 -5.0000774
## [2,] -1.355626 0.22207592 -1.8377997
## [3,]  0.000000 0.53333333 -0.9189385
## [4,]  1.355626 0.22207592 -1.8377997
## [5,]  2.856970 0.01125741 -5.0000774
```

```r
gh5B <- GQdk(1,5)
gh5B <- gh5B[order(gh5B[,2]),]
all.equal(gh5B[,2:1],unname(gh5[,-3]))
```

```
## [1] TRUE
```

Can we do some simple calculations with `GHrule` and the `zeta` values
to reconstruct the 1-D Gauss-Hermite results and convince ourselves we're
getting the same results as from `glmerAGQ` ?


```r
zd <- zetaDevfun(m1)
z2A <- zeta(m1,zvals=gh5[,"z"])
z2 <- t(sapply(gh5[,"z"],zd))
z2B <- sweep(z2,2,zd(0),"-")
all.equal(unname(z2A$devarr),z2B)
```

```
## [1] TRUE
```

```r
GH1d <- function(m,nAGQ=2) {
    gh <- GHrule(nAGQ)
    z <- zeta(m,zvals=gh[,"z"])
    sum(rowSums(z$devarr)*gh[,"w"])
}
GH1vals <- sapply(1:25,GH1d,m=m1)
plot(-GH1vals[-1])
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 


```r
m1GHvec <- sapply(2:25,
                  function(q) {
                      deviance(update(m1,nAGQ=q))
                  })
m1GHvec2 <- sapply(2:25,
                  function(q) {
                      dd <- update(m1,devFunOnly=TRUE,nAGQ=q)
                      dd(unlist(getME(m1,c("theta","fixef"))))
                  })
```


```r
par(las=1,bty="l")
matplot(cbind(diff(m1GHvec2),-diff(GH1vals[-1])),
        type="b",pch=16,lty=1)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 
These patterns don't match ...

`http://tigger.uic.edu/~hedeker/long.html`
