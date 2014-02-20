# setParams <- function(object,params,copy=TRUE) {
# 	if (copy) {
# 		newObj <- object
# 		newObj@pp <- newObj@pp$copy()
# 		newObj@resp <- newObj@resp$copy()
# 		if (!is.null(beta <- params$beta)) newObj@pp$setBeta0(beta)
# 		if (!is.null(theta <- params$theta)) {
# 			## where does theta live and how do I set it?
# 			## (1) in .@theta
# 			## (2) in .@pp$theta
# 			newObj@theta <- theta
# 			newObj@pp$setTheta(theta)
# 		}
# 		return(newObj)
# 	} else stop("modification in place (copy=FALSE) not yet implemented")
# }

if(FALSE){
#	library(lme4pureR)
    ##	library(nloptwrap)
library(ggplot2)
	library(mgcv)
	library(mvtnorm)
	
##	setwd("C:/lme4")
##	library(devtools)
	load_all()
	source('misc/reGenerators_flexLambda.R', echo=TRUE)
	options(error=recover)
}

library(lme4)
library(MASS) ## for mvrnorm
## (rather than mvtnorm::rmvnorm -- try to stick to Recommended pkgs
library(mgcv)
source('reGenerators_flexLambda.R')
(fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
(fm2 <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy))
(fm3 <- lmer(Reaction ~ Days + (0+ Days | Subject) + (1|Subject), sleepstudy))


########################
### basic lmer example:
lmod <- lFormula(Reaction ~ Days + (Days|Subject), sleepstudy)
devfun <- do.call(mkLmerDevfun, lmod)
opt <- optimizeLmer(devfun)
environment(devfun)$pp$theta <- opt$par
mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)

# devf <- pls(lmod,sleepstudy$Reaction)
# opttheta <- bobyqa(c(1, 0, 1), devf, lower=c(0,-Inf,0))
# environment(devf)$beta
# opttheta$par


########################
## diag example
n <- 400
nsubj <- 40
data <- data.frame(subject = factor(gl(nsubj, n/nsubj)), 
				   treatment = factor(rep(rep(c(1,2), each=n/nsubj), times=n/nsubj)))
sd.treat <- c(1, 2)
data <- within(data, {
	mu <- model.matrix(~0+treatment, data)%*%c(0, 5) + 
		model.matrix(~0+subject:treatment, data) %*% c(rnorm(nsubj, sd=sd.treat[1]), 
													   rnorm(nsubj, sd=sd.treat[2]))
	y <- mu + rnorm(n, sd=.5)
}) 
qplot(x=subject, y=y, col=treatment, group=subject, data=data)

lmod <- lFormula(y ~ treatment, data, reGenerators= ~d(~(0+treatment|subject)))			   
devfun <- do.call(mkLmerDevfun, lmod)
opt <- optimizeLmer(devfun)
environment(devfun)$pp$theta <- opt$par
(m <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr))
opt$par*sigma(m); sigma(m)

qplot(x=subject, y=y, col=treatment, group=subject, data=data) + 
	geom_point(aes(y=fitted(m)), shape=2, size=2, col="black")

# devf <- pls(lmod, data$y)
# opttheta <- bobyqa(lmod$reTrms$theta, devf, lower=lmod$reTrms$lower)
# 
# data$tr1 <- 1*(data$treatment=="1")
# data$tr2 <- 1*(data$treatment=="2")
# gm <- gam(y ~ treatment + s(subject, by=tr1, bs="re") + s(subject, by=tr2, bs="re"), 
# 		  data=data, method="REML")
# 
# gm$coefficients[1:2]; environment(devf)$beta
# opttheta$par; (x<-gam.vcomp(gm)[1:2,1])/sqrt(gm$reml.scale)

########################
## iid example
n <- 400
nsubj <- 20
data <- data.frame(subject = factor(gl(nsubj, n/nsubj)), 
				   treatment = factor(rep(c(1,2), each=n/nsubj/2, times=n/nsubj)))
sd.treat <- 2
data <- within(data, {
	mu <- model.matrix(~0+treatment, data)%*%c(0, 5) + 
		model.matrix(~0+subject:treatment, data) %*% rnorm(2*nsubj, sd=sd.treat)
	y <- mu + rnorm(n, sd=.1)
}) 
qplot(x=subject, y=y, col=treatment, group=subject, data=data)


lmod <- lFormula(y ~ treatment, data, reGenerators= ~d(~(0+treatment|subject), iid=TRUE))			   
devfun <- do.call(mkLmerDevfun, lmod)
opt <- optimizeLmer(devfun)
environment(devfun)$pp$theta <- opt$par
(m <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr))
opt$par*sigma(m); sigma(m)

qplot(x=subject, y=y, col=treatment, group=subject, data=data) + 
	geom_point(aes(y=fitted(m)), shape=2, size=2, col="black")

#########################
## variance heterogeneity example
nsubj <- 50
nreps <- 10
nc <- 4
n <- nsubj*nreps*nc
data <- expand.grid(reps=factor(1:nreps),
					dose=1:nc,
					subject=factor(1:nsubj))

sd.dose <- rep(2, nc)# (1:nc)
cor <- .5
S <- {
	C <- diag(nc)
	C[upper.tri(C)] <- C[lower.tri(C)] <- cor
	sd.dose*t(sd.dose*C)
}

data <- within(data, {

    ## ## testing
    ## nsubj <- 5000
    ## reff0 <- as.vector(replicate(nsubj, mvtnorm::rmvnorm(1, sigma=S)[1,]))
    ## round(var(matrix(reff0,ncol=ncol(S),byrow=TRUE)),1)
    ## reff <- c(t(MASS::mvrnorm(nsubj, mu=rep(0,ncol(S)), Sigma=S)))
    ## round(var(matrix(reff,ncol=ncol(S),byrow=TRUE)),1)
    y <- dose + 
        model.matrix(~0+as.factor(dose):subject)%*%
            c(t(MASS::mvrnorm(nsubj, mu=rep(0,ncol(S)), Sigma=S)))+
               rnorm(n, sd=1)
    doseF <- as.factor(dose)
    subjectreps <- subject:reps
})
qplot(x=dose, y=y, group=subject, data=data, geom=c("point", "line")) + facet_wrap(~reps)
lmer(y ~ dose + (0+doseF|subject), data)

dlmod <- lFormula(y ~ dose, data, reGenerators= ~d(~(0+doseF|subject)))			   
ddevfun <- do.call(mkLmerDevfun, dlmod)
dopt <- optimizeLmer(ddevfun, verbose=5)
environment(ddevfun)$pp$theta <- dopt$par
(dm <- mkMerMod(environment(ddevfun), dopt, dlmod$reTrms, fr = dlmod$fr))
dopt$par*sigma(dm)

lmod <- lFormula(y ~ dose, data, reGenerators= ~cs(~(0+doseF|subject)))			   
devfun <- do.call(mkLmerDevfun, lmod)
opt <- optimizeLmer(devfun, verbose=5)
environment(devfun)$pp$theta <- opt$par
(m <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr))
estS <- with(environment(devfun),{
		Lt <- pp$Lambdat
		Lt@x <- pp$thfun(opt$par)
		crossprod(Lt[1:nc, 1:nc])
	})*sigma(m)^2 
sqrt(diag(estS)); cov2cor(estS); sigma(m)
plot(data$y, fitted(m))


hlmod <- lFormula(y ~ dose, data, reGenerators= ~cs(~(0+doseF|subject), het=FALSE))			   
hdevfun <- do.call(mkLmerDevfun, hlmod)
hopt <- optimizeLmer(hdevfun, verbose=5)
environment(hdevfun)$pp$theta <- hopt$par
(hm <- mkMerMod(environment(hdevfun), hopt, hlmod$reTrms, fr = hlmod$fr))
hestS <- with(environment(hdevfun),{
	Lt <- pp$Lambdat
	Lt@x <- pp$thfun(opt$par)
	crossprod(Lt[1:nc, 1:nc])
})*sigma(hm)^2 
sqrt(diag(hestS)); cov2cor(hestS); sigma(hm)
plot(data$y, fitted(m))


nc <- 5 
sd <- 2
rho <- .5
C <- diag(nc); C[upper.tri(C)] <- C[lower.tri(C)] <- rho
(S <- t(C*sd)*sd)
sd*sqrt(diagfactors(rho,1:nc))
sd*rho
(rho-rho^2)/sqrt(1-rho^2) * sd
(rho-rho^2)/sqrt(1-rho^2) * sd


chol(S)
#chol(t(S*sd)*sd)
theta <- c(sd, 0)

diagfactors <- function(x, col=nc) -(((col-1)*x^2- (col-2)*x - 1) /((col-2)*x + 1))
d0 <- which(sd < .Machine$double.eps)
aux <- chol(S[-d0,-d0])
cS <- diag(nc)
diag(cS)[d0] <- 0
cS[-d0, -d0] <- aux  
crossprod(cS)

tmp <- diag(nc)
sd <- theta[1:nc]
tmp[1, ] <- c(sd[1], rho*sd[1]*sd[-1])
diag(tmp)[-1] <- sqrt(diagfactors(rho, 2:nc))*sd[-1]
if(nc>2){
	tmp[2, -(1:2)] <- -(rho*(rho-1)*sd[2]*sd[-(1:2)])/sqrt(-rho*sd[2]^2+sd[2]^2)
	for(i in 3:nc){
		tmp[i, (i+1):nc] <- -(rho-1)*rho*sd[i] /	
	}	
}
tmp[2, ]
if(nc>2) {for(i in 2:nc){
	tmp[i, (i+1):nc] <- -	
}}

updateLambdatx <- local({
	# theta[1:nc]: sd's; theta[nc+1]: transformed corr
	tmp <- diag(nc)
	nl <- nl
	cl <- cl
	function(theta){
		#rho <-  diff(cl) * plogis(theta[nc+1]) + cl[1]
		
		
		
		C[upper.tri(C)] <- C[lower.tri(C)] <- rho
		#cS <- as.vector((chol(t(C*theta[1:nc])*theta[1:nc])))[upper.tri(C)]
		
		rep(cS, times=nl)
	}
})
