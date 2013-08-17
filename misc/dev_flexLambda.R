if(FALSE){
#	library(lme4pureR)
#	library(nloptwrap)
	
	setwd("C:/lme4")
	library(devtools)
	library(ggplot2)
	library(mgcv)
	load_all()
	source('misc/reGenerators_flexLambda.R', echo=TRUE)
	options(error=recover)
}

(fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
(fm2 <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy))
(fm3 <- lmer(Reaction ~ Days + (0+ Days | Subject) + (1|Subject), sleepstudy))

setParams <- function(object,params,copy=TRUE) {
	if (copy) {
		newObj <- object
		newObj@pp <- newObj@pp$copy()
		newObj@resp <- newObj@resp$copy()
		if (!is.null(beta <- params$beta)) newObj@pp$setBeta0(beta)
		if (!is.null(theta <- params$theta)) {
			## where does theta live and how do I set it?
			## (1) in .@theta
			## (2) in .@pp$theta
			newObj@theta <- theta
			newObj@pp$setTheta(theta)
		}
		return(newObj)
	} else stop("modification in place (copy=FALSE) not yet implemented")
}

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
n <- 200
nsubj <- 20
data <- data.frame(subject = factor(gl(n/nsubj, nsubj)), 
				   treatment = factor(rep(rep(c(1,2), each=n/nsubj), times=n/nsubj)))
sd.treat <- c(.5, 2)
data <- within(data, {
	mu <- model.matrix(~0+treatment, data)%*%c(0, 5) + 
		model.matrix(~0+subject:treatment, data) %*% c(rnorm(n/nsubj, sd=sd.treat[1]), 
													   rnorm(n/nsubj, sd=sd.treat[2]))
	y <- mu + rnorm(n, sd=.1)
}) 
qplot(x=subject, y=y, col=treatment, group=subject, data=data)

lmod <- lFormula(y ~ treatment, data, reGenerators= ~d(~(0+treatment|subject)))			   
devfun <- do.call(mkLmerDevfun, lmod)
opt <- optimizeLmer(devfun)
environment(devfun)$pp$theta <- opt$par
(m <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr))
opt$par*sigma(m)


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
n <- 200
nsubj <- 20
data <- data.frame(subject = factor(gl(n/nsubj, nsubj)), 
				   treatment = factor(rep(rep(c(1,2), each=n/nsubj), times=n/nsubj)))
sd.treat <- 2
data <- within(data, {
	mu <- model.matrix(~0+treatment, data)%*%c(0, 5) + 
		model.matrix(~0+subject:treatment, data) %*% c(rnorm(n/nsubj, sd=sd.treat), 
													   rnorm(n/nsubj, sd=sd.treat))
	y <- mu + rnorm(n, sd=.1)
}) 
qplot(x=subject, y=y, col=treatment, group=subject, data=data)


lmod <- lFormula(y ~ treatment, data, reGenerators= ~d(~(0+treatment|subject), iid=TRUE))			   
devfun <- do.call(mkLmerDevfun, lmod)
opt <- optimizeLmer(devfun)
environment(devfun)$pp$theta <- opt$par
(m <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr))
opt$par*sigma(m)

devf <- pls(lmod, data$y)
opttheta <- bobyqa(lmod$reTrms$theta, devf, lower=lmod$reTrms$lower)

gm <- gam(y ~ treatment + s(subject, by=treatment, bs="re"), 
		  data=data, method="REML")

gm$coefficients[1:2]; environment(devf)$beta
opttheta$par; (x<-gam.vcomp(gm)[1:2,1])/sqrt(gm$reml.scale)
