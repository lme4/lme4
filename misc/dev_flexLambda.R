if(FALSE){
	library(lme4pureR)
	library(nloptwrap)
	
	setwd("C:/lme4")
	library(devtools)
	load_all()
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

devf <- pls(lmod,sleepstudy$Reaction)
opttheta <- bobyqa(c(1, 0, 1), devf, lower=c(0,-Inf,0))
environment(devf)$beta
opttheta$par

# opttheta$conv <- "foo"
# lmerdevf <- do.call(mkLmerDevfun, lmod)
# m <- mkMerMod(environment(lmerdevf), opttheta, 
# 			   lmod$reTrms, fr = lmod$fr)
# (setParams(m, list(beta=environment(devf)$beta, theta=opttheta$par)))


########################
#### diag/iid example

n <- 100
data <- data.frame(y=rnorm(n), x1=rnorm(n), x2=rnorm(n),
				   id=gl(20, n/20), id2=sample(gl(10, n/10)))
reGenerators=NULL
frml <- y ~ x1 + (x1|id) + (x1|id2)

fr <- model.frame(formula=as.formula(subbars(frml), env=data),
				  data=data, drop.unused.levels=TRUE)
attr(fr,"formula") <- frml
bars <- findbars(frml[[3]])

y ~ (1|id) + ar(~(.rows|1), order=1)

~ diag(~(f|g))  us(~(f|g))

n <- 20
rho <- .6
.rows <- 1:20
C <- rho^abs(outer(.rows, .rows, "-"))
image(L <- as(chol(C), "Matrix"))

formula <- ~ (0+f|id)

