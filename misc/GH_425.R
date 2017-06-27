library(lme4)
library(optimx)
library(stringi)
library(data.table)
library(glmmTMB)

N <- 910000
## N should be a multiple of 26
N <- 26000

set.seed(1423L)
# highly imbalanced outcome variable
y <- sample.int(2L, size= N, replace=TRUE, prob= c(0.98, 0.02)) - 1L
# product biases
prod <- sample(letters, size= N, replace=TRUE)
#  user biases
my_grps <- stringi::stri_rand_strings(n= round(N/26), length= 10)
grps <- rep(my_grps, each= 26)
x1 <- sample.int(2L, size= N, replace=TRUE, prob= c(0.9, 0.1)) - 1L
x2 <- sample.int(2L, size= N, replace=TRUE, prob= c(0.9, 0.1)) - 1L
x3 <- sample.int(2L, size= N, replace=TRUE, prob= c(0.9, 0.1)) - 1L
x4 <- sample(LETTERS[1:5], size= N, replace=TRUE)

dt <- data.table(y= y,
             prod= prod, grps= grps,
             x1= x1, x2= x2, x3= x3, x4= x4)

## FIXME: may have to add an intercept to make the model match exactly?
t1 <- system.time({
    ## hack to make prod
    glmod <- glFormula(y ~ -1 + x1 + x2 + x3 + x4 + (1|prod) + (1 | grps),
                       data = dt, family=binomial)
    ## what order are the random effects in? (grps first)
    ## names(glmod$reTrms$Ztlist)
    ## 2.  Create the deviance function for optimizing over theta:
    devfun <- do.call(mkGlmerDevfun, glmod)
    ## check contents
    ## ls(environment(devfun))
    ## 3.  Optimize over theta using a rough approximation (i.e. nAGQ = 0):
    devfun2 <- function(theta) {
        return(devfun(c(theta,1000)))
    }
    m1 <- lme4:::optwrap("nloptwrap",fn=devfun2,par=1,lower=0,
                         control=list(),adj=FALSE,verbose=3,
                         calc.derivs=FALSE)
    ## 4.  Update the deviance function for optimizing over theta and beta:
    devfun <- updateGlmerDevfun(devfun, glmod$reTrms)
    ## 5.  Optimize over theta and beta:
    devfun3 <- function(par) {
        return(devfun(c(par[1],1000,par[-(1:2)])))
    }
    nbeta <- sum(environment(devfun)$lower==-Inf)
    ## turn off calc_derivs??
    m2 <- lme4:::optwrap("nloptwrap",fn=devfun3,
                         par=c(m1$par,1000,rep(0,nbeta)),
                         lower=c(0,rep(-Inf,nbeta)),
                         control=list(),adj=TRUE,verbose=3,
                         calc.derivs=FALSE)
    fit1 <- mkMerMod(environment(devfun), m2, glmod$reTrms, fr = glmod$fr)
})
save.image("GH_425.rda")
library(glmmTMB)

form <- y ~ -1 + prod + x1 + x2 + x3 + x4 + (1|grps)
t2 <- system.time(fit2 <- glmmTMB(form, data= dt, family= binomial))
save.image("GH_425.rda")
t3 <- system.time(fit3 <- glmer(form, data= dt, family= binomial,
                                nAGQ=0,
                                control=glmerControl(optimizer="nloptwrap")))
save.image("GH_425.rda")
t4 <- system.time(fit4 <- glmer(form, data= dt, family= binomial,
                                control=glmerControl(optimizer="nloptwrap",
                                                      calc.derivs=FALSE)))

t5 <- system.time(fit5 <- glmer(form, data= dt, family= binomial,
                                control=glmerControl(optimizer="nloptwrap")))

save.image("GH_425.rda")
t6 <- system.time(fit6 <- glmer(form,
                                data= dt, family= binomial,
                                control = glmerControl(optimizer ='optimx',
                               optCtrl=list(method='nlminb'),
                               calc.derivs=FALSE),verbose=1))

save.image("GH_425.rda")
