library(lme4)
cat("lme4 testing level: ", testLevel <- lme4:::testLevel(), "\n")


## for now, use hidden functions [MM: this is a sign, we should *export* them]
getNBdisp <- lme4:::getNBdisp
refitNB   <- lme4:::refitNB

simfun <- function(sd.u=1, NBtheta=0.5,
                   nblock=25,
                   fform=~x,
                   beta=c(1,2),
                   nrep=40,seed) {
    if (!missing(seed)) set.seed(seed)
    ntot <- nblock*nrep
    d1 <- data.frame(x=runif(ntot),f=rep(LETTERS[1:nblock],each=nrep))
    u_f <- rnorm(nblock,sd=sd.u)
    X <- model.matrix(fform,data=d1)
    transform(d1,z=rnbinom(ntot,
                 mu=exp(X %*% beta +u_f[f]),size=NBtheta))
}

if (testLevel > 1) {
    set.seed(102)
    d.1 <- simfun()
    t1 <- system.time(g1 <- glmer.nb(z ~ x + (1|f), data=d.1, verbose=TRUE))
    ## no longer: Model failed to converge with max|grad| = 0.00378242 (tol = 0.001)

    g1
    ## ^^ FIXME: 'Data:' shows '..2' ![eval.parent() etc.?]

    d1 <- getNBdisp(g1)
    (g1B <- refitNB(g1, theta = d1))
    (ddev <- deviance(g1) - deviance(g1B))
    (rel.d <- (fixef(g1) - fixef(g1B)) / fixef(g1))
    stopifnot(abs(ddev) < 2e-6,   # was 6.18e-7 now 1.045e-6
              abs(rel.d) < 0.0004)# now 0


## library(glmmADMB)
## t2 <- system.time(g2 <- glmmadmb(z~x+(1|f),
##                                  data=d1,family="nbinom"))
## glmmADMB_vals <- list(fixef=fixef(g2),
##                       NLL=-logLik(g2),
##                       theta=g2$alpha)
## 0.4487
glmmADMB_vals <-
    list(fixef = c("(Intercept)"=0.92871, x=2.0507),
         NLL = structure(2944.62, class = "logLik", df= 4, nobs= 1000L),
         theta = 0.4487)

##' simplified logLik() so we can compare with "glmmADMB" (and other) results
logLik.m <- function(x) {
    L <- logLik(x)
    attributes(L) <- attributes(L)[c("class","df","nobs")]
    L
}

stopifnot(
          all.equal(   d1,          glmmADMB_vals$theta, tolerance=0.0016)
          ,
          all.equal(fixef(g1B),     glmmADMB_vals$ fixef, tolerance=0.01)# not so close
          ,
          all.equal(logLik.m(g1B), -glmmADMB_vals$ NLL, tolerance=0.001)
          )
}## end if( testLevel > 1 )

if(FALSE) { ## simulation study --------------------

    ## library(glmmADMB) ## avoid R CMD check warning
    simsumfun <- function(...) {
        d <- simfun(...)
        t1 <- system.time(g1 <- glmer.nb(z~x+(1|f),data=d))
        t2 <- system.time(g2 <- glmmadmb(z~x+(1|f),
                                         data=d,family="nbinom"))
        c(t.glmer=unname(t1["elapsed"]),nevals.glmer=g1$nevals,
          theta.glmer=exp(g1$minimum),
          t.glmmadmb=unname(t2["elapsed"]),theta.glmmadmb=g2$alpha)
    }

    ## library(plyr)
    ## sim50 <- raply(50,simsumfun(),.progress="text")
    save("sim50",file="nbinomsim1.RData")
    ## library(reshape)
    ## m1 <- melt(data.frame(run=seq(nrow(sim50)),sim50),id.var="run")
    ## m1 <- data.frame(m1,colsplit(m1$variable,"\\.",c("v","method")))
    ## m2 <- cast(subset(m1,v=="theta",select=c(run,value,method)),
    ##           run~method)

    library(ggplot2)
    ggplot(subset(m1,v=="theta"),aes(x=method,y=value))+
        geom_boxplot()+geom_point()+geom_hline(yintercept=0.5,colour="red")

    ggplot(subset(m1,v=="theta"),aes(x=method,y=value))+
        stat_summary(fun.data=mean_cl_normal)+
            geom_hline(yintercept=0.5,colour="red")

    ggplot(m2,aes(x=glmer-glmmadmb))+geom_histogram()
    ## glmer is slightly more biased (but maybe the MLE itself is biased???)

}## end{simulation study}-------------------------

### epilepsy example:
data(epil,package="MASS")
epil2 <- transform(epil,Visit=(period-2.5)/5,
                   Base=log(base/4),Age=log(age),
                   subject=factor(subject))

## t3 <- system.time(g3  <- glmmadmb(y~Base*trt+Age+Visit+(Visit|subject),
##                                   data=epil2, family="nbinom"))
## glmmADMB_epil_vals <- list(fixef=fixef(g3),
##                            NLL=-logLik(g3),
##                            theta=g3$alpha)

glmmADMB_epil_vals <-
    list(fixef =
         c("(Intercept)"= -1.33, "Base"=0.88392, "trtprogabide"=-0.92997,
           "Age"=0.47514, "Visit"=-0.27016, "Base:trtprogabide"=0.33724),
         NLL = structure(624.551, class = "logLik", df = 9, nobs = 236L),
         theta = 7.4702)

if (testLevel > 3) {
    ## "too slow" for regular testing -- 49 (MM@lynne: 33, then 26) seconds:
    (t4 <- system.time(g4 <- glmer.nb(y~ Base*trt + Age + Visit + (Visit|subject),
                                      data=epil2, verbose=TRUE)))
    ## 1.1-7: Warning in checkConv().. failed .. with max|grad| = 0.0089 (tol = 0.001, comp. 4)

    (Lg4 <- logLik(g4))## logLik() --> ML instead of REML: refitting the model
    attributes(Lg4) <- attributes(Lg4)[c("class","df","nobs")]
    stopifnot(
              all.equal(getNBdisp(g4),   glmmADMB_epil_vals$ theta, tolerance= 0.0022)# was 0.002
              ,
              all.equal(fixef    (g4),   glmmADMB_epil_vals$ fixef, tolerance= 0.004)
              ,
              all.equal(logLik.m (g4), - glmmADMB_epil_vals$ NLL,	tolerance= 0.1) ## was 0.0002
              )
}

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
