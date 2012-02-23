library(lme4Eigen)

simfun <- function(sd.u=1,NBtheta=0.5,
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

set.seed(102)
d1 <- simfun()
t1 <- system.time(g1 <- glmer.nb(z~x+(1|f),data=d1,debug=TRUE))

lme4Eigen:::getNBdisp(g1)
g1B <- lme4Eigen:::refitNB(g1,theta=lme4Eigen:::getNBdisp(g1))
deviance(g1)-deviance(g1B)
(fixef(g1)-fixef(g1B))/fixef(g1)

if (FALSE) {
  library(glmmADMB)
  t2 <- system.time(g2 <- glmmadmb(z~x+(1|f),
                                   data=d1,family="nbinom"))
  glmmADMB_vals <- list(fixef=fixef(g2),
                        NLL=-logLik(g2),
                        theta=g2$alpha)
  ## 0.4487
}

glmmADMB_vals <- structure(list(fixef = 
    structure(c(0.92871, 2.0507), .Names = c("(Intercept)", 
"x")), structure(2944.62, class = "logLik", df = 4, nobs = 1000L), 
    0.4487), .Names = c("fixef", "NLL", "theta"))

stopifnot(all.equal(glmmADMB_vals$theta,lme4Eigen:::getNBdisp(g1),
                     tol=0.0016))

if (FALSE) {
  ## simulation study
  simsumfun <- function(...) {
    d <- simfun(...)
    t1 <- system.time(g1 <- glmer.nb(z~x+(1|f),data=d))
    t2 <- system.time(g2 <- glmmadmb(z~x+(1|f),
                                     data=d,family="nbinom"))
    c(t.glmer=unname(t1["elapsed"]),nevals.glmer=g1$nevals,
      theta.glmer=exp(g1$minimum),
      t.glmmadmb=unname(t2["elapsed"]),theta.glmmadmb=g2$alpha)
  }

  library(plyr)
  sim50 <- raply(50,simsumfun(),.progress="text")
  save("sim50",file="nbinomsim1.RData")
  library(reshape)
  m1 <- melt(data.frame(run=seq(nrow(sim50)),sim50),id.var="run")
  m1 <- data.frame(m1,colsplit(m1$variable,"\\.",c("v","method")))
  m2 <- cast(subset(m1,v=="theta",select=c(run,value,method)),
             run~method)

  library(ggplot2)
  ggplot(subset(m1,v=="theta"),aes(x=method,y=value))+
    geom_boxplot()+geom_point()+geom_hline(yintercept=0.5,colour="red")

  ggplot(subset(m1,v=="theta"),aes(x=method,y=value))+
    stat_summary(fun.data=mean_cl_normal)+
      geom_hline(yintercept=0.5,colour="red")
  
  ggplot(m2,aes(x=glmer-glmmadmb))+geom_histogram()
  ## glmer is slightly more biased (but maybe the MLE itself is biased???)
}


### epilepsy example:
data(epil2,package="glmmADMB")
epil2$subject <- factor(epil2$subject)

if (FALSE) {
  t3 <- system.time(g3  <- glmmadmb(y~Base*trt+Age+Visit+(Visit|subject),
                                    data=epil2, family="nbinom"))
  glmmADMB_epil_vals <- list(fixef=fixef(g3),
                             NLL=-logLik(g3),
                             theta=g3$alpha)
}
glmmADMB_epil_vals <-
  structure(list(fixef = structure(c(-1.33, 0.88392, -0.92997, 
                   0.47514, -0.27016, 0.33724), .Names = c("(Intercept)", "Base", 
                                                  "trtprogabide", "Age", "Visit", "Base:trtprogabide")),
                 NLL = structure(624.551, class = "logLik", df = 9, nobs = 236L), 
                 theta = 7.4702), .Names = c("fixef", "NLL", "theta"))

if (FALSE) {
  ## 49 seconds: too slow for regular testing!
  t4 <- system.time(g4  <- glmer.nb(y~Base*trt+Age+Visit+(Visit|subject),
                                    data=epil2, debug=TRUE))
  stopifnot(all.equal(glmmADMB_epil_vals$theta,lme4Eigen:::getNBdisp(g4),
                      tol=0.0016))
}



