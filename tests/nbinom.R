if (.Platform$OS.type != "windows") {
    library(lme4)
    cat("lme4 testing level: ", testLevel <- lme4:::testLevel(), "\n")


    getNBdisp <- function(x) getME(x,"glmer.nb.theta")
    ## for now, use hidden functions [MM: this is a sign, we should *export* them]
    refitNB   <- lme4:::refitNB

    simfun <- function(sd.u=1, NBtheta=0.5,
                       nblock = 25,
                       fform = ~x,
                       beta = c(1,2),
                       nrep = 40, seed) {
        levelset <- c(LETTERS,letters)
        stopifnot(2 <= nblock, nblock <= length(levelset))
        if (!missing(seed)) set.seed(seed)
        ntot <- nblock*nrep
        d1 <- data.frame(x = runif(ntot),
                         f = factor(rep(levelset[1:nblock], each=nrep)))
        u_f <- rnorm(nblock, sd=sd.u)
        X <- model.matrix(fform, data=d1)
        transform(d1, z = rnbinom(ntot, mu = exp(X %*% beta + u_f[f]), size = NBtheta))
    }

    ##' simplified logLik() so we can compare with "glmmADMB" (and other) results
    logLik.m <- function(x) {
        L <- logLik(x)
        attributes(L) <- attributes(L)[c("class","df","nobs")]
        L
    }

    if (testLevel > 1) withAutoprint({
        set.seed(102)
        d.1 <- simfun()
        t1 <- system.time(g1 <- glmer.nb(z ~ x + (1|f), data=d.1, verbose=TRUE))
        g1
        d1 <- getNBdisp(g1)
        (g1B <- refitNB(g1, theta = d1))
        (ddev <- deviance(g1) - deviance(g1B))
        (reld <- (fixef(g1) - fixef(g1B)) / fixef(g1))
        stopifnot(abs(ddev) < 1e-6,   # was 6.18e-7, 1.045e-6, -6.367e-5, now 0
                  abs(reld) < 1e-6)# 0, then 4.63e-6,  now 0
        ## 2 Aug 2015: ddev==reld==0 on 32-bit Ubuntu 12.04

        if(FALSE) {
            ## comment out to avoid R CMD check warning :
            ## library(glmmADMB)
            t2 <- system.time(g2 <- glmmadmb(z~x+(1|f),
                                             data = d.1, family="nbinom"))
            ## matrix not pos definite in sparse choleski
            t2 # 17.1 sec elapsed
            glmmADMB_vals <- list(fixef= fixef(g2),
                                  LL   = logLik(g2),
                                  theta= g2$alpha)
        } else {
            glmmADMB_vals <-
                list(fixef = c("(Intercept)" = 0.928710, x = 2.05072),
                     LL = structure(-2944.62, class = "logLik", df = 4, nobs = 1000L),
                     theta = 0.4487)
        }


        stopifnot(exprs = {
            all.equal(   d1,         glmmADMB_vals$ theta, tolerance=0.003) #   0.0015907
            all.equal(fixef(g1B),    glmmADMB_vals$ fixef, tolerance=0.02)# was 0.009387 !
            ## Ubuntu 12.04/32-bit: 0.0094
            all.equal(logLik.m(g1B), glmmADMB_vals$ LL,    tolerance=1e-4)# 1.681e-5; Ubuntu 12.04/32-b: 1.61e-5
        })

    })## end if( testLevel > 1 )

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
    data(epil, package="MASS")
    epil2 <- transform(epil,
                       Visit  = (period-2.5)/5,
                       Base   = log(base/4),
                       Age    = log(age),
                       subject= factor(subject))

    if(FALSE) {
        ## comment out to avoid R CMD check warning :
        ## library(glmmADMB)
        t3 <- system.time(g3  <- glmmadmb(y~Base*trt+Age+Visit+(Visit|subject),
                                          data=epil2, family="nbinom")) # t3 : 8.67 sec
        glmmADMB_epil_vals <- list(fixef= fixef(g3),
                                   LL   = logLik(g3),
                                   theta= g3$alpha)
    } else {
        glmmADMB_epil_vals <-
            list(fixef =
                     c("(Intercept)"= -1.33, "Base"=0.8839167, "trtprogabide"= -0.9299658,
                       "Age"= 0.4751434, "Visit"=-0.2701603, "Base:trtprogabide"=0.3372421),
                 LL = structure(-624.551, class = "logLik", df = 9, nobs = 236L),
                 theta = 7.4702)
    }

    if (testLevel > 2) withAutoprint({
        ## "too slow" for regular testing -- 49 (MM@lynne: 33, then 26, then 14) seconds:
        (t4 <- system.time(g4 <- glmer.nb(y ~ Base*trt + Age + Visit + (Visit|subject),
                                          data = epil2, verbose=TRUE)))
        ## 1.1-7 : Warning in checkConv().. failed .. with max|grad| = 0.0089 (tol = 0.001, comp. 4)
        ## 1.1-21: 2 Warnings:  max|grad| = 0.00859, then 0.1176 (0.002, comp. 1)

        stopifnot(exprs = {
            all.equal(getNBdisp(g4),   glmmADMB_epil_vals$ theta, tolerance= 0.03) # 0.0019777
            all.equal(fixef    (g4),   glmmADMB_epil_vals$ fixef, tolerance= 0.04) # 0.003731 (0.00374 on U 12.04)
            ## FIXME: even df differ (10 vs 9) !
            ##      all.equal(logLik.m(g4), - glmmADMB_epil_vals$ LL,	tolerance= 0.0) ## was 0.0002
            all.equal(logLik.m(g4), # for now {this is not *the* truth, just our current approximation of it}:
                      structure(-624.48418, class = "logLik", df = 10, nobs = 236L))
        })
    })


    cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
} ## skip on windows (for speed)
