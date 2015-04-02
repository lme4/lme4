library("lme4")
library("numDeriv")
library("gridExtra")

do_sim <- FALSE

simfun <- function(size=nrow(Penicillin),seed=NULL) {
    if (!is.null(seed)) set.seed(seed)
    sdat <- Penicillin[sample(x=nrow(Penicillin),size=size,replace=TRUE),]
    fm2 <- lmer(diameter ~ 1 + (1|plate) + (1|sample), sdat)
    derivs <- fm2@optinfo$derivs
    ## copied from lme4:::checkConv
    mineig <- min(eigen(fm2@optinfo$derivs$Hessian,only.values=TRUE)$value)
    dd <- as.function(fm2)
    hh <- hessian(dd,getME(fm2,"theta"))
    mineigND <- min(eigen(hh,only.values=TRUE)$value)
    scgrad <- tryCatch(with(derivs, solve(Hessian, gradient)), 
                       error = function(e) e)
    if (inherits(scgrad, "error")) {
        wstr <- "unable to evaluate scaled gradient"
        maxscgrad <- maxmingrad <- NA
    } else {
        maxscgrad <- max(abs(scgrad))
        mingrad <- pmin(abs(scgrad), abs(derivs$gradient))
        maxmingrad <- max(mingrad)
    }
    return(c(maxgrad=max(abs(derivs$gradient)),
             maxscgrad=maxscgrad,
             maxmingrad=maxmingrad,
             mineig=mineig,
             mineigND=mineigND,
             getME(fm2,"theta")))
}


if (do_sim) {
    res <- expand.grid(log10size=seq(2,6,by=0.5),rep=1:20,
                       maxgrad=NA,maxscgrad=NA,maxmingrad=NA,
                       mineig=NA,
                       mineigND=NA,
                       theta1=NA,theta2=NA)
    for (i in 1:nrow(res)) {
        ## want to save 
        cat(i,res$log10size[i],"\n")
        res[i,-(1:2)] <- simfun(size=round(10^res$log10size[i]),seed=1000+i)
        save("res",file="penicillin_conv.RData")
    }
}

## PLOTS

library("ggplot2"); theme_set(theme_bw())
## with apologies for Hadleyverse 2 ...
library("tidyr")
library("dplyr")
load("penicillin_conv.RData")
m <- res %>%
    mutate_each("log10",c(maxgrad,maxscgrad,maxmingrad)) %>%
        gather(var,val,-c(rep,log10size))

ggplot(m,aes(log10size,val))+geom_point()+
    facet_wrap(~var,scales="free")+
        geom_smooth(method="lm",formula=y~1+offset(x))

## PLOT: gradients vs log10size
(gg1 <- ggplot(filter(m,var %in% c("maxgrad","maxscgrad","maxmingrad")),
               aes(log10size,val,colour=var))+geom_point()+
                   geom_smooth(method="lm",formula=y~1+offset(x))+
                       geom_smooth(method="loess",linetype=2)+
                           geom_hline(yintercept=-3,lty=2))


(gg2 <- ggplot(filter(m,var %in% c("mineig","mineigND")),
               aes(log10size,val,colour=var))+geom_point(alpha=0.25,size=5))

library("scales")
(gg3 <- gg2 + scale_y_continuous(limits=c(0,20),oob=squish))
(gg4 <- gg2 + scale_y_continuous(limits=c(0.5,2),oob=squish))

pdf("penicillin_conv.pdf")
grid.arrange(gg1,gg2,gg3,gg4,nrow=2)
dev.off()

filter(res,maxmingrad>1e-3)

