(load("Bachl_data.RData"))# d3
## originally gave 'd1' which was all numeric, then 'd2':
## d2 <- within(d1, {
##     idnr    <- as.factor(idnr)
##     turnid  <- as.factor(turnid)
##     kombiid <- as.factor(kombiid)
## })
## c(object.size(d2)/ object.size(d1))# 0.782
## str(d2)
## str(f.f <- with(d2, droplevels(idnr:turnid)))
## stopifnot(identical(as.integer(d2$kombiid), as.integer(f.f)))
## ## --> really do not need  'kombiid' in data
## ##-> d3 without kombiid:
## d3 <- d2; d3$kombiid <- NULL
## c(object.size(d3)/ object.size(d2)) ## 0.816

## note the y values are exactly from the integers -50:50 :
stopifnot(all(-50:50 == with(d3, sort(unique(rtr2)))))

## and construct the 2-way interaction factor, dropping non-existing levels:
d3[, "id.trn"] <- with(d3, droplevels(idnr:turnid))
str(d3)

library(lme4.0)
t0 <- system.time(m0 <- lmer(rtr2 ~ turnsec + (turnsec | id.trn) +
                             (turnsec | turnid) + (turnsec | idnr),
                             verbose = TRUE, d3))
outf <- "Bachl_out.RData"
save(list=ls(pattern="^[mt]"),file=outf)
detach("package:lme4.0")
library("lme4")
t1 <- system.time(m1 <- lmer(rtr2 ~ turnsec + (turnsec | id.trn) +
                             (turnsec | turnid) + (turnsec | idnr),
                             verbose = TRUE, d3))
save(list=ls(pattern="^[mt]"),file=outf)
t2 <- system.time(m2 <- update(m1,control=lmerControl(optimizer="bobyqa")))
save(list=ls(pattern="^[mt]"),file=outf)

library(nloptr)
defaultControl <- list(algorithm="NLOPT_LN_BOBYQA",
                       xtol_rel = 1e-6, maxeval = 1e5)
nloptwrap2 <- function(fn,par,lower,upper,control=list(),...) {
    for (n in names(defaultControl))
        if (is.null(control[[n]])) control[[n]] <- defaultControl[[n]]
    res <- nloptr(par, eval_f=fn, lb=lower, ub=upper, opts=control, ...)
    with(res,list(par=solution,
                  fval=objective,
                  conv=if (status>0) 0 else status,
                  message=message))
}

t3 <- system.time(m3 <- update(m1,control=lmerControl(optimizer="nloptwrap2")))
save(list=ls(pattern="^[mt]"),file=outf)

t4 <- system.time(m4 <- update(m1,control=lmerControl(optimizer="nloptwrap2",
			optCtrl=list(algorithm = "NLOPT_LN_NELDERMEAD"))))
save(list=ls(pattern="^[mt]"),file=outf)

###---------- Analysis: ------------------------
###           --------
## Start here, once you've run the above once { ESS:  M-[Enter] }
outf <- "Bachl_out.RData" # <- (repeated from above)
if(!(file.exists(outf) && "t4" %in% load(outf))  && !interactive())
    q("no")

source(system.file("testdata/lme-tst-funs.R", package="lme4", mustWork=TRUE))
##-> gimME(), allcoefs() ..

nms.m <- c("lme4.0",
           "i..NM",  "mq.byq", # "mq." : from 'minqa'  package
           "nl.byq", "nl.NM")  # "nl." : from 'nloptr'
tim <- rbind(t0, t1,t2,t3,t4)[, 1:3]
rownames(tim) <- nms.m

round(tim, 1)
##        user.self sys.self elapsed
## lme4.0     302.7      9.3   313.5  <-- lme4.0 clearly fastest
## i..NM      734.7      1.4   739.7
## mq.byq    1155.3      3.0  1164.2
## nl.byq     540.3     27.3   570.5
## nl.NM      725.0      3.0   731.6

mods <- setNames(lapply(ls(patt="^m[0-9]"), get), nms.m)
sapply(mods, object.size)
##    lme4.0    i..NM   mq.byq   nl.byq    nl.NM
## 34'267160  6471536  6471848  6471600  6472336

lapply(mods, getCall)# just to verify our  nms.m  are ok

## Optimizer convergence ok:
stopifnot(0 == sapply(lapply(mods[-1], slot, "optinfo"),
          function(.) .$conv$opt))
## lme4 convergence checks, all are *not* ok {but NM seem worse}
cbind(sapply(lapply(mods[-1], slot, "optinfo"),
             function(.) .$conv$lme4$message[1]))

smods <- lapply(mods, summary)
if(FALSE)## much output:
    print(smods)
lapply(smods, coef)

th.b <- t(sapply(mods, allcoefs))
## or rather also the std.err., t-value:
sb.b <- t(sapply(mods, allcoef.t))

## Look at the distances
D <- dist(th.b)
round(D, 2)
##        lme4.0 i..NM mq.byq nl.byq
## i..NM    1.95
## mq.byq   0.02  1.95
## nl.byq   0.02  1.95   0.00
## nl.NM    3.62  2.71   3.62   3.62

## or, similar
round(100* dist(sb.b))
##        lme4.0 i..NM mq.byq nl.byq
## i..NM     796
## mq.byq      2   797
## nl.byq      2   797      0
## nl.NM     920   405    920    920

## ===> Summary:
## lme4.0 gives practically *the same*  as the two bobyqa()s which
## are even closer to each other;
## The Nelder-Meads are quite different, and only partly close to each other

## Visualize via  MDS -- for the labels need large jitter:
cc <- cmdscale(D)
plot(cc)# see 3 points instead of 5
set.seed(7)
text(jitter(cc, f=400), rownames(cc), col=1:4, adj= c(.5, 1), xpd=NA)
