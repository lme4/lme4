library(lme4)
library(reshape2)
library(plyr)
library(ggplot2)
theme_set(theme_bw())
revguide <- guide_legend(reverse = TRUE)
library(RColorBrewer)
data("Contraception", package = "mlmRev")
bootVars <- load("CBPP_bootbatch.rda")
profVars <- load("CBPP_profbatch.rda")

gm1 <- glmer(incidence/size ~ period + treatment + avg_size + (1 | herd),
             family = binomial,
             data = cbpp2, weights = size)

mfun <- function(x, type, model) {
  m <- reshape2::melt(rename(data.frame(var=rownames(x), x, check.names=FALSE),
                   c("2.5 %"="lwr", "97.5 %"="upr")), id="var")
  data.frame(model=model, type=type, m)
}

combCI0 <- do.call(rbind, c(
  mapply(mfun, list(confint0.boot2, confint0.wald),
         c("boot","Wald"), "gm1", SIMPLIFY=FALSE),
  mapply(mfun, list(confint.boot2, confint.wald),
         c("boot","Wald"), "gm2", SIMPLIFY=FALSE),
  mapply(mfun, list(confint3.prof, confint3.boot2, confint3.wald),
         c("profile","boot","Wald"), "gm3", SIMPLIFY=FALSE)
))

combCI0 <- combCI0[combCI0$var != ".sig01", ]
combCI0 <- combCI0[combCI0$var != "treatmentPartial/null", ]

get_est <- function(mod, model_name) {
  v <- as.data.frame(VarCorr(mod))
  vnames <- paste0("sd_", v$var1, "|", v$grp)
  est <- c(fixef(mod), setNames(v$sdcor, vnames))
  data.frame(model=model_name, var=names(est), est=est, row.names=NULL)
}

combEst2 <- do.call(rbind, list(
  get_est(gm1, "gm1"),
  get_est(gm2, "gm2"),
  get_est(gm3, "gm3")
))
cbpp2 <- transform(cbpp2,obs=factor(seq(nrow(cbpp2))))    
gm2 <- update(gm1,.~.+(1|obs))  ## herd and observation-level REs
gm3 <- update(gm1,.~.-(1|herd)+(1|obs))  ## observation-level REs only

combCI <- merge(dcast(combCI0, var + type + model ~ variable), combEst2, 
                by=c("var","model"), all.x=TRUE)
combCI <- transform(combCI, model_type = interaction(model, type, sep="\n"))
des_ord <- c("sd_(Intercept)|herd", "sd_(Intercept)|obs",
             "avg_size", "treatmentUnknown", "treatmentComplete")
combCI <- subset(combCI,var!="(Intercept)" & !grepl("^period",var))
combCI <- transform(combCI, 
                    var = factor(var, levels = rev(des_ord)))
mnames <- c("herd","herd+obs","obs")

avg_size_2sd <- 2*sd(cbpp2$avg_size)
scfun <- function(v) {
  combCI[combCI$var == "avg_size", v] <-
    combCI[combCI$var == "avg_size", v]*avg_size_2sd
  assign("combCI", combCI, parent.frame())
}
scfun("lwr")
scfun("upr")
scfun("est")

ggplot(combCI, aes(var, est, colour=type, shape=model)) +
  geom_pointrange(aes(ymin=lwr, ymax=upr),
                  position=position_dodge(width=0.6)) +
  labs(x="", y="Estimate (log-odds of seropositivity)") +
  geom_hline(yintercept=0, lty=2) +
  coord_flip() +
  scale_colour_brewer(palette="Dark2", guide = revguide) +
  scale_shape(guide=revguide, labels = mnames)
