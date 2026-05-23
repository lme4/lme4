library("lme4")

m1 <- glmer(incidence/size ~ period + (1|herd), weights = size, cbpp, family = binomial)
c1 <- confint(m1, method = "boot", nsim = 10)
c2 <- confint(m1, method = "boot", nsim = 10, parallel = "multicore")

cbpp2 <- read.csv("cbpp2.csv")
source("glmer_batch_funs.R")

cbpp2 <- transform(
  cbpp2,
  period = factor(period),
  treatment = factor(
    treatment,
    levels = c("Partial/null", "Complete", "Unknown")
  ),
  obs = factor(seq(nrow(cbpp2))),
  avg_size = avg_size / sd(avg_size, na.rm = TRUE)
)

cbpp_df_name <- data.frame(
  mnames = c("herd", "herdobs", "obs"),
  mdesc = c("herd", "herdobs", "obs"))

mforms <- list()
mforms$herd <- incidence / size ~ period + treatment + avg_size + (1 | herd)
mforms$herdobs <- update(mforms$herd, . ~ . + (1 | obs))
mforms$obs <- update(mforms$herd, . ~ . - (1|herd) + (1 | obs))

## don't use lapply because of parent.frame() eval nonsense -- check back once GH#961 is merged?
## this fails:
##   mod_list <- lapply(mforms, glmer, data = cbpp2, family = binomial, weights = size)
## when evaluating 'weights'

cc0 <- glmerControl(optimizer = "bobyqa",
                    ## avoid instability/problems with nAGQ=0 step ...
                    nAGQ0initStep = FALSE,
                    optCtrl = list(maxfun = 1e3),
                    calc.derivs = FALSE)
cc1 <- cc0
cc1$nAGQ0initStep <- TRUE

## sapply(cbpp_mod_list, \(x) -c(logLik(x)))
## sapply(cbpp_mod_list2, \(x) -c(logLik(x)))

cbpp_mod_list <- cbpp_mod_list2 <- list()
for (i in names(mforms)) {
  cbpp_mod_list[[i]] <- glmer(mforms[[i]],
                              data = cbpp2,
                              family = binomial,
                              weights = size,
                              control = cc0),
  cbpp_mod_list2[[i]] <- glmer(mforms[[i]],
                              data = cbpp2,
                              family = binomial,
                              weights = size,
                              control = cc1)

}
