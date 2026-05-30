library("lme4")
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
                    optCtrl = list(maxfun = 1e5))
cc1 <- cc0
cc1$calc.derivs <- FALSE

cbpp_mod_list <- list()
for (i in names(mforms)) {
  cbpp_mod_list[[i]] <- glmer(mforms[[i]],
                              data = cbpp2,
                              family = binomial,
                              weights = size,
                              control = cc0)
}

## turn off derivs for bootstrap CIs
## (need derivs for hacked Wald CIs on REs)
cbpp_mod_list_noderivs <- lapply(cbpp_mod_list,
                                 function(x) update(x, control = cc1))

## b_cifun(cbpp_mod_list2[[2]], nsim = 5)
cbpp_confint_boot <- lapply(cbpp_mod_list_noderivs,
                            b_cifun, nsim = 501)
cbpp_prof <- lapply(cbpp_mod_list, p_fun)

cbpp_confint_prof <- lapply(cbpp_prof, confint)
cbpp_confint_wald <- lapply(cbpp_mod_list, wald_cifun)
cbpp_est <- do.call("rbind", Map(get_est, cbpp_mod_list, cbpp_df_name$mnames))
rownames(cbpp_est) <- NULL

cbpp_combCI <- Map(combfun2,
                   list(cbpp_confint_wald, cbpp_confint_boot, cbpp_confint_prof),
                   c("Wald", "boot", "profile"),
                   MoreArgs = list(df_est = cbpp_est, df_name = cbpp_df_name)) |>
  do.call(what = "rbind")

if (FALSE) {
  
  library(ggplot2); theme_set(theme_bw())
  revguide <- guide_legend(reverse = TRUE)
  des_ord <- c("sd_(Intercept)|herd", "sd_(Intercept)|obs",
               "avg_size", "treatmentUnknown", "treatmentComplete")
  combCI <- subset(cbpp_combCI,var!="(Intercept)" & !grepl("^period",var))
  combCI <- transform(combCI,  var = factor(var, levels = rev(des_ord)))
  combCI <- transform(combCI,
                      vartype = factor(ifelse(grepl("^sd", var), "random effects", "fixed effects")))
  mnames <- cbpp_df_name$mnames
  pd <- position_dodge(width=0.6)
  ggplot(combCI, aes(var, est, colour=type, shape=model)) +
    geom_pointrange(aes(ymin=lwr, ymax=upr),
                    position=pd) +
    labs(x="", y="Estimate (log-odds of seropositivity)") +
    geom_hline(yintercept=0, lty=2) +
    coord_flip() +
    scale_colour_brewer(palette="Dark2", guide = revguide) +
    scale_shape(guide=revguide, labels = mnames) +
    facet_wrap(~ vartype, ncol = 1, scale = "free")

  ggsave("tmp.png")
}

objs <- c("cbpp_confint_prof", "cbpp_confint_boot", "cbpp_confint_wald",
          ## "cbpp_prof",
          "cbpp_est", "cbpp_combCI", "cbpp_df_name")


save(
  list = objs,
  file = "cbpp_batch.rda",
  version = 2
)
sessionInfo()
