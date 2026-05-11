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

cbpp_mod_list <- list()
for (i in names(mforms)) {
  cbpp_mod_list[[i]] <- glmer(mforms[[i]],
                              data = cbpp2,
                              family = binomial,
                              weights = size,
                              control = glmerControl(optimizer = "bobyqa",
                                                     optCtrl = list(maxfun = 1000)))
}

cbpp_confint_boot <- lapply(cbpp_mod_list, b_cifun)
cbpp_prof <- lapply(cbpp_mod_list, p_fun)

if (FALSE) {
  ## experimenting/exploring
  profile.gm1 <- cbpp_confint_prof$herd
  lattice::xyplot(profile.gm1)

  library(ggplot2)
  p2 <- as.data.frame(profile.gm1)
  ggplot(p2, aes(.focal, .zeta)) +
    geom_point() +
    geom_line() +
    facet_wrap(~.par, scale = "free")
}

cbpp_confint_prof <- lapply(cbpp_prof, confint)
cbpp_confint_wald <- lapply(cbpp_mod_list, function(x) confint(x, method = "Wald", signames = FALSE))
cbpp_est <- do.call("rbind", Map(get_est, cbpp_mod_list, cbpp_df_name$mnames))
rownames(cbpp_est) <- NULL

cbpp_combCI <- Map(combfun2,
                   list(cbpp_confint_wald, cbpp_confint_boot, cbpp_confint_prof),
                   c("Wald", "boot", "profile"),
                   MoreArgs = list(df_est = cbpp_est, df_name = cbpp_df_name)) |>
  do.call(what = "rbind")

## ggplot(combCI, aes(est, var, colour = type)) + geom_point(position = position_dodge(0.5)) + facet_wrap(~model)


save(
  list = c("cbpp_confint_prof", "cbpp_confint_boot", "cbpp_confint_wald", "cbpp_prof", "cbpp_est", "cbpp_combCI", "cbpp_df_name"),
  file = "CBPP_batch.rda",
  version = 2
)
sessionInfo()
