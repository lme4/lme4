library("glmmTMB")
cbpp2 <- read.csv("cbpp2.csv")
source("glmmTMB_batch_funs.R")

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
  mdesc  = c("herd", "herdobs", "obs"))

mforms <- list()
mforms$herd    <- incidence / size ~ period + treatment + avg_size + (1 | herd)
mforms$herdobs <- update(mforms$herd, . ~ . + (1 | obs))
mforms$obs     <- update(mforms$herd, . ~ . - (1|herd) + (1 | obs))

cbpp_mod_list_tmb <- list()
for (i in names(mforms)) {
  cbpp_mod_list_tmb[[i]] <- glmmTMB(mforms[[i]],
                                     data    = cbpp2,
                                     family  = binomial,
                                     weights = size)
}

cbpp_confint_wald_tmb <- lapply(cbpp_mod_list_tmb, w_cifun_tmb)
cbpp_confint_ur_tmb   <- lapply(cbpp_mod_list_tmb, u_cifun_tmb)
cbpp_confint_prof_tmb <- lapply(cbpp_mod_list_tmb, p_cifun_tmb)
cbpp_est_tmb <- do.call("rbind", Map(get_est_tmb, cbpp_mod_list_tmb, cbpp_df_name$mnames))
rownames(cbpp_est_tmb) <- NULL

cbpp_combCI_tmb <- Map(combfun2,
                       list(cbpp_confint_wald_tmb, cbpp_confint_ur_tmb, cbpp_confint_prof_tmb),
                       c("Wald", "uniroot", "profile"),
                       MoreArgs = list(df_est = cbpp_est_tmb, df_name = cbpp_df_name)) |>
  do.call(what = "rbind")

save(
  list = c("cbpp_confint_wald_tmb", "cbpp_confint_ur_tmb", "cbpp_confint_prof_tmb",
           "cbpp_est_tmb", "cbpp_combCI_tmb", "cbpp_df_name"),
  file = "cbpp_tmb_batch.rda",
  version = 2
)
sessionInfo()
