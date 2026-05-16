## Parallel to glmer_cbpp_batch.R but using nloptwrap instead of bobyqa.
## nloptwrap wraps the nloptr package (NLOPT_LN_BOBYQA by default) and uses
## a different code path from minqa's bobyqa — useful for diagnosing whether
## profile CI behaviour is optimizer-dependent.
## Note: the original batch uses maxfun = 1000 (very low); this batch uses
## nloptwrap defaults so convergence is not artificially truncated.
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
  mdesc  = c("herd", "herdobs", "obs"))

mforms <- list()
mforms$herd    <- incidence / size ~ period + treatment + avg_size + (1 | herd)
mforms$herdobs <- update(mforms$herd, . ~ . + (1 | obs))
mforms$obs     <- update(mforms$herd, . ~ . - (1|herd) + (1 | obs))

## for loop required (see comment in glmer_cbpp_batch.R re: weights + lapply)
cbpp_mod_list_altopt <- list()
for (i in names(mforms)) {
  cbpp_mod_list_altopt[[i]] <- glmer(mforms[[i]],
                                     data    = cbpp2,
                                     family  = binomial,
                                     weights = size,
                                     control = glmerControl(optimizer = "nloptwrap"))
}

cbpp_confint_boot_altopt <- lapply(cbpp_mod_list_altopt, b_cifun)
cbpp_prof_altopt         <- lapply(cbpp_mod_list_altopt, p_fun)
cbpp_confint_prof_altopt <- lapply(cbpp_prof_altopt, confint)
cbpp_confint_wald_altopt <- lapply(cbpp_mod_list_altopt,
                                   function(x) confint(x, method = "Wald", signames = FALSE))
cbpp_est_altopt <- do.call("rbind", Map(get_est, cbpp_mod_list_altopt, cbpp_df_name$mnames))
rownames(cbpp_est_altopt) <- NULL

cbpp_combCI_altopt <- Map(combfun2,
                          list(cbpp_confint_wald_altopt,
                               cbpp_confint_boot_altopt,
                               cbpp_confint_prof_altopt),
                          c("Wald", "boot", "profile"),
                          MoreArgs = list(df_est  = cbpp_est_altopt,
                                         df_name = cbpp_df_name)) |>
  do.call(what = "rbind")

save(
  list = c("cbpp_confint_prof_altopt", "cbpp_confint_boot_altopt", "cbpp_confint_wald_altopt",
           "cbpp_prof_altopt", "cbpp_est_altopt", "cbpp_combCI_altopt", "cbpp_df_name"),
  file = "cbpp_altopt_batch.rda",
  version = 2
)
sessionInfo()
