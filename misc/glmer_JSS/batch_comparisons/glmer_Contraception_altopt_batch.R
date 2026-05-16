## Parallel to glmer_Contraception_batch.R but using nloptwrap instead of
## the default bobyqa optimizer. See glmer_cbpp_altopt_batch.R for rationale.
library(lme4)
data("Contraception", package = "mlmRev")
source("glmer_batch_funs.R")

Contraception <- transform(Contraception,
                           ch  = factor(Contraception$livch != 0, labels = c("N","Y")),
                           age = age / (2 * sd(age)))

options(contrasts = c("contr.sum", "contr.poly"))

contr_df_name <- data.frame(
  mnames = c("ic", "bc", "bc_int", "bc_int_dist_urb_var", "bc_int_dist_urb_nest", "bc_int_dist_urb_cross"),
  mdesc  =
    c("int_child    + age  + (1 | district)",
      "binary_child + age  + (1 | district)",
      "binary_child × age + (1 | district)",
      "binary_child × age + (1 + urban | district)",
      "binary_child × age + (1 | district/urban)",
      "binary_child × age + (1 | district:urban"))

mforms <- list()
mforms$ic                    <- use ~ age + I(age^2) + urban + livch + (1|district)
mforms$bc                    <- update(mforms$ic, . ~ . - livch + ch)
mforms$bc_int_urb            <- update(mforms$bc, . ~ . + age:ch)
mforms$bc_int_dist_urb_var   <- update(mforms$bc_int_urb, . ~ . - (1|district) + (1 + urban | district))
mforms$bc_int_dist_urb_nest  <- update(mforms$bc_int_urb, . ~ . - (1|district) + (1 | district/urban))
mforms$bc_int_dist_urb_cross <- update(mforms$bc_int_urb, . ~ . - (1|district) + (1 | urban:district))

contr_mod_list_altopt <- lapply(mforms, function(f)
  glmer(f, data = Contraception, family = binomial,
        control = glmerControl(optimizer = "nloptwrap")))

contr_est_altopt <- do.call("rbind", Map(get_est, contr_mod_list_altopt, contr_df_name$mnames))
rownames(contr_est_altopt) <- NULL

contr_confint_wald_altopt <- lapply(contr_mod_list_altopt,
                                    function(x) confint(x, method = "Wald", signames = FALSE))
contr_prof_altopt         <- lapply(contr_mod_list_altopt, p_fun)
contr_confint_prof_altopt <- lapply(contr_prof_altopt, confint)
contr_confint_boot_altopt <- lapply(contr_mod_list_altopt, b_cifun)

contr_combCI_altopt <- Map(combfun2,
                           list(contr_confint_wald_altopt,
                                contr_confint_boot_altopt,
                                contr_confint_prof_altopt),
                           c("Wald", "boot", "profile"),
                           MoreArgs = list(df_est  = contr_est_altopt,
                                          df_name = contr_df_name)) |>
  do.call(what = "rbind")

save(
  list = c("contr_confint_prof_altopt", "contr_confint_boot_altopt", "contr_confint_wald_altopt",
           "contr_prof_altopt", "contr_est_altopt", "contr_combCI_altopt", "contr_df_name"),
  file = "Contraception_altopt_batch.rda",
  version = 2
)
sessionInfo()
