library(glmmTMB)
data("Contraception", package = "mlmRev")
source("glmmTMB_batch_funs.R")

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

contr_mod_list_tmb <- list()
for (i in names(mforms)) {
  contr_mod_list_tmb[[i]] <- glmmTMB(mforms[[i]], data = Contraception, family = binomial)
}

contr_est_tmb <- do.call("rbind", Map(get_est_tmb, contr_mod_list_tmb, contr_df_name$mnames))
rownames(contr_est_tmb) <- NULL

contr_confint_wald_tmb <- lapply(contr_mod_list_tmb, w_cifun_tmb)
contr_confint_ur_tmb   <- lapply(contr_mod_list_tmb, u_cifun_tmb)
contr_confint_prof_tmb <- lapply(contr_mod_list_tmb, p_cifun_tmb)

contr_combCI_tmb <- Map(combfun2,
                        list(contr_confint_wald_tmb, contr_confint_ur_tmb, contr_confint_prof_tmb),
                        c("Wald", "uniroot", "profile"),
                        MoreArgs = list(df_est = contr_est_tmb, df_name = contr_df_name)) |>
  do.call(what = "rbind")

save(
  list = c("contr_confint_wald_tmb", "contr_confint_ur_tmb", "contr_confint_prof_tmb",
           "contr_combCI_tmb", "contr_est_tmb", "contr_df_name"),
  file = "Contraception_tmb_batch.rda",
  version = 2
)
sessionInfo()
