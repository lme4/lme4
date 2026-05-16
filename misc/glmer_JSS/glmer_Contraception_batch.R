library(lme4)
data("Contraception", package = "mlmRev")
source("glmer_batch_funs.R")

Contraception <- transform(Contraception,
                           ch = factor(Contraception$livch != 0, labels = c("N","Y")),
                           age = age/(2*sd(age)))

options(contrasts = c("contr.sum", "contr.poly"))

contr_df_name <- data.frame(
  mnames = c("ic", "bc", "bc_int", "bc_int_dist_urb_var", "bc_int_dist_urb_nest", "bc_int_dist_urb_cross"),
  mdesc =
    c("int_child    + age  + (1 | district)",
    "binary_child + age  + (1 | district)",
    "binary_child × age + (1 | district)",
    "binary_child × age + (1 + urban | district)",
    "binary_child × age + (1 | district/urban)",
    "binary_child × age + (1 | district:urban"))

mforms <- list()
mforms$ic <- use ~ age + I(age^2) + urban + livch + (1|district)
mforms$bc <- update(mforms$ic, . ~ . - livch + ch)
mforms$bc_int_urb <- update(mforms$bc, . ~ . + age:ch)
mforms$bc_int_dist_urb_var <- update(mforms$bc_int_urb, . ~ . - (1|district) + (1 + urban | district))
mforms$bc_int_dist_urb_nest <- update(mforms$bc_int_urb, . ~ . - (1|district) + (1 | district/urban))
mforms$bc_int_dist_urb_cross <- update(mforms$bc_int_urb, . ~ . - (1|district) + (1 | urban:district))

contr_mod_list <- lapply(mforms, glmer, data = Contraception, family = binomial)

contr_est <- do.call("rbind", Map(get_est, contr_mod_list, contr_df_name$mnames))
rownames(contr_est) <- NULL

contr_confint_wald <- lapply(contr_mod_list, function(x) confint(x, method = "Wald", signames = FALSE))
contr_prof <- lapply(contr_mod_list, p_fun)
contr_confint_prof <- lapply(contr_prof, confint)
contr_confint_boot <- lapply(contr_mod_list, b_cifun)

contr_combCI <- Map(combfun2,
                   list(contr_confint_wald, contr_confint_boot, contr_confint_prof),
                   c("Wald", "boot", "profile"),
                   MoreArgs = list(df_est = contr_est, df_name = contr_df_name)) |>
  do.call(what = "rbind")

## structured covariances: fit/store separately?
## (much faster than profile/boot/etc, if they work ...)
contr.diag <- glmer(use ~ diag(1 + age|district), Contraception, binomial)
contr.hetdiag <- glmer(use ~ diag(1 + age|district, hom = FALSE), 
                       data = Contraception, family = binomial)
contr.cs <- glmer(use ~ cs(1 + age|district), Contraception, binomial)
contr.hetcs <- glmer(use ~ cs(1 + age|district, hom = FALSE), Contraception, binomial)

## leave contr_prof out of save_list until we have better idea how to 'butcher' (reduce size/exclude environments when serializing)
## (we can get by without it

## inspired by butcher:::axe_env()

## axe_env.thpr <- function(x, verbose = FALSE, ...) {
##   a <- attributes(x)
##   for (d in c("forward", "backward")) {
##     for (i in seq_along(a[[d]])) {
##       cat(i, "\n")
##       attr(attr(a[[d]][[i]], "formula"), ".Environment") <- baseenv()
##     }
##   }
##   attributes(x) <- a
##   x
## }
## x_clean <- axe_env.thpr(x)


save(
  list = c("contr_confint_prof", "contr_confint_boot", "contr_confint_wald", "contr_est", "contr_combCI", "contr_df_name"),
  file = "Contraception_batch.rda",
  version = 2
)

