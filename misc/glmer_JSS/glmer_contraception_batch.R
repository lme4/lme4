library(lme4)
data("Contraception", package = "mlmRev")

get_est <- function(mod, model_name) {
  v <- as.data.frame(VarCorr(mod))
  vnames <- with(v,
                 ifelse(is.na(var2),
                        paste0("sd_", var1, "|", grp),
                        paste0("cor_", var1, ".", var2, "|", grp)))
  est <- c(fixef(mod), setNames(v$sdcor, vnames))
  data.frame(model=model_name, var=names(est), est=est, row.names=NULL)
}

b_cifun <- function(x) {
  cat(deparse(getCall(x)$formula), "\n")
  confint(x, method = "boot", seed = 101, nsim = 501, signames = FALSE,
          parallel = "multicore", ncpus = 10)
}

p_fun <- function(x) {
  cat(deparse(getCall(x)$formula), "\n")
  profile(x, signames = FALSE,
          devtol = 1e-2,
          maxpts = 250,
          parallel = "multicore", ncpus = 10)
}

Contraception <- transform(Contraception,
                           ch = factor(Contraception$livch != 0, labels = c("N","Y")),
                           age = age/(2*sd(age)))

df_name <- data.frame(
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

mod_list <- lapply(mforms, glmer, data = Contraception, family = binomial)
names(mod_list) <- df_name$mnames

contr_est <- do.call("rbind", Map(get_est, mod_list, df_name$mnames))
rownames(contr_est) <- NULL

contr_confint_wald <- lapply(mod_list, function(x) confint(x, method = "Wald", signames = FALSE))
contr_prof <- lapply(mod_list, p_fun)
contr_confint_prof <- lapply(contr_prof, confint)
contr_confint_boot <- lapply(mod_list, b_cifun)


stop()

## structured covariances: fit/store separately?
## (much faster than profile/boot/etc, if they work ...)
contr.diag <- glmer(use ~ diag(1 + age|district), Contraception, binomial)
contr.hetdiag <- glmer(use ~ diag(1 + age|district, hom = FALSE), 
                       data = Contraception, family = binomial)
contr.cs <- glmer(use ~ cs(1 + age|district), Contraception, binomial)
contr.hetcs <- glmer(use ~ cs(1 + age|district, hom = FALSE), Contraception, binomial)

vars <- c(ls(pattern = "contr.*"), "df_name")
save(list = vars, file = "Contraception_batch.rda", version = 2)
