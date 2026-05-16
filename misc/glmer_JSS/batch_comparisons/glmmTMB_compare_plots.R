## Comparison plots: glmer vs glmmTMB confidence intervals
## Run glmer_cbpp_batch.R, glmer_Contraception_batch.R,
##     glmmTMB_cbpp_batch.R, glmmTMB_Contraception_batch.R first.
##
## CI type mapping: glmer {Wald, boot, profile} <-> glmmTMB {Wald, uniroot, profile}
## glmmTMB uniroot = profile-likelihood root-finding (comparable to lme4 profile).
## Scalar RE CIs (profile/uniroot) are back-transformed from log(sd) to sd scale
## by backtransform_theta() in glmmTMB_batch_funs.R.
## Multi-RE group off-diagonal parameters (Cholesky scale) are not shown.
library(ggplot2)

revguide <- guide_legend(reverse = TRUE)

## ---- CBPP ---------------------------------------------------------------

load("cbpp_batch.rda")       # cbpp_combCI, cbpp_df_name
load("cbpp_tmb_batch.rda")   # cbpp_combCI_tmb (cbpp_df_name overwritten, same values)

cbpp_combCI$package     <- "glmer"
cbpp_combCI_tmb$package <- "glmmTMB"

cbpp_combined <- rbind(
  cbpp_combCI[,     c("var","lwr","upr","est","model","model_id","type","package")],
  cbpp_combCI_tmb[, c("var","lwr","upr","est","model","model_id","type","package")]
)

des_ord <- c("sd_(Intercept)|herd", "sd_(Intercept)|obs",
             "avg_size", "treatmentUnknown", "treatmentComplete")

cbpp_comp <- subset(cbpp_combined, var %in% des_ord)
cbpp_comp <- transform(cbpp_comp,
  vartype = factor(ifelse(grepl("^sd", var), "random effects", "fixed effects")),
  var     = factor(var, levels = rev(des_ord)),
  model   = factor(model, levels = cbpp_df_name$mdesc),
  package = factor(package, levels = c("glmer", "glmmTMB"))
)

mnames_cbpp <- cbpp_df_name$mnames
pd <- position_dodge(width = 0.6)

ggplot(cbpp_comp, aes(var, est, colour = type, shape = model)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr), position = pd) +
  labs(x = "", y = "Estimate (log-odds of seropositivity)") +
  geom_hline(yintercept = 0, lty = 2) +
  coord_flip() +
  scale_colour_brewer(palette = "Dark2", guide = revguide) +
  scale_shape(guide = revguide, labels = mnames_cbpp) +
  facet_grid(vartype ~ package, scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(strip.text.y = element_text(size = 7))

## ---- Contraception ------------------------------------------------------

load("Contraception_batch.rda")       # contr_combCI, contr_df_name
load("Contraception_tmb_batch.rda")   # contr_combCI_tmb (contr_df_name overwritten, same values)

contr_combCI$package     <- "glmer"
contr_combCI_tmb$package <- "glmmTMB"

contr_combined <- rbind(
  contr_combCI[,     c("var","lwr","upr","est","model","model_id","type","package")],
  contr_combCI_tmb[, c("var","lwr","upr","est","model","model_id","type","package")]
)

contr_comp <- subset(contr_combined, var != "(Intercept)")
contr_comp <- transform(contr_comp,
  vartype = factor(ifelse(grepl("^(sd|cor)", var), "random effects", "fixed effects")),
  model   = factor(model, levels = rev(contr_df_name$mdesc)),
  package = factor(package, levels = c("glmer", "glmmTMB"))
)

## Split by vartype so each panel can have its own y-axis limits via coord_flip(ylim).
## A single facet_grid cannot set different limits per row without extra packages.
## coord_flip(ylim = ...) clips (zooms) without dropping data points.

contr_base <- function(dat) {
  ggplot(dat, aes(var, est, colour = type, shape = model)) +
    geom_pointrange(aes(ymin = lwr, ymax = upr), position = pd) +
    geom_hline(yintercept = 0, lty = 2) +
    scale_colour_brewer(palette = "Dark2", guide = revguide) +
    scale_shape(guide = revguide) +
    facet_grid(vartype ~ package) +
    labs(x = "") +
    theme_bw() +
    theme(strip.text.y = element_text(size = 7))
}

p_contr_fe <- contr_base(subset(contr_comp, vartype == "fixed effects")) +
  coord_flip(ylim = c(-2.5, 0.5)) +
  labs(y = "Estimate (log-odds of contraceptive use)")

p_contr_re <- contr_base(subset(contr_comp, vartype == "random effects")) +
  coord_flip(ylim = c(0, 0.8)) +
  labs(y = "") +
  theme(strip.text.x = element_blank(),   # hide redundant package labels
        legend.position = "none")

gridExtra::grid.arrange(p_contr_fe, p_contr_re, ncol = 1, heights = c(3, 1))
