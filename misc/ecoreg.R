download.file("https://github.com/bbolker/mmd_utils/raw/refs/heads/master/ecoreg.rds", dest = "ecoreg.rds")
ee <- readRDS("ecoreg.rds")

form <- mmamm_log ~
    (NPP_log_sc + Feat_log_sc + NPP_cv_sc + Feat_cv_sc)^2 +
    (NPP_log_sc + Feat_log_sc + NPP_cv_sc + Feat_cv_sc | biome)

library(lme4)
ee1 <- lmer(form, ee, REML = FALSE)
isSingular(ee1)


form5 <- mmamm_log ~
  (NPP_log_sc + Feat_log_sc + NPP_cv_sc + Feat_cv_sc)^2 +
  diag(NPP_log_sc + Feat_log_sc + NPP_cv_sc + Feat_cv_sc | biome)

lmer(form5, ee, REML = FALSE)

form6 <- mmamm_log ~
  (NPP_log_sc + Feat_log_sc + NPP_cv_sc + Feat_cv_sc)^2 +
  cs(NPP_log_sc + Feat_log_sc + NPP_cv_sc + Feat_cv_sc | biome)

mm <- lmer(form6, ee, REML = FALSE)
