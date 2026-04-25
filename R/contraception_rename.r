# This is for renaming the batch into more explanable names.

# ── Batch 3: boot CIs with signames=TRUE ──────────────────────────────────────
load("./misc/glmer_JSS/glmer_profbatch3.RData")   # loads confint.boot.cm1_2, etc.

# Rename each object to flag signames=TRUE
confint.boot.cm1_2.sigT <- confint.boot.cm1_2
confint.boot.cm2_2.sigT <- confint.boot.cm2_2
confint.boot.cm3_2.sigT <- confint.boot.cm3_2
confint.boot.cm4_2.sigT <- confint.boot.cm4_2
confint.boot.cm5_2.sigT <- confint.boot.cm5_2
confint.boot.cm6_2.sigT <- confint.boot.cm6_2

save(list = ls(pattern = "confint\\.boot\\..*\\.sigT"),
     file = "./misc/glmer_JSS/Contra_boot_sigTRUE.RData")

rm(list = ls())


# ── Batch 4: prof CIs (signames=FALSE) + wald CIs (signames=TRUE) ─────────────
load("./misc/glmer_JSS/glmer_profbatch4.RData")

# prof was always signames=FALSE — no suffix change needed, but be explicit
confint.prof.cm3_2.sigF <- confint.prof.cm3_2
confint.prof.cm4_2.sigF <- confint.prof.cm4_2
confint.prof.cm5_2.sigF <- confint.prof.cm5_2
confint.prof.cm6_2.sigF <- confint.prof.cm6_2

save(list = ls(pattern = "confint\\.prof\\..*\\.sigF"),
     file = "./misc/glmer_JSS/Contra_prof_sigFALSE.RData")

# wald here was signames=TRUE
confint.wald.cm1_2.sigT <- confint.wald.cm1_2
confint.wald.cm2_2.sigT <- confint.wald.cm2_2
confint.wald.cm3_2.sigT <- confint.wald.cm3_2
confint.wald.cm4_2.sigT <- confint.wald.cm4_2
confint.wald.cm5_2.sigT <- confint.wald.cm5_2
confint.wald.cm6_2.sigT <- confint.wald.cm6_2

save(list = ls(pattern = "confint\\.wald\\..*\\.sigT"),
     file = "./misc/glmer_JSS/Contra_wald_sigTRUE.RData")

rm(list = ls())


# ── Batch 5: boot + wald CIs, both signames=FALSE (current run) ───────────────
load("./misc/glmer_JSS/glmer_profbatch5.RData")

confint.boot.cm1_2.sigF <- confint.boot.cm1_2
confint.boot.cm2_2.sigF <- confint.boot.cm2_2
confint.boot.cm3_2.sigF <- confint.boot.cm3_2
confint.boot.cm4_2.sigF <- confint.boot.cm4_2
confint.boot.cm5_2.sigF <- confint.boot.cm5_2
confint.boot.cm6_2.sigF <- confint.boot.cm6_2

confint.wald.cm1_2.sigF <- confint.wald.cm1_2
confint.wald.cm2_2.sigF <- confint.wald.cm2_2
confint.wald.cm3_2.sigF <- confint.wald.cm3_2
confint.wald.cm4_2.sigF <- confint.wald.cm4_2
confint.wald.cm5_2.sigF <- confint.wald.cm5_2
confint.wald.cm6_2.sigF <- confint.wald.cm6_2

save(list = ls(pattern = "confint\\.boot\\..*\\.sigF"),
     file = "./misc/glmer_JSS/Contra_boot_sigFALSE.RData")

save(list = ls(pattern = "confint\\.wald\\..*\\.sigF"),
     file = "./misc/glmer_JSS/Contra_wald_sigFALSE.RData")
