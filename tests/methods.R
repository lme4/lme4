library(lme4)
library(testthat)
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
expect_equal(colnames(model.frame(fm1)),c("Reaction","Days","Subject"))
expect_equal(colnames(model.frame(fm1,fixed.only=TRUE)),c("Reaction","Days"))
expect_equal(formula(fm1),Reaction ~ Days + (Days | Subject))
expect_equal(formula(fm1,fixed.only=TRUE),Reaction ~ Days)

## ugly example: model frame with compound elements
fm2 <- lmer(log(Reaction) ~ splines::ns(Days,3) +
            + I(1+Days^3) + (Days|Subject), sleepstudy)
expect_equal(names(model.frame(fm2)),
             c("log(Reaction)", "splines::ns(Days, 3)",
               "I(1 + Days^3)", "Days", "Subject"))
expect_equal(names(model.frame(fm2,fixed.only=TRUE)),
             c("log(Reaction)", "splines::ns(Days, 3)",
               "I(1 + Days^3)"))

