# testing "||" notation for independent ranefs
# setwd("C:/birs/lme4")
# library(devtools)
# rm(list=ls())
# load_all()


library(testthat)

#basic intercept + slope
expect_equivalent(
  lFormula(Reaction ~ Days + (Days||Subject), sleepstudy)$reTrms,
  lFormula(Reaction ~ Days + (1|Subject) + (0 + Days|Subject), sleepstudy)$reTrms,
)  

expect_equivalent(
  fitted(lmer(Reaction ~ Days + (Days||Subject), sleepstudy)),
  fitted(lmer(Reaction ~ Days + (1|Subject) + (0 + Days|Subject), sleepstudy)),
)  


#works with nested
expect_equivalent(findbars(y ~ (x || id / id2)),
                  findbars(y ~ (1 | id  / id2) + (0 + x | id / id2)))

#works with multiple
expect_equivalent(findbars(y ~ (x1 + x2  || id / id2) + (x3 | id3) + (x4 || id4)),
                  findbars(y ~ (1 | id / id2) + (0 + x1 | id / id2) +
                             (0 + x2 | id / id2) + (x3 | id3) + (1 | id4) + 
                             (0 + x4| id4)))

# leaves superfluous || alone:
expect_equivalent(findbars(y ~ z + (0 + x || id)),
                  findbars(y ~ z + (0 + x  | id)))

# plays nice with parens in fixed formula:
expect_equivalent(findbars(y ~ (z + x)^2 + (x || id)),
                  findbars(y ~ (z + x)^2 + (1 | id) + (0 + x | id)))

# works for interactions:
expect_equivalent(findbars(y ~ (x1*x2 || id)),
                  findbars(y ~ (1 | id) + (0+x1 | id) + (0 + x2 | id) +
                             (0 + x1:x2 | id)))
      
## update now works as expected:
(m <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy))
update(m, .~.-(0 + Days | Subject))
   
