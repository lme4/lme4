#context("testing '||' notation for independent ranefs")

form1 <- Reaction ~ Days + (Days||Subject)
form2 <- Reaction ~ Days + (1|Subject) + (0 + Days|Subject)
m1 <- lmer(form1, sleepstudy)
m2 <- lmer(form2, sleepstudy)

test_that("basic intercept + slope '||' works", {

  ## no longer expect exactly the same structures (diag() vs two terms)
	## expect_equivalent(
	## 	lFormula(form1, sleepstudy)$reTrms,
	## 	lFormula(form2, sleepstudy)$reTrms
	## )
	
	expect_equivalent(
		fitted(m1),
		fitted(m2)
	)

  ## fitted values from glmmTMB
  ## library(glmmTMB)
  ## m3 <- glmmTMB(form1, sleepstudy, REML = TRUE, start = list(theta = log(c(25,6))))
  ## dput(head(fitted(m3), 4))

  expect_equivalent(head(fitted(m2), 4),
                    c(252.917801422647, 272.708576602026, 292.499351781404, 312.290126960783),
                    tolerance = 5e-5)

})


