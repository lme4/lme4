library("lme4")
library("testthat")

context("testing '||' notation for independent ranefs")

test_that("basic intercept + slope '||' works", {
	expect_equivalent(
		lFormula(Reaction ~ Days + (Days||Subject), sleepstudy)$reTrms,
		lFormula(Reaction ~ Days + (1|Subject) + (0 + Days|Subject), sleepstudy)$reTrms
	)
	
	expect_equivalent(
		fitted(lmer(Reaction ~ Days + (Days||Subject), sleepstudy)),
		fitted(lmer(Reaction ~ Days + (1|Subject) + (0 + Days|Subject), sleepstudy))
	)
})
