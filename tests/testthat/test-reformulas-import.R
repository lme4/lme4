## these tests may fail if run a second time in the same session, because of rlang only-report-once-per-session setting
f <- ~ 1 + (1|f)
test_that("deprecation warnings from lme4 formula processing", {
  expect_warning(subbars(f), "has moved")
  expect_warning(nobars(f), "has moved")
  expect_warning(findbars(f), "has moved")
  expect_warning(mkReTrms(findbars(f), fr = data.frame(f = factor(1:10))), "has moved")
  expect_warning(expandDoubleVerts(f), "has moved")
})
