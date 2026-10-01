## GH #289: the toenail Laplace fit disagreed with other implementations.
test_that("the toenail fit agrees with the Laplace results reported in GH #289", {
  m <- glmer(
    outcome ~ treatment + visit + (1 | patientID),
    toenail,
    family = binomial
  )
  ## References are the glmmML results in the issue, using default glmer controls.
  ## Before the refresh, the random-intercept SD was 4.68738 and logLik -626.1709.
  expect_equal(unname(getME(m, "theta")), 4.708440, tolerance = 1e-4)
  expect_equal(
    unname(fixef(m)),
    c(-1.070209, -0.7005596, -0.9126392),
    tolerance = 1e-4
  )
  expect_equal(as.numeric(logLik(m)), -626.1627, tolerance = 1e-6)
})
