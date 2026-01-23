library("testthat")

test_that("checking singular fit for covariance", {
  
  set.seed(1)
  test.us <- new("Covariance.us", nc = 4L, simulate = TRUE)
  expect_equal(isSingular(test.us), FALSE)
  
  test.diag <- new("Covariance.diag", nc = 4L, simulate = TRUE)
  expect_equal(isSingular(test.diag), FALSE)
  
  test.cs <- new("Covariance.cs", nc = 4L, hom = FALSE, simulate = TRUE)
  test.homcs <- new("Covariance.cs", nc = 4L, hom = TRUE, simulate = TRUE)
  expect_equal(isSingular(test.cs), FALSE)
  expect_equal(isSingular(test.homcs), FALSE)
  
  test.ar1 <- new("Covariance.ar1", nc = 4L, hom = TRUE, simulate = TRUE)
  test.hetar1 <- new("Covariance.ar1", nc = 4L, hom = FALSE, simulate = TRUE)
  expect_equal(isSingular(test.ar1), FALSE)
  expect_equal(isSingular(test.hetar1), FALSE)
  
  test.reCovs <- c(test.us, test.diag, test.cs, test.homcs, test.ar1, test.hetar1)
  
  ## tests should fail once one of the diagional entries is 0
  sapply(test.reCovs, function(i){
    new_theta <- getTheta(i)
    new_theta[1] <- 0
    test <- setTheta(i, new_theta)
    expect_equal(isSingular(test), TRUE)
  })
  
  test.reCovs.lim <- test.reCovs[-c(1:2)]
  
  ## tests should fail once rho is equal to -1 or 1
  sapply(test.reCovs.lim, function(i){
    test <- setVC(i, vcomp = getVC(i)$vcomp, ccomp = 1)
    expect_equal(isSingular(test), TRUE)
    test2 <- setVC(i, vcomp = getVC(i)$vcomp, ccomp = -1)
    expect_equal(isSingular(test2), TRUE)
  })
})

test_that("checking singular fit for merMod", {

  set.seed(101)
  
  n_groups <- 20
  n_per_group <- 20
  n <- n_groups * n_per_group
  
  dat <- data.frame(
    group1 = rep(1:n_groups, each = n_per_group),
    group2 = rep(1:n_groups, each = n_per_group),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  form <- y ~ 1 + (1 + x1|group1) + (1 + x2|group2)
  dat$y <- simulate(form[-2], ## one-sided formula
                    newdata = dat,
                    family = gaussian,
                    newparams = list(beta = c(-2),
                                     theta = c(2, 3, 4, 2, 3, 4),
                                     sigma = 2))[[1]]
  mod <- lmer(form, data = dat)
  ## should be no issues here, ensuring isSingular is not broken
  expect_equal(isSingular(mod), FALSE)
  
  form2 <- y2 ~ 1 + (1 + x1|group1) + (1 + x2|group2)
  dat$y2 <- simulate(form[-2], 
                    newdata = dat,
                    family = gaussian,
                    newparams = list(beta = c(-2),
                                     theta = c(0, 0, 0, 2.5, 1.5, 0.8),
                                     sigma = 2))[[1]]
  mod2 <- suppressWarnings(lmer(form2, data = dat))
  ## second one, should be an issue with group 1 as we set thetas = 0
  expect_equal(isSingular(mod2), TRUE)
  
  ## also need to consider the glmer case
  form <- y ~ 1 + (1 + x1 | group1) + (1 + x2 | group2)
  
  dat$y <- simulate(form[-2], newdata = dat,
                    family = poisson(link = "log"), 
                    newparams = list(theta = c(0, 0, 0, 1.0, 0.5, 0.3),
                                     beta = 2))[[1]]
  
  form2 <- y2 ~ 1 + (1 + x1 | group1) + (1 + x2 | group2)
  
  dat$y2 <- simulate(form2[-2], newdata = dat,
                     family = Gamma(link = "log"), 
                     newparams = list(theta = c(0, 0, 0, 1.0, 0.5, 0.3),
                                      beta = 2,
                                      sigma = 2))[[1]]
  
  
  mod_pois <- suppressWarnings(
    glmer(form, family = poisson(link = "log"), data = dat))
  
  mod_gam <- suppressWarnings(
    glmer(form2, family = Gamma(link = "log"), data = dat))
  
  expect_equal(isSingular(mod_pois), TRUE)
  expect_equal(isSingular(mod_gam), TRUE)
})
