library(lme4)
print(packageVersion("Matrix"))
m1 <- lmer(Reaction 1~ Days + (Days|Subject), data = sleepstudy)
predict(m1, newdata = sleepstudy, se.fit = TRUE)
simulate(~Days + (Days|Subject),
         newdata = sleepstudy,
         newparams = list(beta = c(250, 1),
                          theta = c(1, 1, 1),
                          sigma = 1),
         family = gaussian)

set.seed(101)
dd <- expand.grid(id = factor(1:1000), year = 1:15)
d_age <- data.frame(id = factor(1:1000),
                    age0 = sample(20:50, size = 1000, replace= TRUE))
dd2 <- merge(dd, d_age, by = "id") |>
    transform(age = age0 + year,
              score = sample(1:50, size = nrow(dd), replace = TRUE))
form <- adl~I(year^2)+score*year + age*year + (1+year|id)
colnames(X <- model.matrix(nobars(form[-2]), data = dd2))
beta_vec <- c(20, 0.1, 0.5, 1, 1, 0.5, 0.5)
stopifnot(ncol(X) == length(beta_vec))
dd2$adl <- simulate(form[-2],
                    newdata = dd2,
                    newparams = list(beta = beta_vec,
                                     theta = c(2, 0.1, 0),
                                     sigma = 1),
                    family = gaussian)[[1]]

fit <- lmer(form, data = dd2)

pred.dd <- expand.grid(age = 50:100, year = 1:15,
                       score = seq(min(dd2$score), max(dd2$score), length.out=100), 
                       id = sample(dd2$id, size = 10, replace = FALSE))

pred <- predict(fit, newdata = pred.dd) ## OK
pred <- predict(fit, newdata = pred.dd, se.fit = TRUE)

pred <- predict(fit, newdata = pred.dd[sample(nrow(pred.dd), size = 1000,
                                            replace = FALSE),], se.fit = TRUE)
