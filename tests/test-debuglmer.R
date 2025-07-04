set.seed(123) # 
n_subjects <- 20
n_times <- 5


sim_data <- data.frame(
  Subject = factor(rep(1:n_subjects, each = n_times)),
  Time = rep(1:n_times, times = n_subjects)
)

rand_effects <- data.frame(
  ID = 1:n_subjects,
  Intercept_re = rnorm(n_subjects, 0, 5),
  Time_re = rnorm(n_subjects, 0, 2)
)


sim_data <- merge(sim_data, rand_effects, by.x = "Subject", by.y = "ID")


sim_data$y <- (10 + sim_data$Intercept_re) +
              (2 + sim_data$Time_re) * sim_data$Time +
              rnorm(nrow(sim_data), 0, 1.5)
          
  
lmer(y ~ Time + us(Time | Subject), data = sim_data, REML = FALSE)

