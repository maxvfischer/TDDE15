##################################
###########   Lab 4   ############
##################################

#################### Functions ####################

## Kernel (covariate)
squared_exponential <- function(X, X_tick, sigma, l) {
  n1 <- length(X)
  n2 <- length(X_tick)
  K <- matrix(NA,n1,n2)
  for (i in 1:n2){
    K[,i] <- sigma^2*exp(-0.5*( (X-X_tick[i])/l)^2 )
  }
  return(K)
}

# Plot 
plotGP <- function(x = NULL, y = NULL, x_star, GP_mean, GP_var) {
  
  # Calculate lower and upper 95 % pointwise confidence interval
  GP_CI <- data.frame(upper = as.vector(GP_mean) + 1.96 * sqrt(GP_var),
                    lower = as.vector(GP_mean) - 1.96 * sqrt(GP_var))
  
  # Upper and lower y-limits based on max/min mean +- std (pointwise)
  y_lim <- c(min(GP_CI$lower) - 0.5, max(GP_CI$upper) + 0.5)
  
  # Generate initial plot with labels
  plot(x = x_star, 
       y = GP_mean, 
       type = 'l', 
       ylim = y_lim,
       xlab = 'x-values',
       ylab = 'Posterior of f')
  
  # Pointwise confidence interval (grayed)
  polygon(x = c(rev(x_star), x_star), 
          y = c(rev(GP_CI$upper), GP_CI$lower), 
          col = 'grey80', border = NA,
          ylim = y_lim,
          xlab = 'x-values',
          ylab = 'Posterior of f')
  
  # Lines: Mean, upper CI and lower CI
  lines(x = x_star, y = GP_mean, col = 'red')
  lines(x = x_star, y = GP_CI$upper, lty = 'dashed', col = 'red')
  lines(x = x_star, y = GP_CI$lower, lty = 'dashed', col = 'red')
  
  # Observations
  if (!is.null(x) && !is.null(y)) {
    points(x = x, y = y, type = 'p')
  }
}

#################### Task 2.1.1 ####################
## Input:
# x: Vector of training inputs.
# y: Vector of training targets/outputs.
# XStar: Vector of inputs where the posterior distribution is evaluated, i.e. X ∗ .
# hyperParam: Vector with two elements, σ_f and l.
# sigmaNoise: Noise standard deviation σ_n .
posteriorGP <- function(x, y, XStar,
                        hyperParam, sigmaNoise) {
  n <- length(x)
  sigma_f <- hyperParam[1]
  l <- hyperParam[2]
  
  # Covariance matrices
  # k(X, X)
  # k(X, X*)
  # k(X*, X*)
  K.X_X = squared_exponential(X = x, X_tick = x,
                           sigma = sigma_f, l = l)
  K.X_XStar <- squared_exponential(X = x, X_tick = XStar,
                                               sigma = sigma_f, l = l)
  K.XStar_XStar <- squared_exponential(X = XStar, X_tick = XStar,
                                                   sigma = sigma_f, l = l)
  
  # Diagonal matrix of sigma noise
  diagonal_variance_noise <- (sigmaNoise^2) * diag(length(x))

  # Cholesky factorization
  # (K + sigma_identity) = L %*% t(L)
  L_upper <- chol(K.X_X + diagonal_variance_noise)
  L_lower <- t(L_upper)

  # Predictive mean
  # Same as on lecture 10 where
  # mean_f = K(XStar, X)%*%inverse([K(X, X) + sigma_identity])*y
  alpha_denom <- solve(a = L_lower, b = y)
  alpha <- solve(a = t(L_lower), b = alpha_denom)
  f_mean <- t(K.X_XStar) %*% alpha
  
  # Predictive variance
  v <- solve(a = L_lower, b = K.X_XStar)
  f_variance <- diag(K.XStar_XStar - t(v) %*% v)
  return (list('pred_mean' = f_mean, 'pred_var' = f_variance))
}

#################### Task 2.1.2 ####################
sigma_f <- 1
l <- 0.3
hyperParam <- c(sigma_f, l)
x <- 0.4
y <- 0.719
x_star <- seq(from = -1, to = 1, by = 0.01)
sigma_n <- 0.1

# x: train
# y: train
# XStar: test
GP_posterior <- posteriorGP(x = x, y = y, XStar = x_star,
                             hyperParam = hyperParam,
                             sigmaNoise = sigma_n)

GP_CI <- data.frame(upper = as.vector(GP_posterior$pred_mean) + 1.96 * sqrt(GP_posterior$pred_var),
                    lower = as.vector(GP_posterior$pred_mean) - 1.96 * sqrt(GP_posterior$pred_var))
plotGP(x = x, y = y, x_star = x_star, 
       GP_mean = GP_posterior$pred_mean,
       GP_var = GP_posterior$pred_var)

#################### Task 2.1.3 ####################
sigma_f <- 1
l <- 0.3
hyperParam <- c(sigma_f, l)
x <- c(0.4, -0.6)
y <- c(0.719, -0.044)

x_star <- seq(from = -1, to = 1, by = 0.01)
sigma_n <- 0.1

GP_posterior <- posteriorGP(x = x, y = y, XStar = x_star,
                            hyperParam = hyperParam,
                            sigmaNoise = sigma_n)

plotGP(x = x, y = y, x_star = x_star,
       GP_mean = GP_posterior$pred_mean,
       GP_var = GP_posterior$pred_var)

#################### Task 2.1.4 ####################
sigma_f <- 1
l <- 0.3
hyperParam <- c(sigma_f, l)
x <- c(-1.0, -0.6, -0.2, 0.4, 0.8)
y <- c(0.768, -0.044, -0.940, 0.719, -0.664)

x_star <- seq(-1, 1, 0.01)
sigma_n <- 0.1

GP_posterior <- posteriorGP(x = x, y = y, XStar = x_star,
                            hyperParam = hyperParam,
                            sigmaNoise = sigma_n)

plotGP(x = x, y = y, x_star = x_star,
       GP_mean = GP_posterior$pred_mean,
       GP_var = GP_posterior$pred_var)

#################### Task 2.1.5 ####################
sigma_f <- 1
l <- 5
hyperParam <- c(sigma_f, l)
x <- c(-1.0, -0.6, -0.2, 0.4, 0.8)
y <- c(0.768, -0.044, -0.940, 0.719, -0.664)

x_star <- seq(-1, 1, 0.01)
sigma_n <- 0.1

GP_posterior <- posteriorGP(x = x, y = y, XStar = x_star,
                            hyperParam = hyperParam,
                            sigmaNoise = sigma_n)

plotGP(x = x, y = y, x_star = x_star,
       GP_mean = GP_posterior$pred_mean,
       GP_var = GP_posterior$pred_var)

###################### Task 2 ######################

# Install and import package
if (!require(kernlab)) {
  install.packages('kernlab')
}

library(kernlab)

# Import dataset
temp_tullinge <- read.csv("https://raw.githubusercontent.com/STIMALiU/AdvMLCourse/master/GaussianProcess/Code/TempTullinge.csv", header=TRUE, sep=";")

# Variables
time <- seq(from = 1, to = 365*6, by = 5) # Every fifth day since start
year_days <- seq(from = 1, to = 365, by = 5) # Everu fifth day of the year
day <- c(year_days, year_days, year_days, year_days, year_days, year_days) # year_days concatinated 6 times
temp_time <- temp_tullinge$temp[time] # Only the temperatures of every fifth day

#################### Task 2.2.1 ####################
#################### Functions #####################

# Squared exponential kernel wrapped in a function
squared_exponential_nested <- function(sigma_f = 1, l = 1) {
  rval <- function(x, y = NULL) {
    n1 <- length(x)
    n2 <- length(y)
    K <- matrix(NA, n1, n2)
    r <- abs(x - y)
    for (i in 1:n2) {
      K[, i] <- (sigma_f^2)*exp(-0.5*(r/l))
    }
    return (K)
  }
  class(rval) <- 'kernel'
  return(rval)
}

sigma_f <- 1
l <- 0.3
x <- c(1, 3, 4)
x_star <- c(2, 3, 4)

sqrd_exp <- squared_exponential_nested(sigma_f = sigma_f, l = l)
eval <- sqrd_exp(x = 1, y = 2)
K <- kernelMatrix(kernel = sqrd_exp, x = x, y = x_star)

#################### Task 2.2.2 ####################
# temp = f(time) + ε
# ε ~ N(0, σ_n^2)
# f ~ GP(0, k(time, time'))
set.seed(12345)
sigma_f <- 20
l <- 0.2

# Estimate noise variable (sigma_n)
fit <- lm(temp_time ~ time + time^2)
sigma_n <- sd(fit$residuals)

# Fit Gaussian process
GPfit_time <- gausspr(x = time,
                 y = temp_time, 
                 kernel = squared_exponential_nested, 
                 kpar = list(sigma_f = sigma_f, l = l),
                 var = sigma_n^2)

# Compute posterior mean of every training data
GPpred_time <- predict(GPfit_time, time)

# Plot
plot(x = time, y = temp_time, ylab = 'Temperature', xlab = 'Number of days')
lines(x = time, y = GPpred_time, col = 'red')
legend("topright", 
       legend = c("Data", "Posterior mean"), 
       col = c("black", "red"), 
       pch = c(1, NA),
       lty = c(NA, 1))

#################### Task 2.2.3 ####################

# Scale both time, GPpred and XStar, or just GPpred?
posterior_f <- posteriorGP(x = (time),
                           y = scale(GPpred_time),
                           XStar = (seq(1, 365*6, 1)),
                           hyperParam = c(sigma_f, l),
                           sigmaNoise = sigma_n)

# Plot
plotGP(x_star = time,
       GP_mean = GPpred_time,
       GP_var = posterior_f$pred_var[time])
points(x = time, y = temp_time)
legend('topright',
       legend = c("Data", "Posterior mean", "95% CI"),
       col = c("black", "red", "grey80"),
       pch = c(1, NA, 15),
       lty = c(NA, 1, NA))

#################### Task 2.2.4 ####################
# temp = f(day) + ε
# ε ~ N(0, σ_n^2)
# f ~ GP(0, k(day, day'))
# temp_time will work for day aswell, as it's the temp of every 6th day
set.seed(12345)
sigma_f <- 20
l <- 0.2

# Estimate noise variable (sigma_n)
fit <- lm(temp_time ~ day + day^2)
sigma_n <- sd(fit$residuals)

# Fit Gaussian process
GPfit_day <- gausspr(x = day,
                 y = temp_time,
                 kernel = squared_exponential_nested,
                 kpar = list(sigma_f = sigma_f, l = l),
                 var = sigma_n^2)

# Predict f mean with fitted Gaussian process
GPpred_day <- predict(GPfit_day, day)

# Plot data-points, time and day
plot(x = time, y = temp_time, ylab = 'Temperature', xlab = 'Days', ylim = c(min(temp_time), max(temp_time) + 10))
lines(x = time, y = GPpred_time, col = 'red')
lines(x = time, y = GPpred_day, col = 'blue')
legend('topright',
       legend = c('Data', 'Time', 'Day'),
       col = c('black', 'red', 'blue'),
       pch = c(1, NA, NA),
       lty = c(NA, 1, 1))

#################### Task 2.2.5 ####################

#################### Functions #####################

# Generalized pereodic kernel
generalized_periodic_kernel <- function(sigma_f, l1, l2, d) {
  rval <- function(x, x_star) {
    return ((sigma_f^2) * exp( -(2 * sin(pi * abs(x - x_star) / d)^2) / l1^2) * exp(-0.5 * (abs(x - x_star)^2) / l2^2))
  }
  class(rval) <- 'kernel'
  return(rval)
}
set.seed(12345)
sigma_f <- 20
l1 <- 1
l2 <- 10
d <- 365/sd(time)

# Calculate noise variable (sigma_n)
fit <- lm(temp_time ~ time)
sigma_n <- sd(fit$residuals)

# Fit Gaussian process
GPfit_periodic <- gausspr(x = time,
                 y = temp_time,
                 kernel = generalized_periodic_kernel,
                 kpar = list(sigma_f = sigma_f, l1 = l1, l2 = l2, d = d),
                 var = sigma_n)

# Predict mean of f with Gaussian process
GPpred_periodic <- predict(GPfit_periodic, time)

# Plot
plot(x = time, y = temp_time, ylab = 'Temperature', xlab = 'Days', ylim = c(min(temp_time), max(temp_time) + 10))
lines(x = time, y = GPpred_time, col = 'red')
lines(x = time, y = GPpred_day, col = 'blue')
lines(x = time, y = GPpred_periodic, col = 'green')
legend('topright',
       legend = c('Data', 'Time', 'Day', 'Periodic'),
       col = c('black', 'red', 'blue', 'green'),
       pch = c(1, NA, NA, NA),
       lty = c(NA, 1, 1, 1))

#################### Task 2.3 ####################
# Install necessary packages
if (!require('AtmRay')) {
  install.packages('AtmRay')
}
library(AtmRay)

# Import data
data <- read.csv("https://raw.githubusercontent.com/STIMALiU/AdvMLCourse/master/GaussianProcess/Code/banknoteFraud.csv", header=FALSE, sep=",")
names(data) <- c("varWave","skewWave","kurtWave","entropyWave","fraud")
data[,5] <- as.factor(data[,5])

# Sample 1000 training points from data
set.seed(12345)
train_indices <- sample(x = 1:dim(data)[1], size = 1000, replace = FALSE)
data.train <- data[train_indices, ]
data.test <- data[-train_indices, ]

#################### Task 2.3.1 ####################

# Fit a Gaussian Process classification on training data
# Fit fraud by varWave and skewWave
GPfit_classification <- gausspr(fraud ~ varWave + skewWave,
                                data = data.train)

# Grid values
grid.varWave <- seq(from = min(data.train$varWave),
                    to = max(data.train$varWave),
                    length.out = 100)
grid.skewWave <- seq(from = min(data.train$skewWave),
                     to = max(data.train$skewWave),
                     length.out = 100)
# x: Each value in y has a corresponding row in x with all x values,
#    thus, each row in x will be the same
# y: Each value in x has a corresponding column in y with all y values,
#    thus, each column in y will be the same
gridPoints <- meshgrid(grid.varWave, grid.skewWave)
# All combinations of varWave and skewWave
# For each value in x (varWave), there is combinations with all y values (skewWave)
gridPoints <- cbind(c(gridPoints$x), c(gridPoints$y))
# Data frame the above
gridPoints <- data.frame(gridPoints) 
names(gridPoints) <- c('varWave', 'skewWave')

# Prediction by Gaussian process
# Probability of each combination of varWave and skewWave
# in gridPoints being 0 or 1 (not fraud or fraud)
GPgrid_probs <- predict(GPfit_classification, 
                                 gridPoints,
                                 type = 'probabilities')

# Indices of fraud/non-fraud observations
fraud_indices <- which(data.train$fraud %in% c(1))
non_fraud_indices <- which(data.train$fraud %in% c(0))

# Plot contours
# Create grid by varWave and skewWave
# GPgrid_probs[, 1] contains the probabilities of fraud for each
# combination of varWave and skewWave in the grid
contour(x = grid.varWave, 
        y = grid.skewWave,
        z = matrix(GPgrid_probs[, 2], 100, byrow = TRUE),
        xlab = 'varWave',
        ylab = 'skewWave',
        main = 'P(fraud)')
# Plot fraud points
points(x = data.train$varWave[fraud_indices], 
       y = data.train$skewWave[fraud_indices],
       col = 'blue')
# Plot non-fraud points
points(x = data.train$varWave[non_fraud_indices],
       y = data.train$skewWave[non_fraud_indices],
       col = 'red')
legend('topright',
       legend = c('Fraud', 'Not fraud'),
       col = c('blue', 'red'),
       pch = c(1, 1))

# Confusion matrix and accuracy on train data
GPpred_classification_train <- predict(GPfit_classification, data.train)
confusion_matrix_train <- table(GPpred_classification_train, data.train$fraud)
accuracy_train <- sum(diag(confusion_matrix_train))/sum(confusion_matrix_train)

#################### Task 2.3.2 ####################
# Confusion matrix and accuracy on test data
GPpred_classification_test <- predict(GPfit_classification, data.test)
confusion_matrix_test <- table(GPpred_classification_test, data.test$fraud)
accuracy <- sum(diag(confusion_matrix_test))/sum(confusion_matrix_test)

#################### Task 2.3.3 ####################
# Fit the Gaussian process classifier on all variables
GPfit_classification_all_var <- gausspr(fraud ~., data = data.train)

# Predict by test data
GPpred_classification_test <- predict(GPfit_classification_all_var, data.test)

# Confusion matrix and accuracy
confusion_matrix_all_var <- table(GPpred_classification_test, data.test$fraud)
accuracy <- sum(diag(confusion_matrix_all_var))/sum(confusion_matrix_all_var)
