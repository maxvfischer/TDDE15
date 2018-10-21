########################################
### Lab 4 - Gaussian Process Regression and Classification
### By Max Fischer
### TDDE15 - Advanced Machine Learning
### Linköpings university
########################################
############### Functions ##############

# Name: squared_exponential
# Input:
#   x: Observations
#   x_star: Observations
#   sigma_f: Standard diviation of f
#   l: Smoothness factor
squared_exponential <- function(x, x_star, sigma_f, l) {
  n1 <- length(x)
  n2 <- length(x_star)
  K <- matrix(NA,n1,n2)
  for (i in 1:n2){
    K[,i] <- sigma_f^2*exp(-0.5*( (x-x_star[i])/l)^2 )
  }
  return(K)
}

# Name: nested_squared_exponential
# Input:
#   x: Observations
#   y: Observations
#   sigma_f: Standard diviation of f
#   l: Controls correlation between same day in different years
nested_squared_exponential <- function(sigma_f, l) {
  rval <- squared_exponential <- function(x, y = NULL) {
    n1 <- length(x)
    n2 <- length(y)
    K <- matrix(NA,n1,n2)
    for (i in 1:n2){
      K[,i] <- sigma_f^2*exp(-0.5*( (x-y[i])/l)^2 )
    }
    return(K)
  }
  
  class(rval) <- 'kernel'
  return (rval)
}

# Name: squared_exponential
# Input:
#   x: Observations
#   y: Observations
#   sigma_f: Standard diviation of f
#   l1: Controls periodic part of the correlation
#   l1: Controls correlation between same day in different years
general_periodic_kernel <- function(sigma_f, l1, l2, d) {
  rval <- function (x, y = NULL) {
    r <- abs(x - y)
    return (
      (sigma_f^2)*exp(-2*(sin(pi*r/d)^2)/(l1^2))*exp(-0.5*(r/l2)^2)
    )
  }
  
  class(rval) <- 'kernel'
  return (rval)
}

# Name: posterior_GP
# Input:
#   x: Observations
#   y: Observations
#   x_star: Values to predict posterior mean of f over
#   kernel: Covariance function
#   sigma_n: Standard diviation of the measured data
#   sigma_f: Standard diviation of f
#   l: Controls correlation between same day in different years
posterior_GP <- function(x, y, x_star, kernel, sigma_n, sigma_f, l) {
  
  # Number of observations
  n <- length(x)
  
  # Calculate the covariance matricies:
  # k(X, X), k(X, X*), k(X*, X*)
  K_x_x <- squared_exponential(x = x,
                               x_star = x,
                               sigma_f = sigma_f,
                               l = l)
  K_x_xstar <- squared_exponential(x = x,
                                   x_star = x_star,
                                   sigma_f = sigma_f,
                                   l = l)
  K_xstar_xstar <- squared_exponential(x = x_star,
                                       x_star = x_star,
                                       sigma_f = sigma_f,
                                       l = l)
  
  # Compute the Choleski factorization of 
  # k(X, X) + sigma_n^2
  # (covariance matrix of y)
  #
  # As chol returns the upper triangular part and
  # we need the lower, we transpose it
  L_upper <- chol(K_x_x + (sigma_n^2)*diag(n))
  L_lower <- t(L_upper)
  
  # Compute alpha, used to compute the
  # posterior mean of f
  alpha_b <- solve(a = L_lower,
                             b = y)
  alpha <- solve(a = t(L_lower),
                 b = alpha_b)
  
  # Compute posterior mean of f
  posterior_mean_f <- t(K_x_xstar) %*% alpha
  
  # Compute posterior covariance matrix of f
  v <- solve(a = L_lower,
             b = K_x_xstar)
  posterior_covariance_matrix_f <- K_xstar_xstar - t(v) %*% v
  
  # As we only want the variance of f, we extract the
  # diagonal in the covariance matrix of f
  posterior_variance_f <- diag(posterior_covariance_matrix_f)
  
  return (list(mean = posterior_mean_f, variance = posterior_variance_f))
}

# Name: posterior_GP
# Input:
#   mean: Mean to be plotted along the y-axis
#   interval: Values to be plotted along the x-axis
#   variance: Variance of the mean
#   observations: Measurements done
visualize_GP <- function(mean, interval, variance = NULL, observations) {
  
  if (!is.null(variance)) {
    # Compute confidence interval
    CI <- data.frame(upper = mean + 1.96*sqrt(variance),
                     lower = mean - 1.96*sqrt(variance))
    
    # Compute visually nice y-lim
    ylim <- c(min(CI$lower) - 1,
              max(CI$upper) + 1)
    
    plot(x = interval,
         y = mean,
         type = 'l',
         col = 'red',
         ylab = 'Posterior mean',
         xlab = 'Interval',
         ylim = ylim)
    
    # Draw confidence interval on plot
    polygon(x = c(rev(interval), interval),
            y = c(rev(CI$upper), CI$lower),
            col = rgb(0, 0, 0, 0.3)) 
    
    # Add observations as points
    points(x = observations$x,
           y = observations$y,
           col = 'blue',
           pch = 16)
    
    # Add legend to top right corner
    legend('topright',
           legend = c('Mean of f', '95 % CI', 'Observations'),
           col = c('red', rgb(0, 0, 0, 0.3), 'blue'),
           lty = c(1, NA, NA),
           pch = c(NA, 15, 16))
    
  } else {
    
    # Compute visually nice y-lim
    ylim <- c(min(observations$y) - 1,
              max(observations$y) + 1)
    
    plot(x = interval,
         y = mean,
         type = 'l',
         col = 'red',
         ylab = 'Posterior mean',
         xlab = 'Interval',
         ylim = ylim)
    
    # Add observations as points
    points(x = observations$x,
           y = observations$y,
           col = 'blue',
           pch = 16)
    
    # Add legend to top right corner
    legend('topright',
           legend = c('Mean of f', 'Observations'),
           col = c('red', 'blue'),
           lty = c(1, NA),
           pch = c(NA, 16))
  }
  
}

################ Setup ################
par(mfrow = c(1, 1))

############## Task 2.1.1 #############
# See functions:
# posterior_GP
# squared_exponential

############## Task 2.1.2 #############

# Set hyperparameters
sigma_f <- 1
l <- 0.3

# Set standard deviation of measurement
sigma_n <- 0.1

# Measurements
observations <- data.frame(x = c(0.4),
                           y = c(0.719))

# Set interval to get posterior from
interval <- seq(from = -1,
                to = 1,
                length.out = 100)

# Get posterior mean and variance of f
posterior_f <- posterior_GP(x = observations$x,
                            y = observations$y,
                            x_star = interval,
                            sigma_f = sigma_f,
                            l = l,
                            sigma_n = sigma_n)

# Visalize posterior mean and CI of posterior mean
visualize_GP(mean = posterior_f$mean,
             interval = interval,
             variance = posterior_f$variance,
             observations = observations)

############## Task 2.1.3 #############

# Add new observation
observations <- data.frame(x = c(0.4, -0.6),
                           y = c(0.719, 0.044))

# Get posterior mean and variance of f
posterior_f <- posterior_GP(x = observations$x,
                            y = observations$y,
                            x_star = interval,
                            sigma_f = sigma_f,
                            l = l,
                            sigma_n = sigma_n)

# Visalize posterior mean and CI of posterior mean
visualize_GP(mean = posterior_f$mean,
             interval = interval,
             variance = posterior_f$variance,
             observations = observations)

############## Task 2.1.4 #############

# Add new observation
observations <- data.frame(x = c(-1.0, -0.6, -0.2, 0.4, 0.8),
                           y = c(0.768, -0.044, -0.940, 0.719, -0.664))

# Get posterior mean and variance of f
posterior_f <- posterior_GP(x = observations$x,
                            y = observations$y,
                            x_star = interval,
                            sigma_f = sigma_f,
                            l = l,
                            sigma_n = sigma_n)

# Visalize posterior mean and CI of posterior mean
visualize_GP(mean = posterior_f$mean,
             interval = interval,
             variance = posterior_f$variance,
             observations = observations)

############## Task 2.1.5 #############

# Update hyperparameters
sigma_f <- 1
l <- 1

# Get posterior mean and variance of f
posterior_f <- posterior_GP(x = observations$x,
                            y = observations$y,
                            x_star = interval,
                            sigma_f = sigma_f,
                            l = l,
                            sigma_n = sigma_n)

# Visalize posterior mean and CI of posterior mean
visualize_GP(mean = posterior_f$mean,
             interval = interval,
             variance = posterior_f$variance,
             observations = observations)

################ Setup ################

# Install packages if not already installed
if (!require(kernlab)) {
  install.packages('kernlab')
}

# Import packages
library(kernlab)

# Import data
temp_tullinge <- read.csv("https://raw.githubusercontent.com/STIMALiU/AdvMLCourse/master/GaussianProcess/Code/TempTullinge.csv", header=TRUE, sep=";")

# Create variables
time <- seq(from = 1,
            to = 365*6,
            by = 5)
day <- rep(
  x = seq(from = 1,
          to = 365,
          by = 5),
  times = 6
  )

# Extract temperatures for every fifth day
temp_time <- temp_tullinge$temp[time]

############## Task 2.2.1 #############

# Create data variables
x <- 1
x_star <- 2

# Instantiate kernel
kernel <- nested_squared_exponential(sigma_f = 1, l = 0.3)

# Evaluate kernel on x = 1, x_star = 2
variance <- kernel(x = x,
                     y = x_star)

# Create data variables
x <- c(1, 3, 4)
x_star <- c(2, 3, 4)

covariance_matrix <- kernelMatrix(x = x,
                                  y = x_star,
                                  kernel = kernel)

############## Task 2.2.2 #############

# Generate standard deviation of measurements 
# by computing the standard deviation of the
# residuals from a linear quadratic regression,
# fitted on: temp ~ time + time^2
fit <- lm(temp_time ~ time + time^2)
sigma_n <- sd(fit$residuals)

# Set hyperparameters
sigma_f <- 20
l <- 0.2

# Fit Gaussian Process regression
GP.fit <- gausspr(x = time,
                  y = temp_time,
                  kernel = nested_squared_exponential,
                  kpar = list(sigma_f = sigma_f, l = l),
                  var = sigma_n^2)

# Predict via fitted Gaussian Process regression
GP.mean_task2 <- predict(GP.fit, time)

# Visualize prediction
visualize_GP(mean = GP.mean_task2,
             interval = time,
             observations = data.frame(x = time, y = temp_time))

############## Task 2.2.3 #############

# Instantiate kernel
kernel <- nested_squared_exponential(sigma_f = 1, l = 0.3)

# Use own implemented function to generate variances
# of the mean posterior
GP.pred <- posterior_GP(x = scale(time),
                        y = scale(temp_time),
                        x_star = scale(seq(from = 1, to = 365*6, by = 1)),
                        kernel = kernel,
                        sigma_n = sigma_n,
                        sigma_f = sigma_f,
                        l = l)

# Visualize prediction
visualize_GP(mean = GP.mean_task2,
             interval = time,
             observations = data.frame(x = time, y = temp_time),
             variance = GP.pred$variance[time])

############## Task 2.2.4 #############

# Set hyperparameters
sigma_f <- 20
l <- 0.2

# Fit Gaussian Process regression to day variable
# Day variable treats same day of the year as the
# same variable.
GP.fit <- gausspr(x = day,
                  y = temp_time,
                  kernel = nested_squared_exponential,
                  kpar = list(sigma_f = sigma_f, l = l),
                  var = sigma_n^2)

# Predict via fitted Gaussian Process regression
GP.mean_task4 <- predict(GP.fit, day)

# Visualize previous prediction
visualize_GP(mean = GP.mean_task2,
             interval = time,
             observations = data.frame(x = time, y = temp_time),
             variance = GP.pred$variance[time])

# Add prediction done with day variable to the vizualization
lines(x = time, y = GP.mean_task4, col = 'green')

############## Task 2.2.5 #############
# See function: general_periodic_kernel

# Set hyperparameters
sigma_f <- 20
l1 <- 1
l2 <- 20
d <- 365/sd(time)

# Fit Gaussian Process regression to dtime variable
# and with the general periodic kernel
GP.fit <- gausspr(x = time,
                  y = temp_time,
                  kernel = general_periodic_kernel,
                  kpar = list(sigma_f = sigma_f, l1 = l1, l2 = l2, d = d),
                  var = sigma_n^2)

# Predict via fitted Gaussian Process regression
GP.mean_task5 <- predict(GP.fit, time)

# Visualize previous prediction
visualize_GP(mean = GP.mean_task2,
             interval = time,
             observations = data.frame(x = time, y = temp_time),
             variance = GP.pred$variance[time])

# Add prediction done with day variable to the vizualization
lines(x = time, y = GP.mean_task4, col = 'green')

# Add prediction done with time variable and general periodic kernel
lines(x = time, y =GP.mean_task5, col = 'black')

############## Task 2.3 #############

################ Setup ################

# Import data
data <- read.csv("https://raw.githubusercontent.com/STIMALiU/AdvMLCourse/master/GaussianProcess/Code/banknoteFraud.csv", header=FALSE, sep=",")

# Set names to columns
names(data) <- c("varWave","skewWave","kurtWave","entropyWave","fraud")
data[,5] <- as.factor(data[,5])

# Set seed
set.seed(12345)

# Split data into train and test set
train_indices <- sample(1:dim(data)[1], size = 1000,
                        replace = FALSE)
data.train <- data[train_indices, ]
data.test <- data[-train_indices, ]

############## Task 2.3.1 #############

# Fit Gaussian Process classifier
GP.fit <- gausspr(fraud ~ varWave + skewWave, data = data.train)

# Create grid
x1 <- seq(from = min(data.train$varWave),
          to = max(data.train$varWave),
          length = 100)
x2 <- seq(from = min(data.train$skewWave),
          to = max(data.train$skewWave),
          length = 100)
gridPoints <- meshgrid(x1, x2)
gridPoints <- cbind(c(gridPoints$x), c(gridPoints$y))
gridPoints <- data.frame(gridPoints)
names(gridPoints) <- c("varWave", "skewWave")

# Predict via fitted Gaussian Process classification
# on grid
GP.pred_grid <- predict(GP.fit, gridPoints, type="probabilities")

# Get indices of fraud
fraud_indices <- which(data.train$fraud == 1)

# Render contour of varWave and skewWave
contour(x = x1,
        y = x2,
        z = matrix(GP.pred_grid[,2], 100, byrow = TRUE), 
        20,
        xlab = "varWave", ylab = "skewWave", main = 'Prob(fraud)')

# Add data points of fraud/non-fraud by varWave and skewWave
points(x = data.train$varWave[fraud_indices],
       y = data.train$skewWave[fraud_indices],
       col = "blue")
points(x = data.train$varWave[-fraud_indices],
       y = data.train$skewWave[-fraud_indices],
       col = "red")

# Predict via fitted Gaussian Process classification
# on training data
GP.pred_train <- predict(GP.fit, data.train)

# Compute confusion matrix
confusion_matrix_train <- table(GP.pred_train, data.train$fraud)

# Compute accuracy
accuracy_train <- sum(diag(confusion_matrix_train))/sum(confusion_matrix_train)

############## Task 2.3.2 #############

# Predict via fitted Gaussian Process classification
# on testing data
GP.pred_test <- predict(GP.fit, data.test)

# Compute confusion matrix
confusion_matrix_test <- table(GP.pred_test, data.test$fraud)

# Compute accuracy
accuracy_test <- sum(diag(confusion_matrix_train))/sum(confusion_matrix_train)

############## Task 2.3.3 #############

# Fit Gaussian Process classifier on all variables
GP.fit_all_var <- gausspr(fraud ~., data = data.train)

# Predict via fitted Gaussian Process classification
# on testing data
GP.pred_all_var <- predict(GP.fit_all_var, data.test)

# Compute confusion matrix
confusion_matrix_test_all_var <- table(GP.pred_all_var, data.test$fraud)

# Compute accuracy
accuracy_test_all_var <- sum(diag(confusion_matrix_test_all_var))/sum(confusion_matrix_test_all_var)
