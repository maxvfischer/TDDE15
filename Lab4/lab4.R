##################################
###########   Lab 4   ############
##################################

#################### Functions ####################
squared_exponential <- function(X, X_tick, sigma, l) {
  n1 <- length(X)
  n2 <- length(X_tick)
  K <- matrix(NA,n1,n2)
  for (i in 1:n2){
    K[,i] <- sigma^2*exp(-0.5*( (X-X_tick[i])/l)^2 )
  }
  return(K)
}


#################### Task 1 ####################
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
  L_upper <- chol(x = (K.X_X + diagonal_variance_noise))
  L_lower <- t(L_upper)

  # Predictive mean
  # Same as on lecture 10 where
  # mean_f = K(XStar, X)%*%inverse([K(X, X) + sigma_identity])*y
  alpha_denom <- solve(a = L_lower, b = y)
  alpha <- solve(a = t(L_lower), b = alpha_denom)
  f_mean <- t(K.X_XStar) %*% alpha
  
  # Predictive variance
  v <- solve(a = L_lower, b = K.X_XStar)
  f_variance <- K.XStar_XStar - t(v) %*% v
  return (list('pred_mean' = f_mean, 'pred_var' = f_variance))
}

#################### Task 2 ####################
sigma_f <- 1
l <- 0.3
hyperParam <- c(sigma_f, l)
x <- 0.4
y <- 0.719
x_star <- seq(from = -1, to = 1, by = 0.01)
sigma_n <- 0.1

GP_posterior <- posteriorGP(x = x, y = y, XStar = x_star,
                             hyperParam = hyperParam,
                             sigmaNoise = sigma_n)

# Calculate lower and upper 95 % pointwise confidence interval
GP_CI <- data.frame(upper = as.vector(GP_posterior$pred_mean) + 1.96 * as.vector(sqrt(GP_posterior$pred_var)),
                    lower = as.vector(GP_posterior$pred_mean) - 1.96 * as.vector(sqrt(GP_posterior$pred_var)))
plot(x = x_star, y = GP_posterior$pred_mean, type = 'l')
lines(x = x_star, y = GP_CI$upper)
lines(x = x_star, y = GP_CI$lower)
######################





