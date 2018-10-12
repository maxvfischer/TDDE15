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

## Input:
# x: Vector of training inputs.
# y: Vector of training targets/outputs.
# XStar: Vector of inputs where the posterior distribution is evaluated, i.e. X ∗ .
# hyperParam: Vector with two elements, σ_f and l.
# sigmaNoise: Noise standard deviation σ_n .
posteriorGP <- function(x, y, XStar,
                        hyperParam, sigmaNoise) {
  n <- length(x)
  sigma <- hyperParam[1]
  l <- hyperParam[2]
  
  # Covariance matrix of X and X
  # X is of 1-dim 
  K <- squared_exponential(X = x, X_tick = x,
                           sigma = sigma, l = l)
  
  # Covariance matrix of X (train) and X* (test)
  K_star <- squared_exponential(X = x, X_tick = XStar,
                                sigma = sigma, l = l)
  
  # Diagonal matrix with sigma in diagonal
  sigma_noice_identity <- (sigmaNoise^2)*diag(n)
  
  # Cholesky factorization
  # (K + sigma_identity) = L %*% t(L)
  L_upper <- chol(x = (K + sigma_noice_identity))
  L_lower <- t(L_upper)
  
  # Predictive mean
  # Same as on lecture 10 where
  # mean_f = K(XStar, X)%*%inverse([K(X, X) + sigma_identity])*y
  alpha_denom <- solve(a = L_lower, b = y)
  alpha <- solve(a = t(L_lower), b = alpha_denom)
  f_mean <- t(K_star) %*% alpha
  
  # Predictive variance
  f_variance <- solve(a = L_lower, b = K_star)
  
  return (list('pred_mean': f_mean, 'pred_var': f_variance))
}

prior_sigma_f <- 1
prior_l_f <- 0.3
