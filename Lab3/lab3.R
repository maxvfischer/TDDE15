##################################
###########   Lab 3   ############
##################################

## Part 1

sample_emission <- function(z_t) {
  model <- sample(1:3, 1, replace = TRUE)
  
  if (model == 1)  {
    return (rnorm(n = 1, mean = z_t, sd = 1))
  } else if (model == 2) {
    return (rnorm(n = 1, mean = z_t - 1, sd = 1))
  } else if (model == 3) {
    return (rnorm(n = 1, mean = z_t + 1, sd = 1))
  }
}

sample_transition <- function(z_t_1) {
  model <- sample(0:2, 1, replace = TRUE)
  return (rnorm(n = 1, mean = z_t_1 + model, sd = 1))
}

# Number of samples
T = 100
x <- numeric(0)
z <- numeric(0)

# Sample states and observations
for (i in 1:T) {
  if (i == 1) {
    z[i] <- runif(n = 1, min = 0, max = 100) 
  } else {
    z[i] <- sample_transition(z[i-1])  
  }
  x[i] <- sample_emission(z[i])
}

