##################################
###########   Lab 3   ############
##################################

#################### Functions ####################
emission_model <- function(z_t, sd) {
  model <- sample(1:3, 1, replace = TRUE)
  
  if (model == 1)  {
    return (rnorm(n = 1, mean = z_t, sd = sd))
  } else if (model == 2) {
    return (rnorm(n = 1, mean = z_t - 1, sd = sd))
  } else if (model == 3) {
    return (rnorm(n = 1, mean = z_t + 1, sd = sd))
  }
}

transition_model <- function(z_t_1) {
  model <- sample(0:2, 1, replace = TRUE)
  return (rnorm(n = 1, mean = z_t_1 + model, sd = 1))
}

emission_density <- function(x_t, z_t, sd) {
  model1 <- dnorm(x = x_t, mean = z_t, sd = sd)
  model2 <- dnorm(x = x_t, mean = z_t - 1, sd = sd)
  model3 <- dnorm(x = x_t, mean = z_t + 1, sd = sd)
  return ((1/3)*(model1 + model2 + model3))
}

generate_data <- function(T, sd) {
  X_true <- rep(0, T) # Observations
  Z_true <- rep(0, T) # Hiddem states
  
  for (t in 1:T) {
    
    # Generate true hidden states
    if (t == 1) {
      Z_true[t] <- runif(n = 1, min = 0, max = 100) # z_1
    } else {
      Z_true[t] <- transition_model(Z_true[t-1])  
    }
    
    # Generate true observations
    X_true[t] <- emission_model(Z_true[t], sd)
  }
  return (list("Z_true" = Z_true, "X_true" = X_true))
}

#################### Task 1 ####################
## Generate true hidden states and observations
T <- 100 # Number of time steps
sd <- 1
data <- generate_data(T, sd)

## Particle filter algorithm
N <- 100 # Number of particles
Z_particles <- matrix(0, nrow = T, ncol = N) # Matrix with particles during each time step
Z_particles[1, ] <- as.matrix(t(runif(n = N, min = 0, max = 100)))# Initial particles (hypothesis of robot states)
Z_particles_temp <- rep(0, N) # Array that will be filled with new hypothesis
weights <- rep(0, N)

# For each time step
for (t in 1:(T-1)) {
  
  # Generate N particles (hypothesis of hidden state)
  # via transition model p(z_t | z_t-1)
  # Given that we are here (z_t-1), this is the probabilities of
  # where we will be in the next time step
  # 
  # We use a mixed transition model:
  # p(z_t|z_t-1) = 1/3(N(z_t|µ=z_t, σ=1) + N(z_t|µ=z_t + 1, σ=1) + N(z_t|µ=z_t+2, σ=1))
  # Uniform to draw which distribution we're going to draw from
  for (n in 1:N) {

    # Generate new particle (state hypothesis) from old one. p(z_t | z_t_1)
    Z_particles_temp[n] <- transition_model(as.vector(Z_particles[t,n]))
    
    # Generate weight for the new hypothesis.
    # How likely is the observation in this time step, 
    # given the new hypothesis state (robot position)
    weights[n] <- emission_density(data$X_true[t], Z_particles_temp[n], 1)
  }
  
  # Sample 100 particles (hypothesis), with respect to their
  # weight (relative "probability" among the particles)
  Z_particles[t+1,] <- sample(x = Z_particles_temp,
                               size = 100,
                               prob = weights,
                               replace = TRUE)
}

# Plot settings
par(mfrow=c(2, 2))

# First time step
hist(Z_particles[1, ], xlab = 'Particles', main = 't = 1', breaks = 20)
abline(v = mean(Z_particles[1, ]), col = 'red')
print(Z_particles[1, ]) # Particles
print(data$Z_true[1]) # True location

# 30 th time step
hist(Z_particles[30, ], xlab = 'Particles', main = 't = 30', breaks = 20)
abline(v = mean(Z_particles[30, ]), col = 'red')
print(Z_particles[30, ]) # Particles
print(data$Z_true[30]) # True location

# 60th time step
hist(Z_particles[60, ], xlab = 'Particles', main = 't = 60', breaks = 20)
abline(v = mean(Z_particles[60, ]), col = 'red')
print(Z_particles[60, ]) # Particles
print(data$Z_true[60]) # True location

# Last time step
hist(Z_particles[100, ], xlab = 'Particles', main = 't = 100', breaks = 20)
abline(v = mean(Z_particles[100, ]), col = 'red')
print(Z_particles[100, ]) # Particles
print(data$Z_true[100]) # True location

#################### Task 2 ####################
## Generate true hidden states and observations
T <- 100 # Number of time steps
sigmas <- c(5, 50)
data_5 <- generate_data(T, sigmas[1])
data_50 <- generate_data(T, sigmas[2])

## Particle filter algorithm
N <- 100 # No of time steps
Z_5_particles <- matrix(0, nrow = T, ncol = N)
Z_5_particles[1, ]<- runif(n = N, min = 0, max = 100) # Initial particles (hypothesis of robot states)
Z_5_particles_temp <- rep(0, N) # Array that will be filled with new hypothesis

Z_50_particles <- matrix(0, nrow = T, ncol = N)
Z_50_particles[1, ] <- runif(n = N, min = 0, max = 100)
Z_50_particles_temp <- rep(0, N)

weights_5 <- rep(0, N)
weights_50 <- rep(0, N)

# For each time step
for (t in 1:(T-1)) {
  
  # Generate N particles (hypothesis of hidden state)
  # via transition model p(z_t | z_t-1)
  # Given that we are here (z_t-1), this is the probabilities of
  # where we will be in the next time step
  # 
  # We use a mixed transition model:
  # p(z_t|z_t-1) = 1/3(N(z_t|µ=z_t, σ=1) + N(z_t|µ=z_t + 1, σ=1) + N(z_t|µ=z_t+2, σ=1))
  # Uniform to draw which distribution we're going to draw from
  for (n in 1:N) {
    
    # Generate new particle (state hypothesis) from old one. p(z_t | z_t_1)
    Z_5_particles_temp[n] <- transition_model(as.vector(Z_5_particles[t, n]))
    Z_50_particles_temp[n] <- transition_model(as.vector(Z_50_particles[t, n]))
    
    # Generate weight for the new hypothesis.
    # How likely is the observation in this time step, 
    # given the new hypothesis state (robot position)
    weights_5[n] <- emission_density(data_5$X_true[n], Z_5_particles_temp[n], sigmas[1])
    weights_50[n] <- emission_density(data_50$X_true[n], Z_50_particles_temp[n], sigmas[2])
  }
  
  # Sample 100 particles (hypothesis), with respect to their
  # weight (relative "probability" among the particles)
  Z_5_particles[t+1, ] <- sample(x = Z_5_particles_temp,
                            size = 100,
                            prob = weights_5,
                            replace = TRUE)
  
  Z_50_particles[t+1, ] <- sample(x = Z_50_particles_temp,
                            size = 100,
                            prob = weights_50,
                            replace = TRUE)
}

## Sigma = 5

# Plot setup
par(mfrow=c(2, 2))

# First time step
hist(Z_5_particles[1, ], xlab = 'Particles', main = 't = 1', breaks = 20)
abline(v = mean(Z_5_particles[1, ]), col = 'red')
print(Z_5_particles[1, ]) # Particles
print(data_5$Z_true[1]) # True location

# 30 th time step
hist(Z_5_particles[30, ], xlab = 'Particles', main = 't = 30', breaks = 20)
abline(v = mean(Z_5_particles[30, ]), col = 'red')
print(Z_5_particles[30, ]) # Particles
print(data_5$Z_true[30]) # True location

# 60th time step
hist(Z_5_particles[60, ], xlab = 'Particles', main = 't = 60', breaks = 20)
abline(v = mean(Z_5_particles[60, ]), col = 'red')
print(Z_5_particles[60, ]) # Particles
print(data_5$Z_true[60]) # True location

# Last time step
hist(Z_5_particles[100, ], xlab = 'Particles', main = 't = 100', breaks = 20)
abline(v = mean(Z_5_particles[100, ]), col = 'red')
print(Z_5_particles[100, ]) # Particles
print(data_5$Z_true[100 ]) # True location

## Sigma = 50

# Plot setup
par(mfrow=c(2, 2))

# First time step
hist(Z_50_particles[1, ], xlab = 'Particles', main = 't = 1', breaks = 20)
abline(v = mean(Z_50_particles[1, ]), col = 'red')
print(Z_50_particles[1, ]) # Particles
print(data_50$Z_true[1]) # True location

# 30 th time step
hist(Z_50_particles[30, ], xlab = 'Particles', main = 't = 30', breaks = 20)
abline(v = mean(Z_50_particles[30, ]), col = 'red')
print(Z_50_particles[30, ]) # Particles
print(data_50$Z_true[30]) # True location

# 60th time step
hist(Z_50_particles[60, ], xlab = 'Particles', main = 't = 60', breaks = 20)
abline(v = mean(Z_50_particles[60, ]), col = 'red')
print(Z_50_particles[60, ]) # Particles
print(data_50$Z_true[60]) # True location

# Last time step
hist(Z_50_particles[100, ], xlab = 'Particles', main = 't = 100', breaks = 20)
abline(v = mean(Z_50_particles[100, ]), col = 'red')
print(Z_50_particles[100, ]) # Particles
print(data_50$Z_true[100 ]) # True location

#################### Task 3 ####################
## Generate true hidden states and observations
T <- 100 # Number of time steps
sd <- 1
data <- generate_data(T, sd)

## Particle filter algorithm
N <- 100 # Number of particles
Z_particles <- matrix(0, nrow = T, ncol = N)
Z_particles[1, ] <- runif(n = N, min = 0, max = 100) # Initial particles (hypothesis of robot states)
Z_particles_temp <- rep(0, N) # Array that will be filled with new hypothesis
weights <- rep(1, N)

# For each time step
for (t in 1:(T-1)) {
  
  # Generate N particles (hypothesis of hidden state)
  # via transition model p(z_t | z_t-1)
  # Given that we are here (z_t-1), this is the probabilities of
  # where we will be in the next time step
  # 
  # We use a mixed transition model:
  # p(z_t|z_t-1) = 1/3(N(z_t|µ=z_t, σ=1) + N(z_t|µ=z_t + 1, σ=1) + N(z_t|µ=z_t+2, σ=1))
  # Uniform to draw which distribution we're going to draw from
  for (n in 1:N) {
    
    # Generate new particle (state hypothesis) from old one. p(z_t | z_t_1)
    Z_particles_temp[n] <- transition_model(as.vector(Z_particles[t, n]))
  }
  
  # Sample 100 particles (hypothesis), with respect to their
  # weight (relative "probability" among the particles)
  Z_particles[t+1, ] <- sample(x = Z_particles_temp,
                              size = 100,
                              prob = weights,
                              replace = TRUE)
}

# Plot setup
par(mfrow=c(2, 2))

# First time step
plot(Z_particles[1, ], ylab = 'Location', xlab = 'Particles', main = 't = 1')
abline(h = mean(Z_particles[1, ]), col = 'red')
print(Z_particles[1, ]) # Particles
print(data$Z_true[1]) # True location

# 30 th time step
plot(Z_particles[30, ], type='l', ylab = 'Location', xlab = 'Particles', main = 't = 30')
abline(h = mean(Z_particles[30, ]), col = 'red')
print(Z_particles[30, ]) # Particles
print(data$Z_true[30]) # True location

# 60th time step
plot(Z_particles[60, ], type='l', ylab = 'Location', xlab = 'Particles', main = 't = 60')
abline(h = mean(Z_particles[60, ]), col = 'red')
print(Z_particles[60, ]) # Particles
print(data$Z_true[60]) # True location

# Last time step
plot(Z_particles[100, ], type='l', ylab = 'Location', xlab = 'Particles', main = 't = 100')
abline(h = mean(Z_particles[100, ]), col = 'red')
print(Z_particles[100, ]) # Particles
print(data$Z_true[100 ]) # True location
