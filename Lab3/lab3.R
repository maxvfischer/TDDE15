##################################
###########   Lab 3   ############
##################################

## Part 1

emission_model <- function(z_t) {
  model <- sample(1:3, 1, replace = TRUE)
  
  if (model == 1)  {
    return (rnorm(n = 1, mean = z_t, sd = 1))
  } else if (model == 2) {
    return (rnorm(n = 1, mean = z_t - 1, sd = 1))
  } else if (model == 3) {
    return (rnorm(n = 1, mean = z_t + 1, sd = 1))
  }
}

transition_model <- function(z_t_1) {
  model <- sample(0:2, 1, replace = TRUE)
  return (rnorm(n = 1, mean = z_t_1 + model, sd = 1))
}

emission_density <- function(z_t) {
  model1 <- rnorm(n = 1, mean = z_t, sd = 1)
  model2 <- rnorm(n = 1, mean = z_t - 1, sd = 1)
  model3 <- rnorm(n = 1, mean = z_t + 1, sd = 1)
  return ((1/3)*(model1 + model2 + model3))
}

## Generate true hidden states and observations
N <- 100 # Number of states
X_true <- numeric(0) # Observations
Z_true <- numeric(0) # Hiddem states

for (n in 1:N) {
  
  # Generate true hidden states
  if (n == 1) {
    Z_true[n] <- runif(n = 1, min = 0, max = 100) # z_1
  } else {
    Z_true[n] <- transition_model(Z_true[n-1])  
  }
  
  # Generate true observations
  X_true[n] <- emission_model(Z_true[n])
}

## Particle filter algorithm
T <- 100 # No of time steps
Z_particles_t_1 <- runif(n = 100, min = 0, max = 100) # Initial particles (hypothesis of robot states)
Z_particles_t <- numeric(0) # Array that will be filled with new hypothesis
weights <- numeric(0)
# For each time step
for (t in 1:T) {
  
  # Generate N particles (hypothesis of hidden state)
  # via transition model p(z_t | z_t-1)
  # Given that we are here (z_t-1), this is the probabilities of
  # where we will be in the next time step
  # 
  # We use a mixed transition model:
  # p(z_t|z_t-1) = 1/3(N(z_t|µ=z_t, σ=1) + N(z_t|µ=z_t + 1, σ=1) + N(z_t|µ=z_t+2, σ=1))
  for (n in 1:N) {
    
    # Generate new particle (state hypothesis) from old one. p(z_t | z_t_1)
    Z_particles_t[n] <- transition_model(Z_particles_t_1[n])
    
    # Generate weight for the new hypothesis.
    # How likely is the observation in this time step, 
    # given the new hypothesis state (robot position)
    weights[n] <- emission_density(Z_particles_t[n])
  }
  
  # Sample 100 particles (hypothesis), with respect to their
  # weight (relative "probability" among the particles)
  Z_particles_t_1 <- sample(x = Z_particles_t,
                               size = 100,
                               prob = weights)

}






