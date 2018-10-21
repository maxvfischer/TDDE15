########################################
### Lab 3 - State-Space Models
### By Max Fischer
### TDDE15 - Advanced Machine Learning
### Linköpings university
########################################
############### Functions ##############

# Name: kalman
# Input:
#   T: Number of time steps
#   A: How the mean of z_t will be affected by z_(t-1)
#      Used in transition model
#   B: How much the control variable will affect
#      z_t. Used in transition model
#   C: Scale the mean (z_t) in the emission model
#   Q: Covariate matrix in the emission model
#   R: Covariate matrix in the transition model
#   obs: All observations
kalman <- function(T, A, B, C, Q, R, obs, u, init_var = 1) {
  my <- NULL
  sigma <- NULL
  my_unweighted <- NULL
  sigma_unweighted <- NULL
  kalman_gain <- NULL
  
  # Our best guess is that my_1 is our first observation
  my[1] <- obs[1]
  
  # We don't know what sigma_1 is, so lets just set it to
  # a random number. We chose 1
  sigma[1] <- init_var
  
  for (t in 2:T) {
    # Calculate the unweighted prediction of the mean
    my_unweighted[t] <- A[t]*my[t - 1] + B[t]*u[t]
    
    # Calculate the unweighted prediction of the covariate matrix
    sigma_unweighted[t] <- A[t]*sigma[t - 1]*t(A[t]) + R[t]
    
    # Kalman gain.
    # Used to weight between our unweighted prediction and the
    # observation
    kalman_gain[t] <- sigma_unweighted[t]*t(C[t]) * inv(C[t]*sigma_unweighted[t]*t(C[t] + Q[t]))
    
    # Calculate the weighted mean, thus our prediction of the
    # hidden state
    my[t] <- my_unweighted[t] + kalman_gain[t]*(obs[t] - C[t]*my_unweighted[t])
    
    # Calculate the weighted covariance matrix, thus our prediction
    # of the predition error
    sigma[t] <- (I - kalman_gain[t]*C[t])*sigma_unweighted[t]
  }
  
  return (list(my = my, sigma = sigma))
}

# Name: sample_emission_model
# Input:
#   z_t_1: Hidden state in previous time step
#   sd: Standard deviation of the models
#       in the mixed model
sample_transition_model <- function(z_t_1, sd = 1) {
  
  # Sample from which model in the mixed model
  # we're sampling from
  model <- sample(x = 0:2, size = 1)
  
  # Return the sample
  return (rnorm(n = 1,
                mean = z_t_1 + model,
                sd = sd))
}

# Name: sample_emission_model
# Input:
#   z_t: Hidden state in the current time step
#   sd: Standard deviation of the models
#       in the mixed model
sample_emission_model <- function(z_t, sd = 1) {
  
  # Sample from which model in the mixed model
  # we're sampling from
  model <- sample(x = -1:1, size = 1)
  
  # Return the sample
  return (rnorm(n = 1,
                mean = z_t + model,
                sd = sd))
}

density_emission_model <- function(x_t, z_t, sd = 1) {
  models <- sapply(-1:1, function(x) {
    dnorm(x = x_t,
          mean = z_t + x,
          sd = sd)
  })
  
  return (sum(models)/3)
}

# Name: sample_data
# Input:
#   T: Number of samples (time steps)
sample_data <- function(T, sd_emission = 1, sd_transition = 1) {
  
  state_variable <- NULL
  observed_variable <- NULL
  
  # Inizialize z_1
  # Defined as Uniform(0, 100)
  state_variable[1] <- runif(n = 1,
                             min = 0,
                             max = 100)
  
  # Sample first observed variable from the hidden variable
  observed_variable[1] <- sample_emission_model(z_t = state_variable[1],
                                                sd = sd_emission)
  
  # For every time step t = 1,...,T
  # Sample one hidden variable and one observed variable
  for (t in 2:T) {
    state_variable[t] <- sample_transition_model(z_t_1 = state_variable[t - 1])
    observed_variable[t] <- sample_emission_model(z_t = state_variable[t],
                                                  sd = sd_emission)
  }
  
  return(list(z = state_variable, x = observed_variable))
}

# Name: particle_filter
# Input:
#   M: Number of particles
#   obs: Observations for each time step
particle_filter <- function(M, obs, sd_emission = 1, sd_transition = 1, fix_weights = 0) {
  
  T <- length(obs)
  
  # Create matrix that will contain all
  # particles for each time step
  particles <- matrix(NA,
                      nrow = T,
                      ncol = M)
  temp_particles <- NULL
  
  # Generate initial M particles,
  # Uniform(0, 100)
  particles[1, ] <- runif(n = M,
                     min = 0,
                     max = 100)
  
  weights <- NULL
  
  for (t in 2:T) {
    for (m in 1:M) {
      temp_particles[m] <- sample_transition_model(z_t_1 = particles[t - 1, m])
      
      if (fix_weights == 0) {
        weights[m] <- density_emission_model(x_t = obs[t],
                                             z_t = temp_particles[m],
                                             sd = sd_emission)
      } else {
        weights[m] <- fix_weights
      }

    }

    particles[t, ] <- sample(x = temp_particles,
                           size = M,
                           replace = TRUE,
                           prob = weights)
  }
  
  return (particles)
}

# Name: visualize_particles
# Input:
#   particles: Particles for each time step
#   obs: Observations for each time step
#   t: Current time step
visualize_particles <- function(particles, true, t) {
  M <- dim(particles)[2]
  
  hist(x = particles[t, ], 
       breaks = 20,
       main = paste('t = ', t), 
       xlab = 'Particles')
  points(x = particles[t, ],
         y = rep(0, M),
         pch = 16,
         col = rgb(0, 0, 0, 0.3))
  abline(v = mean(particles[t, ]), col = 'red')
  abline(v = true[t], col = 'blue')
  legend('topright',
         legend = c("Particles", "Particle mean", "Observation"),
         col = c(rgb(0, 0, 0, 1), 'red', 'blue'),
         pch = c(16, NA, NA),
         lty = c(NA, 1, 1))
}

################ Setup ################

################ Task 1 ################
T <- 100
M <- 100

# Generate sample data
samples <- sample_data(T)

# Generate particles for each time step
particles <- particle_filter(M = M,
                             obs = samples$x)

# Visualization of particles, sd = 1, t = {1, 40, 80, 100}
par(mfrow = c(2, 2))
visualize_particles(particles = particles, true = samples$z, t = 1)
visualize_particles(particles = particles, true = samples$z, t = 40)
visualize_particles(particles = particles, true = samples$z, t = 80)
visualize_particles(particles = particles, true = samples$z, t = 100)

################ Task 2 ################

# Generate sample data
samples_5 <- sample_data(T,
                         sd_emission = 5)
samples_50 <- sample_data(T,
                           sd_emission = 50)

# Generate particles for each time step
particles_5 <- particle_filter(M = M,
                             obs = samples_5$x,
                             sd_emission = 5)
particles_50 <- particle_filter(M = M,
                               obs = samples_50$x,
                               sd_emission = 50)

# Visualization of particles, sd = 5, t = {1, 40, 80, 100}
par(mfrow = c(2, 2))
visualize_particles(particles = particles_5, true = samples$z, t = 1)
visualize_particles(particles = particles_5, true = samples$z, t = 40)
visualize_particles(particles = particles_5, true = samples$z, t = 80)
visualize_particles(particles = particles_5, true = samples$z, t = 100)

# Visualization of particles, sd = 50, t = {1, 40, 80, 100}
par(mfrow = c(2, 2))
visualize_particles(particles = particles_50, true = samples$z, t = 1)
visualize_particles(particles = particles_50, true = samples$z, t = 40)
visualize_particles(particles = particles_50, true = samples$z, t = 80)
visualize_particles(particles = particles_50, true = samples$z, t = 100)

################ Task 3 ################

# Generate sample data
samples_5 <- sample_data(T,
                         sd_emission = 1)

# Generate particles for each time step
particles_5 <- particle_filter(M = M,
                               obs = samples_5$x,
                               sd_emission = 1,
                               fix_weights = 1)

# Visualization of particles, sd = 1,fixed weights = 1, t = {1, 40, 80, 100}
par(mfrow = c(2, 2))
visualize_particles(particles = particles_50, true = samples$z, t = 1)
visualize_particles(particles = particles_50, true = samples$z, t = 40)
visualize_particles(particles = particles_50, true = samples$z, t = 80)
visualize_particles(particles = particles_50, true = samples$z, t = 100)