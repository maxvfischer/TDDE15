########################################
### Lab 2 - Hidden Markov Models
### By Max Fischer
### TDDE15 - Advanced Machine Learning
### Linköpings university
########################################
############### Functions ##############

# Name: get_alpha
# Input:
#   possible_states: States that a hidden variable can be
#   emiss_model: Emission model (matrix)
#   trans_model: Transition model (matrix)
#   obs: Observation in time step t
#   alpha_prev: Alpha in time step t-1
get_alpha <- function(possible_states, emiss_model, trans_model, obs, alpha_prev) {
  
  alpha <- NULL
  
  # For each possible state (sector)
  # This is every state in X_t
  for (state in possible_states) {
      
    # Given a state, how likely is our
    # observation
    emission <- emiss_model[state, obs]
    
    # For every possible state (t_1), how
    # likely is each of the new states
    #
    # This is then weighted with alpha from
    # t_1 of each state. Alpha is some kind
    weighted_transition <- sum(
      sapply(possible_states, function(z) {
        alpha_prev[z] * trans_model[z, state]
      })
    )
    
    # Alpha for every possible state in Z_t
    alpha[state] <- emission * weighted_transition
  }
  return (alpha)
}

# Name: get_beta
# Input:
#   possible_states: States that a hidden variable can be
#   emiss_model: Emission model (matrix)
#   trans_model: Transition model (matrix)
#   obs: Observation in time step t+1
#   alpha_prev: Alpha in time step t+1
get_beta <- function(possible_states, trans_model, emiss_model, obs_next, beta_next) {
  
  beta <- NULL
  
  # For each possible state (sector)
  # This is every state in X_t
  for (state in possible_states) {
    beta[state] <- sum(
      sapply(possible_states, function(z) {
        beta_next[z] * emiss_model[z, obs_next] * trans_model[state, z]
      })
    )
  }
  return (beta)
}

# Name: forward_backward
# Input:
#   possible_states: States that a hidden variable can be
#   trans_model: Transition model (matrix)
#   emiss_model: Emission model (matrix)
#   obs_vars: Observed variables
forward_backward <- function(possible_states, trans_model, emiss_model, obs_vars) {
  T <- length(obs_vars)
  alpha <- matrix(NA, nrow = T, ncol = length(possible_states))
  beta <- matrix(NA, nrow = T, ncol = length(possible_states))
  
  # Forward
  # Initialize alpha: 
  # p(Z_0) (place robot in any state with 1/10 probability)
  first_obs <- obs_vars[1]
  initial <- rep(0.1, 10)
  alpha[1, ] <- emiss_model[, first_obs] * initial # alpha[1] = p(x_0|Z_0)*p(Z_0)
  
  # For every time step t = 2,...,T
  for (t in 2:T) {
    # Generate alpha for time step t
    alpha[t, ] <- get_alpha(possible_states = possible_states,
                        emiss_model = emiss_model,
                        trans_model = trans_model,
                        obs = obs_vars[t],
                        alpha_prev = alpha[t - 1, ])
  }
  
  # Backward
  # Initalize beta
  beta[T, ] <- 1
  
  # For every time step t = (T-1),...,1
  for (t in (T-1):1) {
    # Generate beta for time step t
    beta[t, ] <- get_beta(possible_states = possible_states,
                          emiss_model = emiss_model,
                          trans_model = trans_model,
                          obs = obs_vars[t+1],
                          beta_next = beta[t+1, ])
  }
  return (list(alpha = alpha, beta = beta))
}

# Name: filtering
# Input:
#   alphas: Alpha probabilities generated in the Forward Backward algorithm
filtering <- function(alphas) {
  filtering <- matrix(NA,
                      nrow = dim(alphas)[1],
                      ncol = dim(alphas)[2])
  
  for (t in 1:dim(alphas)[1]){
    filtering[t, ] <- alphas[t, ]/sum(alphas[t, ])
  }

  return(filtering)
}

# Name: smoothing
# Input:
#   alphas: Alpha probabilities generated in the Forward Backward algorithm
#   betas: Beta probabilities generated in the Forward Backward algorithm
smoothing <- function(alphas, betas) {
  smoothing <- matrix(NA,
                      nrow = dim(alphas)[1],
                      ncol = dim(alphas)[2])
  
  for (t in 1:dim(alphas)[1]) {
    smoothing[t, ] <- (alphas[t, ] * betas[t, ]) / (sum(alphas[t, ] * betas[t, ]))
  }
  
  return (smoothing)
}

# Name: viterbi_algo
# Input:
#   possible_states: States that a hidden variable can be
#   trans_model: Transition model (matrix)
#   emiss_model: Emission model (matrix)
#   obs_vars: Observed variables
viterbi_algo <- function(possible_states, trans_model, emiss_model, obs_vars) {
  T <- length(obs_vars)
  
  omega <- matrix(NA,
                  nrow = T,
                  ncol = length(possible_states))
  psi <- matrix(NA,
                  nrow = T,
                  ncol = length(possible_states))
  Z <- NULL
  
  # Initialize omega_0
  initial <- rep(0.1, 10)
  first_obs <- obs_vars[1]
  omega_init <- log(initial) + log(emiss_model[, first_obs])

  for (t in 0:(T-1)) {
    obs <- obs_vars[t + 1]
    
    for (state in possible_states) {
      trans_omega <- sapply(possible_states, function(z) {
        if (t == 0) {
          log(trans_model[z, state]) + omega_init
        } else {
          log(trans_model[z, state]) + omega[t, z]
        }

      })
      max <- max(trans_omega)
      
      omega[t + 1, state] <- emiss_model[state, obs] + max
      psi[t + 1, state] <- which.max(trans_omega)
    }
  }
  
  Z[T] <- which.max(omega[T, ])
  
  for (t in (T-1):0) {
    Z[t] <- psi[t+1, Z[t+1]]
  }
  
  return (Z)
}

# Name: viterbi_algo
# Input:
#   prediction: Vector of predictions
#   true: Vector of true states
accuracy <- function(prediction, true) {
  confusion_matrix <- table(prediction, true)
  accuracy <- sum(diag(confusion_matrix))/sum(confusion_matrix)
  return (accuracy)
}

################ Setup ################
# Install and load necessary packages
if (!require(HMM)) {
  install.packages("HMM")
}

library(HMM)

if (!require(entropy)) {
  install.packages("entropy")
}

library(entropy)
################ Task 1 ################
# A robot walks around a ring. The ring
# is divided into 10 sectors. The robot
# is in one sector at any given time step
# and it's equal probable for the robot
# to stay in the state as it is to move
# to the next state.

# The robot has a tracking device.
# If the robot is in sector i, the tracking
# device will report that the robot is
# in the sectors [i - 2, i + 2] with
# equal probability.

# Create transition matrix from description above
# Each row consists of: p(z^t|z^t-1), t = 1, ..., 10
transition_vector <- c(0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0.5, 0.5, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0.5, 0.5, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0.5, 0.5, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0.5, 0.5, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0.5, 0.5, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0.5, 0.5, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.5,
                       0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0.5)
transition_matrix <- matrix(data = transition_vector,
                              nrow = 10,
                              ncol = 10)

# Create emission matrix from description above
emission_vector <- c(0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0.2, 0.2,
                     0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0.2,
                     0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0,
                     0, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0,
                     0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0,
                     0, 0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0,
                     0, 0, 0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0,
                     0, 0, 0, 0, 0, 0.2, 0.2, 0.2, 0.2, 0.2,
                     0.2, 0, 0, 0, 0, 0, 0.2, 0.2, 0.2, 0.2,
                     0.2, 0.2, 0, 0, 0, 0, 0, 0.2, 0.2, 0.2)
emission_matrix <- matrix(data = emission_vector,
                          nrow = 10,
                          ncol = 10)

# Initialize a Hidden Markov Model (HMM)
# States are the hidden variables
HMM.states <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

# Symbols are the observable variables
HMM.symbols <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

HMM <- initHMM(States = HMM.states,
                    Symbols = HMM.symbols,
                    startProb = rep(0.1, 10),
                    transProbs = transition_matrix,
                    emissionProbs = emission_matrix)

################ Task 2 ################
# Simulate 100 time steps from the model
# above
T <- 100
HMM.sim <- simHMM(hmm = HMM,
                  length = T)

################ Task 3 ################

# Generate alpha and beta probabilities
# from forward-backward algorithm
FB <- forward_backward(possible_states = 1:10,
                       trans_model = transition_matrix,
                       emiss_model = emission_matrix,
                       obs_vars = HMM.sim$observation)

# Filtering probabilities
filtering_100 <- filtering(FB$alpha)

# Smoothing probabilities
smoothing_100 <- smoothing(alphas = FB$alpha, betas = FB$beta)

# Viterbi - Most probable path
viterbi_100 <- viterbi_algo(possible_states = 1:10,
                        trans_model = transition_matrix,
                        emiss_model = emission_matrix,
                        obs_vars = HMM.sim$observation)

################ Task 4 ################

# Predict based on the highest probability
filtering_prediction_100 <- apply(filtering_100, MARGIN = 1, which.max)
smoothing_prediction_100 <- apply(smoothing_100, MARGIN = 1, which.max)

# Calculate accuracy from prediction
accuracy_filtering_100 <- accuracy(prediction = filtering_prediction,
                               true = HMM.sim$states)
accuracy_smoothing_100 <- accuracy(prediction = smoothing_prediction,
                               true = HMM.sim$states)

################ Task 5 ################
# .....

################ Task 6 ################

# Simulate 200 observations and hidden state
# from the Hidden Markov Model
HMM.sim_200 <- simHMM(hmm = HMM,
                      length = 200)

# Generating alpha and beta predictions
FB_200 <- forward_backward(possible_states = 1:10,
                           trans_model = transition_matrix,
                           emiss_model = emission_matrix,
                           obs_vars = HMM.sim_200$observation)

# Filtering probabilities
filtering_200 <- filtering(FB_200$alpha)

# Predict based on the highest probability
filtering_prediction_200 <- apply(filtering_200, MARGIN = 1, which.max)

# Entroy of 100 and 200 observations
entropy_100 <- entropy.empirical(filtering_prediction_100)
entropy_200 <- entropy.empirical(filtering_prediction_200)

################ Task 7 ################

# Generate probabilities of the hidden state
# for time step 101
t_101 <- transition_matrix %*% filtering_100[100, ]
