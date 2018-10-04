##################################
###########   Lab 2   ############
##################################

# Check if package is installed
# If not, install
if (!require(HMM)) {
  install.packages("HMM")
}

# Import library
library(HMM)

# Model the behavior of a robot that walks around a ring.
#
# At any given time, the robot is in one of the sectors and
# decides with equal probability to stay in that sector or
# move to the next sector.
# 
# You don't have direct observation of the robot, but it's
# equipped with a tracking device that you can access.
#
# The device is not accurate thought, it the robit is in the
# sector i, the device will report that the robit is in the
# sectors [i-2, i+2] with equal probability.

## Part 1
# Build a hidden Markov model (HMM) for the scenario above

# Create ring transition matrix
transition_probs_vector <- c(0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0.5, 0.5, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0.5, 0.5, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0.5, 0.5, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0.5, 0.5, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0.5, 0.5, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0.5, 0.5, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.5,
                      0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0.5)
transition_probs <- matrix(transition_probs_vector, 10, byrow=TRUE)

# Create ring emission matrix
emission_probs_vector <- c(0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0.2, 0.2,
                    0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0.2,
                    0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0,
                    0, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0,
                    0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0,
                    0, 0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0,
                    0, 0, 0, 0, 0.2, 0.2, 0.2, 0.2, 0, 0,
                    0, 0, 0, 0, 0, 0.2, 0.2, 0.2, 0.2, 0.2,
                    0.2, 0, 0, 0, 0, 0, 0.2, 0.2, 0.2, 0.2,
                    0.2, 0.2, 0, 0, 0, 0, 0, 0.2, 0.2, 0.2)
emission_probs <- matrix(emission_probs_vector, 10, byrow=TRUE)

# Declare states and symbols (observable "states")
states <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
symbols <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")

# Initialize HMM
hmm <- initHMM(States = states,
               Symbols = symbols,
               transProbs = transition_probs,
               emissionProbs = emission_probs)

## Task 2
# Simulate the HMM for 100 time steps
sim_hmm <- simHMM(hmm = hmm, length = 100)

## Task 3

# Filtering function
filtering <- function(hmm, obs) {
  log_alpha <- forward(hmm = hmm, observation = obs)
  
  # As forward generate log, we need to de-log
  alpha <- exp(log_alpha)
  
  # Use prop.table to normalize the probabilities
  filtering_prob <- prop.table(x = alpha, margin = 2)
  
  return (filtering_prob)
}

# Smoothning function
smoothning <- function(hmm, obs) {
  log_alpha <- forward(hmm, observation = obs)
  log_beta <- backward(hmm, observation = obs)
  
  # As forawrd and backward generate log, we need to de-log
  alpha <- exp(log_alpha)
  beta <- exp(log_beta)
  
  # Create un-normalized smoothning
  alpha_beta <- alpha*beta
  
  # Normalize smoothning
  smoothning_prob <- prop.table(x = alpha_beta, margin = 2)
  
  return (smoothning_prob)
}

# Calculate probability matrix from Filtering and Smoothning.
# Column: Time step
# Row: State
filtering_prob <- filtering(hmm = hmm, obs = sim_hmm$observation)
smoothning_prob <- smoothning(hmm = hmm, obs = sim_hmm$observation)

# Most probable path by Viterbi
viterbi <- viterbi(hmm = hmm, obs = sim_hmm$observation)

## Task 4

# Indices of most probable states (filtering and smoothning)
filtering_most_prob_index <- apply(filtering_prob, MARGIN = 2, which.max)
smoothning_most_prob_index <- apply(smoothning_prob, MARGIN = 2, which.max)

# Filtering accuracy
filtering_equal_to <- ifelse(states[filtering_most_prob_index] == sim_hmm$states, 1, 0)
filtering_accuracy <- sum(filtering_equal_to)/length(sim_hmm$states)

# Smoothning accuracy
smoothning_equal_to <- ifelse(states[smoothning_most_prob_index] == sim_hmm$states, 1, 0)
smoothning_accuracy <- sum(smoothning_equal_to)/length(sim_hmm$states)

# Viterbi accuracy
viterbi_equal_to <- ifelse(viterbi == sim_hmm$states, 1, 0)
viterbi_accuracy <- sum(viterbi_equal_to)/length(sim_hmm$states)

## Task 5

# Simulate 200 samples
sim_hmm_new <- simHMM(hmm = hmm, length = 200)

# Generate probability matrix from Filtering and Smoothning.
# Also most probable path by Viterbi
filtering_prob_new <- filtering(hmm = hmm, obs = sim_hmm_new$observation)
smoothning_prob_new <- smoothning(hmm = hmm, obs = sim_hmm_new$observation)
viterbi_new <- viterbi(hmm = hmm, observation = sim_hmm_new$observation)

# Indices of the most probable state (Filtering and Smoothning)
filtering_new_most_prob_index <- apply(filtering_prob_new, MARGIN = 2, FUN = which.max)
smoothning_new_most_prob_index <- apply(smoothning_prob_new, MARGIN = 2, FUN = which.max)

# Filtering accuracy
filtering_new_equal_to <- ifelse(states[filtering_new_most_prob_index] == sim_hmm_new$states, 1, 0)
filtering_new_accuracy <- sum(filtering_new_equal_to)/length(sim_hmm_new$states)

# Smoothning accuracy
smoothning_new_equal_to <- ifelse(states[smoothning_new_most_prob_index] == sim_hmm_new$states, 1, 0)
smoothning_new_accuracy <- sum(smoothning_new_equal_to)/length(sim_hmm_new$states)

# Viterbi accuracy
viterbi_new_equal_to <- ifelse(viterbi_new == sim_hmm_new$states, 1, 0)
viterbi_new_accuracy <- sum(viterbi_new_equal_to)/length(sim_hmm_new$states)

## Task 6

# Check if package "entropy" is installed
# If not -> Install
if (!require("entropy")) {
  install.packages(entropy)
}

# Include entropy library
library(entropy)

# Entropy of each time step (Filtering)
filtering_entropy <- apply(filtering_prob, MARGIN = 2, FUN = function(x) entropy.empirical(x))
plot(filtering_entropy, type='l',
     ylab='Entropy (Filtered)',
     xlab='Time step (t)',
     col='blue')

## Task 7

# Extract probabilities in t = 100
t_100 <- filtering_prob[,100]

# Generate probabilities in t = 101 by matrix multiplying
# transition matrix with the probability vector of t = 100
t_101 <- hmm$transProbs %*% t_100