########################################
### Lab 1 - Graphical Models
### By Max Fischer
### TDDE15 - Advanced Machine Learning
### Linköpings university
########################################
############### Functions ##############
predict_from_network <- function(junc_tree, data, obs_variables, pred_variable) {
  for (i in 1:dim(data)[1]) {
    X <- NULL
    for (j in obs_variables) {
      X[j] <- if(data[i, j] == "yes") "yes" else "no"
    }
    
    # Set evidence in junction tree for observation i
    # We have observations of all variables except:
    # S: If a person smokes or not
    evidence <- setEvidence(object = junc_tree,
                            nodes = obs_variables,
                            states = X)
    
    # Do prediction of S from junction tree with the above evidence
    prob_dist_fit <- querygrain(object = evidence,
                                nodes = pred_variable)$S
    prediction_fit[i] <- if (prob_dist_fit["yes"] >= 0.5) "yes" else "no"
  }
  
  return (prediction_fit)
}

################ Setup ################
# Install package if not installed
if (!require(bnlearn)) {
  install.packages("bnlearn")
}

# Import package
library(bnlearn)

# Import data
# Variables:
# A, S, T, L, B, E, X, D
data("asia")

################ Task 1 ################
# It is possible to show that the hill-climbing algorithm
# can returns non-equivilant Bayesian Network structures.
# In our case, we show it by using two different
# score-functions in the HC-algorithm

par(mfrow = c(1, 2))

# Bayesian Network 1
set.seed(12345)
hc_network1 <- hc(x = asia,
                  start = BN_init_1, restart = 5)
hc_network1 <- cpdag(hc_network1) # Finds the equivilance class of the BN
plot(hc_network1)

# Bayesian Network 2
set.seed(12345)
hc_network2 <- hc(x = asia,
                  score = 'aic')
hc_network2 <- cpdag(hc_network2)
plot(hc_network2)

# Check equality
print(all.equal(hc_network1, hc_network2))

################ Task 2 ################
# Install packages if not installed
if (!require(gRain)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("RBGL")
}

library(gRain)

# Split Asia dataset into 80 % training and 20 % testing
# Learn a BN on the training set.

# Split dataset into train and test
train_indices <- sample(x = seq(1, dim(asia)[1], 1), 
                        size = dim(asia)[1]*0.8,
                        replace = FALSE)
asia.train <- asia[train_indices,]
asia.test <- asia[-train_indices,]

# Learn BN-network structure
BN.structure <- hc(x = asia.train,
                   restart = 5)
BN.structure_true <- model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")

# Fit parameters to train data
BN.fit <- bn.fit(x = BN.structure,
                 data = asia.train)
BN.fit_true <- bn.fit(x = BN.structure_true,
                      data = asia.train)

# Convert fit to gRain-object
BN.fit_gRain <- as.grain(BN.fit)
BN.fit_true_gRain <- as.grain(BN.fit_true)

# Compile BN
# Creating a junction tree (Lauritzen-Spiegelhalter algorithm) and establishing clique potentials
junc_tree <- compile(BN.fit_gRain)
junc_tree_true <- compile(BN.fit_true_gRain)

# Predict S from Bayesian Network and test data observations
prediction_fit <- predict_from_network(junc_tree = junc_tree,
                                       data = asia.test,
                                       obs_variables = c("A", "T", "L", "B", "E", "X", "D"),
                                       pred_variable = c("S"))
prediction_fit_true <- predict_from_network(junc_tree = junc_tree_true,
                                            data = asia.test,
                                            obs_variables = c("A", "T", "L", "B", "E", "X", "D"),
                                            pred_variable = c("S"))

# Calculate confusion matricies
confusion_matrix_fit <- table(prediction_fit, asia.test$S)
print(confusion_matrix_fit)
confusion_matrix_true <- table(prediction_fit_true, asia.test$S)
print(confusion_matrix_true)

################ Task 3 ################
# Predict S only by observations of the Markov blanket.
# It should result in the same answer, as observations of 
# S's Markov blanket makes S independent of all the other
# variables

# Extract Markov blanket from fitted BN structure generated from HC
# and from fitted BN structure generated from true model
MB_fit <- mb(x = BN.fit,
             node = c("S"))
MB_fit_true <- mb(x = BN.fit_true,
                  node = c("S"))

prediction_fit_MB <- predict_from_network(junc_tree = junc_tree,
                                          data = asia.test,
                                          obs_variables = MB_fit,
                                          pred_variable = c("S"))
prediction_fit_true_MB <- predict_from_network(junc_tree = junc_tree_true,
                                               data = asia.test,
                                               obs_variables = MB_fit_true,
                                               pred_variable = c("S"))

# Calculate confusion matricies
confusion_matrix_fit <- table(prediction_fit, asia.test$S)
confusion_matrix_fit_true <- table(prediction_fit_true, asia.test$S)

################ Task 4  ################
# Learn the structure and parameters on the train data using a Naive Bayes Bayesian Network
# In Naive Bayes, we assume that the predictive variables are independent given the
# class variable (in this case S).
# p(A|S)*p(T|S)*p(L|S)*p(B|S)*p(E|S)*p(X|S)*p(D|S)*p(S)
# As a Bayesian Network it is represented by S being the parent to all other nodes.
# If S is observed, all the predictive variables are independent.

# Create Native Bayes Bayesian Network
naive_bayes_structure <- model2network("[S][A|S][T|S][L|S][B|S][E|S][X|S][D|S]")

# Fit parameters of network to train data
BN.fit_naive_bayes <- bn.fit(x = naive_bayes_structure,
                             data = asia.test)

# Convert fit to gRain-object
BN.fit_naive_bayes_grain <- as.grain(BN.fit_naive_bayes)

# Generate juncion tree and clique potentials
junc_tree_naive_bayes <- compile(BN.fit_naive_bayes_grain)

prediction_fit_naive_bayes <- predict_from_network(junc_tree = junc_tree_naive_bayes,
                                                   data = asia.test,
                                                   obs_variables = c("A", "T", "L", "B", "E", "X", "D"),
                                                   pred_variable = c("S"))

prediction_fit_true <- predict_from_network(junc_tree = junc_tree_true,
                                                   data = asia.test,
                                                   obs_variables = c("A", "T", "L", "B", "E", "X", "D"),
                                                   pred_variable = c("S"))

# Calculate confusion matricies
confusion_matrix_naive_bayes <- table(prediction_fit_naive_bayes, asia.test$S)
confusion_matrix_fit_true <- table(prediction_fit_true, asia.test$S)