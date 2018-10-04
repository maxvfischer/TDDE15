##################################
###########   Lab 1   ############
##################################

## Packages
print(require(gRain))
if (!require(bnlearn)) {
  install.packages(bnlearn)
}

if (!require(gRain)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("RBGL")
  install.packages(gRain)
}

library(bnlearn)
library(gRain)

## Part 1

# Import dataset
# D (dyspnoea), a two-level factor with levels yes and no.
# T (tuberculosis), a two-level factor with levels yes and no.
# L (lung cancer), a two-level factor with levels yes and no.
# B (bronchitis), a two-level factor with levels yes and no.
# A (visit to Asia), a two-level factor with levels yes and no.
# S (smoking), a two-level factor with levels yes and no.
# X (chest X-ray), a two-level factor with levels yes and no.
# E (tuberculosis versus lung cancer/bronchitis), a two-level factor with levels yes and no.
data('asia')

# Hill climbing
# Restart set to 5. When a local maxima is found, a number of random changes is done to the graph
# and then the algorithm continues. It often leads to different graphs.
iter = 1000
no_of_equal = 0
percent_equal = numeric(0)

# Multiple iterations are done and the percentage of equal graphs
# are calculated during each iteration
for (i in 1:iter) {
  HC1 = hc(asia, restart=5)
  HC2 = hc(asia, restart=5)
  if (all.equal(HC1, HC2) == TRUE) {
    no_of_equal = no_of_equal + 1
  }
  percent_equal[i] <- no_of_equal/i
}

# The convergence of how many percent of the iterations that are equal are ploted
plot(percent_equal, type='l', xlab='Iteration', ylab='Percentage equal graphs')
abline(h = mean(percent_equal), col="red", lty=2)
legend(750, 1, legend=c("Percentage equal", "Mean"),
       col=c("black", "red"), lty=1:2, cex=0.8)

## Part 2

# Create train set (80%) and test set (20%)
train_size = dim(asia)[1]*0.8
train_indices = sample(nrow(asia), 4000)
asia.train <- asia[train_indices, ]
asia.test <- asia[-train_indices, ]
S <- c("S")

# Create Bayesian Network from train data
train.hc <- hc(asia.train) 
train.fit <- bn.fit(train.hc, asia.train) # Fit to train data
train.as_grain <- as.grain(train.fit) # Create gRain object
train.compiled <- compile(train.as_grain)

# Create Bayesian Network from true network
true = model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]") # True Bayesian network
true.fit <- bn.fit(true, asia.train)
true.as_grain <- as.grain(true.fit)
true.compiled <- compile(true.as_grain)

prediction_fit <- NULL
prediction_true <- NULL

# Loop for each observation in test set
for (i in 1:dim(asia.test)[1]) {
  
  # Create correct formated vector for each observation
  z <- NULL
  for (j in c("A", "T", "L", "B", "E", "X", "D")) {
    if (asia.test[i, j] == "no") {
      z <- c(z, "no")
    }
    else {
      z <- c(z, "yes")
    }
  }
  
  # Set evidence for train data and do prediction
  hc3 <- setEvidence(train.compiled, nodes=c("A", "T", "L", "B", "E", "X", "D"), states=z)
  x <- querygrain(hc3, c("S"))
  prediction_fit <- if(x$S["no"] > x$S["yes"]) c(prediction_fit, "no") else c(prediction_fit, "yes")
  
  # Set evidence for true Bayesian network and do prediction
  hc4 <- setEvidence(true.compiled, nodes=c("A", "T", "L", "B", "E", "X", "D"), states=z)
  k <- querygrain(hc4, c("S"))
  prediction_true <- if(x$S["no"] > x$S["yes"]) c(prediction_true, "no") else c(prediction_true, "yes")
}

# Confusion matrix
confusion_matrix_fit <- table(prediction_fit, asia.test$S)
confusion_matrix_true <- table(prediction_true, asia.test$S)

# Task 3
train.hc <- hc(asia.train)
train.mb <- mb(train.hc, "S")
train.fit <- bn.fit(train.hc, asia.train)
train.as_grain <- as.grain(train.fit)
train.compiled <- compile(train.as_grain)

# Create Bayesian Network from true network
true = model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]") # True Bayesian network
true.mb <- mb(true, "S")
true.fit <- bn.fit(true, asia.train)
true.as_grain <- as.grain(true.fit)
true.compiled <- compile(true.as_grain)

prediction_fit <- NULL
prediction_true <- NULL

# Loop for each observation in test set
for (i in 1:dim(asia.test)[1]) {
  
  # For each observation in the test set, create correct formated vector of the 
  # fitted model's Markov blanket
  z <- NULL
  for (j in train.mb) {
    if (asia.test[i, j] == "no") {
      z <- c(z, "no")
    }
    else {
      z <- c(z, "yes")
    }
  }
  
  # For each observation in the test set, create correct formated vector of the 
  # true model's Markov blanket
  z <- NULL
  for (j in true.mb) {
    if (asia.test[i, j] == "no") {
      z <- c(z, "no")
    }
    else {
      z <- c(z, "yes")
    }
  }
  
  # Set evidence for train data and do prediction
  hc3 <- setEvidence(train.compiled, nodes=train.mb, states=z)
  x <- querygrain(hc3, c("S"))
  prediction_fit <- if(x$S["no"] > x$S["yes"]) c(prediction_fit, "no") else c(prediction_fit, "yes")
  
  # Set evidence for true Bayesian network and do prediction
  hc4 <- setEvidence(true.compiled, nodes=true.mb, states=z)
  k <- querygrain(hc4, c("S"))
  prediction_true <- if(x$S["no"] > x$S["yes"]) c(prediction_true, "no") else c(prediction_true, "yes")
}

# Confusion matrix
confusion_matrix_fit <- table(prediction_fit, asia.test$S)
confusion_matrix_true <- table(prediction_true, asia.test$S)

# Task 4
