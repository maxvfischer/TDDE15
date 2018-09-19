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
for (i in seq(0, iter, 1)) {
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
