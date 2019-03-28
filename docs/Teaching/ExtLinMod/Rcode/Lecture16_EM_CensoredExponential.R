n = 1000
Y <- rexp(n, 1/10)

cc <- 10
(true.mle <- sum(Y)/n)

###########################
# Censor the observations
Y[Y>cc] <- NA 
mean(Y, na.rm=T) ## mean of uncensored observations

r = sum(!is.na(Y))

# MLE using the full likelihood that takes the missing data mechanism into account
(sum(Y, na.rm=T) + (n-r)*cc)/r


##################################
# EM
##################################

theta.em <- c()
theta.em[1] <- 5 # initial guess

i <- 1
change <-100
while(change >10^-10){
  Y.em <- Y
  Y.em[is.na(Y.em)] <- theta.em[i] + cc
  theta.em[i+1] <- sum(Y.em)/n
  change <- abs(theta.em[i+1]-theta.em[i])
  i <- i+1
  }

# note the fast convergence
theta.em

# these should agree  
theta.em[length(theta.em)] ## estimate from the EM
(sum(Y, na.rm=T) + (n-r)*cc)/r



