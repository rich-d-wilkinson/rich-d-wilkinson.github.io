#
# Let's start by creating some data from the model
# phi *N(mu, sigma^2) + (1-phi)*N(nu, tau^2)
#

# specify true values
phi.t = 0.3
mu.t=-1
nu.t=3
sigma.t = 0.5 # Assume known
tau.t = 3

n=10^3 # number of datapoints

Y <- rbinom(n, 1,phi.t)

X <- rnorm(n, mu.t, sigma.t)*Y + rnorm(n, nu.t, tau.t)*(1-Y)

hist(X, prob=T)
density.func <- function(x) phi.t*dnorm(x, mu.t, sigma.t) + (1-phi.t)*dnorm(x, nu.t, tau.t)
lines(x=seq(-10,20,0.1), y=density.func(seq(-10,20,0.1)))

############################################################
#
#
# Let's now try the EM algorithm to see whether we can recover these parameters.
#
#
############################################################

calc_pi <- function(phi, mu, nu, x){
  numerator1 <- phi
  numerator2 <- dnorm(x, mu, sigma.t)
  denom <- phi*dnorm(x,mu, sigma.t) + (1-phi)* dnorm(x,nu, tau.t)
  return(numerator1*numerator2/denom)
}

# write a loop to calculate all the p values
calc_p <- function(phi,mu,nu){
  tmp <- sapply(X, function(x) calc_pi(phi,mu,nu,x))
  return(tmp)
}
 

calc_theta <- function(p){
  phi <- sum(p)/n
  mu <- sum(p*X)/sum(p)
  nu = sum((1-p)*X)/sum(1-p)
  return(c(phi, mu,nu))
}

EM <- function(n_step, phi0, mu0, nu0){
  phi <- c()
  phi[1] <- phi0
  mu <- c()
  mu[1] <- mu0
  nu <- c()
  nu[1] <- nu0
  
  for(nn in 1:n_step){
    pp <- calc_p(phi[nn], mu[nn], nu[nn])
    theta.new <- calc_theta(pp)
    phi[nn+1] <- theta.new[1]
    mu[nn+1] <- theta.new[2]
    nu[nn+1] <- theta.new[3]
  }
  return(cbind(phi, mu, nu))
}


EM(100, 0.5, 0,0)


#######################
#
# Real example using the iris dataset
#
###########################

plot(iris[, 1:4], col = iris$Species, pch = 18, main = "Fisher's Iris Dataset")
library(ggplot2)
qplot(iris$Sepal.Length, colour = iris$Species)

library(dplyr)
virg = filter(iris, Species=='virginica')
setosa = filter(iris, Species=='setosa')

#Let's suppose we've lost the species label
X <- c(virg$Sepal.Length, setosa$Sepal.Length)

hist(X, xlab='sepal length')

# Let's now model this as a mixture distribution. I'm going to cheat a little here, 
# and assume we know the variances of the two groups, as I've only presented
# the version of the algorithm where the variances are known.
sigma.t =sd(virg$Sepal.Length) 
tau.t = sd(setosa$Sepal.Length)

# As a starting pointLet's assume the two species are equally prevalent in the data
# and that the two means are 4 and 6
EM(100, 0.5, 4,8)

# These estimates are very close to the means from the data sets split 
#into the two species
mean(virg$Sepal.Length)
mean(setosa$Sepal.Length)

# With a little more work, we could fit three Gaussians for all three species of iris
# and estimate the variances as well.

#####################################################

# Thankfully there is are R packges for doing this:
library(mclust)
mcl.fit <- Mclust(iris$Sepal.Length, G=3)
summary(mcl.fit, parameters=TRUE)
plot(mcl.fit, what='classification')
plot(mcl.fit, what='density')
plot(mcl.fit, what='uncertainty')

# we can also allow different variances:
mcl.fit2 <- Mclust(iris$Sepal.Length, G=3, modelNames='V')
summary(mcl.fit2, parameters=TRUE)
plot(mcl.fit2, what='classification')
plot(mcl.fit2, what='density')


# We can do the same thing but using more than just the sepal length.
# Select 4 continuous variables and look for three distinct groups.
mcl.model <- Mclust(iris[, 1:4], 3)
summary(mcl.model, parameters=TRUE)
# Plot our results.
plot(mcl.model, what = "classification", main = "Mclust Classification")
# compare with
plot(iris[, 1:4], col = iris$Species, pch = 18, main = "Fisher's Iris Dataset")
