
#
# Lets estimate P(X < -1 ) where X ~ Cauchy
#

n <- 100
x <- rcauchy(n)

sum(x< -1 )/n
pcauchy(-1)

# What is the error in this estimator?

MCestimates <- replicate(10^5, {
  x <- rcauchy(n)
  sum(x< -1 )/n
})


hist(MCestimates, probability = T)
abline(v=pcauchy(-1), col=2)

# MC error
sd(MCestimates)

# Check the normal assumption holds
lines(seq(0,1,0.01), dnorm(seq(0,1,0.01), mean = mean(MCestimates), sd = sd(MCestimates)), col=3)


## plot error rate as function of n

i<-1
error <- c()
n.values <- c(100,200,400,800,1600,3200, 6400)
for(n in n.values){
  print(n)
  
  error[i] <- sd(replicate(10^4, {
    x <- rcauchy(n)
    sum(x< -1 )/n
  }))
  i <-i+1
  }
plot(n.values, error)


## check that ther log(error) decays with rate -1/2
plot(n.values, error, log='xy', xlab='log(n)', ylab='log(RMSE)')
lm(log(error)~log(n.values))
abline(a=log10(sd(as.numeric(rcauchy(100000)< -1))), b=-1/2)


##################################################

