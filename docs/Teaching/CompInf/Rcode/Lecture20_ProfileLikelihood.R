


###################
## Normal example
#
# Suppose X_i ~ N(mu, sigma^2)
#
#
#############################################

mu = 0
sigma2=2
n=100

X <- rnorm(n, mu, sqrt(sigma2))

# Profile likelihood for mu

profile_loglike_mu <- function(mu){

  sigma2hat <- mean((X-mu)^2)
  -n/2*  log(sigma2hat)-n/2
}

mu.vals <- seq(-3,3,0.01)
prof.vals <- sapply(mu.vals, profile_loglike_mu)
plot(mu.vals, prof.vals, type='l')
 # try changing sigma2


l2 <- function(sigma){ # profile likelihood for sigma
  -n*log(sigma) -1/(2*sigma^2)*sum((x-mean(x))^2)
}

sigma <- seq(0,20, 0.1)

plot(sigma, sapply(sigma,l2), type='l')
(opt.out2 =optimise(l2, c(0,20), maximum = TRUE))
abline(v=opt.out2$max, col=2)






###################################################################
# Weibull example
#################################################################


## define the data
alldata<-c(6,6,6,6,7,9,10,10,11,13,16,17,19,20,22,23,25,32,32,34,35)
uncensored<-c(6,6,6,7,10,13,16,22,23)
d<-9


## function to compute the negative loglikelihood
## use with nlmin to find the m.l.e.
## Parameters are log(alpha) and log(beta)
 
weibullloglhood<-function(logtheta,uncensored,alldata){
  alpha<-exp(logtheta[1])
  beta<-exp(logtheta[2])
  -(9*log(alpha)+alpha*9*log(beta)+(alpha-1)*sum(log(uncensored))-beta^alpha*sum(alldata^alpha))  
}


## obtain the m.l.e
theta<-exp(optim(c(0,0),weibullloglhood,uncensored=uncensored,alldata=alldata)$par)
theta

## function to compute the profile loglikelihood
profileloglhood<-function(alpha){
    nalpha<-length(alpha)   
    M<-matrix(alldata,nrow=21,ncol=nalpha)^matrix(alpha,nrow=21,ncol=nalpha,byrow=T)
    9*log(alpha)+alpha*9*log((d/apply(M,2,sum))^(1/alpha))+(alpha-1)*sum(log(uncensored))-d
}

## plot the profile loglikehood
alpha<-seq(from=0.1,to=3,length=100)
plot(alpha,profileloglhood(alpha),type="l")
(opt.out <- optimize(profileloglhood, lower=0, upper=100, maximum=TRUE))
alpha.hat <-opt.out$maximum
plot(alpha,2*(profileloglhood(alpha.hat)-profileloglhood(alpha)),type="l", main='deviance(alpha)', ylab='profile deviance')
abline(h=qchisq(0.95,df=1), lty=2, col=2 )

## function to calculate profile deviance and 
## subtract 3.841. Minimise this to obtain confidence limits
profiledeviance<-function(alpha){
  (2*(profileloglhood(1.35)-profileloglhood(alpha))-3.841)^2
}

#find where this is zero
plot(alpha,profiledeviance(alpha), type='l', ylim=c(0,20))



## Find the lower and upper confidence limits
(lower = nlminb(0.5,profiledeviance,lower=0.1,upper=1.35)$par)
(upper = nlminb(2,profiledeviance,lower=1.35,upper=5)$par)
plot(alpha,2*(profileloglhood(alpha.hat)-profileloglhood(alpha)),type="l")
abline(h=qchisq(0.95,df=1), lty=2, col=2 )
abline(v=lower, lty=2, col=3)
abline(v=upper, lty=2, col=3)

## Compare with 95% interval based on asymptotic normality
## of the m.l.e. See chapter 4.

Varmatrix<-function(alpha,beta){
    v11<- -9/alpha^2-beta^alpha*(log(beta))^2*sum(alldata^alpha)-2*beta^alpha*log(beta)*sum(alldata^alpha*log(alldata)-beta^alpha*sum(alldata^alpha*(log(alldata))^2))
    v22<- 1/beta^2*(beta^alpha*alpha*(1-alpha)*sum(alldata^alpha)-alpha*9)
    v12<- 1/beta*(9-beta^alpha*log(beta)*alpha*sum(alldata^alpha)-beta^alpha*sum(alldata^alpha)-beta^alpha*alpha*sum(alldata^alpha*log(alldata)) )
V<- -solve(matrix(c(v11,v12,v12,v22),nrow=2,ncol=2))

}

V<-Varmatrix(theta[1],theta[2])
c(theta[1]-1.96*V[1,1]^0.5,theta[1]+1.96*V[1,1]^0.5)

abline(v=0.9964,lty=2,col=2)
abline(v=1.71,lty=2,col=2)


#############################################
#
# Machine component failure
#
#





profile_loglike_beta < 

