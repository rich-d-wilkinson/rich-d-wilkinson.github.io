library(magic)
library(mvtnorm)
library(lme4)

# Rail data example. Model is

# y_ij = mu + b_i + epsilon_{ij}
load('MAS473.RData')
# b_i ~ N(0,psisq)
# epsilon_{ij} ~ N(0,sigmasq)
attach(raildata)
library(ggplot2)
qplot(Rail, travel, geom='boxplot')

fm1.ml<-lmer(travel~1+(1|Rail),raildata,REML=F)
fm1.reml<-lmer(travel~1+(1|Rail),raildata)

logLik(fm1.ml)
logLik(fm1.reml)


###############################################
# Fit model using ordinary maximum likelihood #
###############################################

(fm1.ml<-lmer(travel~1+(1|Rail),raildata,REML=F))

# Now try to get the same parameter estimates by numerically maximising the log likelihood

# Define a function to calculate (minus) the log likelihood

minus.log.ordinary.likelihood.raildata<-function(theta,y){
  
  # theta[1] = log sigmasq 
  # theta[2] = log psisq
  
  sigmasq<-exp(theta[1])
  psisq<-exp(theta[2]) # (force sigmasq and psisq to be positive)
  
  V1<-matrix(psisq,3,3)
  diag(V1)<-diag(V1)+sigmasq
  V<-adiag(V1,V1,V1,V1,V1,V1) # Variance covariance matrix of data
  X<- model.matrix(travel~1,raildata)
  
  -dmvnorm(t(y),X*theta[3],V,log=T) # -log likelihood
}

y<-matrix(travel,18,1)
theta<-c(log(10),log(500),50)



# minimise - log likelihood. Compare estimates with fm1.ml

theta.mle<-optim(c(log(10),log(500),50), minus.log.ordinary.likelihood.raildata,y=matrix(travel,18,1))

exp(theta.mle$par[1:2]) # estimates of sigmasq and psisq
theta.mle$par[3] # estimate of mu
summary(fm1.ml)
# likely to be small discrepancies due to optimisation routine. 
# Can try different starting values to get global maximum


########################
# Fit model using REML #
########################

(fm1.reml<-lmer(travel~1+(1|Rail),raildata))

minus.log.reml.likelihood.raildata<-function(theta,y){
  sigmasq<-exp(theta[1])
  psisq<-exp(theta[2])
  V1<-matrix(psisq,3,3)
  diag(V1)<-diag(V1)+sigmasq
  V<-adiag(V1,V1,V1,V1,V1,V1)
  X<- model.matrix(travel~1,raildata)
  Vinv<-solve(V)
  
  betahat<-solve(t(X)%*% Vinv %*% X)%*% t(X) %*% Vinv %*% y
  
 -( -0.5*log(det(V)) - 0.5*log(det(t(X)%*%Vinv%*%X)) - 0.5*(18-1)*log(2*pi) - 0.5*t(y-X%*%betahat)%*%Vinv%*%(y-X%*%betahat))
    

}
# Example to show construction of V, betahat and log REML criterion
sigmasq<-10
psisq<-100
V1<-matrix(psisq,3,3)
diag(V1)<-diag(V1)+sigmasq
(V<-adiag(V1,V1,V1,V1,V1,V1))

X<- model.matrix(travel~1,raildata)
y<-matrix(travel,18,1)
Vinv<-solve(V)
(betahat<-solve(t(X)%*% Vinv %*% X)%*% t(X) %*% Vinv %*% y)
 -0.5*log(det(V)) - 0.5*log(det(t(X)%*%Vinv%*%X)) - 0.5*(18-1)*log(2*pi) - 0.5*t(y-X%*%betahat)%*%Vinv%*%(y-X%*%betahat)

# Optimise log REML likelihood
theta.mle.reml<-optim(log(c(10,500)), minus.log.reml.likelihood.raildata,y=matrix(travel,18,1))
exp(theta.mle.reml$par[1:2]) # estimates of sigmasq and psisq

2*theta.mle.reml$value




