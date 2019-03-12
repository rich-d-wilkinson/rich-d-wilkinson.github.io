

#####################
### Section 4.2.1 ###
#####################

# Example on page 45.

# 1. Generate some data. 
set.seed(1)
n<-10 # the sample size
x<-seq(from=1,to=10,length=n)
y<-0.2*x+rnorm(n) # use this for full model to be true
#y<-rnorm(10) # use this for reduced model to be true
plot(x,y,pch=16)

# 2. Fit the two linear models
lm.full<-lm(y~x)
lm.reduced<-lm(y~1)
abline(lm.reduced)
abline(lm.full)

# 3. Usual F-test to test H_0: lm.reduced is the true model
anova(lm.reduced,lm.full)

# 4. Compute maximum likelihood estimates for parameters in reduced model
(alpha.mle<-mean(y))
(reduced.sigmasq.mle<-mean((y-alpha.mle)^2))

# compare estimates with those in the fitted model - will they be the same?
lm.reduced$coefficients[1]
summary(lm.reduced)$sigma^2

# 5. Evaluate the log-likelihood for the reduced model
sum(log(dnorm(y,alpha.mle,reduced.sigmasq.mle^0.5)))

# direct way to get the log likelihood in R. 
logLik(lm.reduced)

# Compute log likelihood for full model
logLik(lm.full)

# 6. Perform the GLRT with chi-sq approximation
(obs.test.statistic<- -2*(logLik(lm.reduced) - logLik(lm.full)))
attributes(obs.test.statistic)<-NULL # get rid of the 'clutter'!

# Calculate the p-value. How does it compare with the p-value in step 3?
1-pchisq(obs.test.statistic,1)

curve(dchisq(x,df=1),from=0,to=10, lwd=2)
abline(v=obs.test.statistic, col="red", lwd=2)


# Plots to visualise why full model gives the higher likelihood
# Full model fits data more closely, so smaller m.l.e for sigma^2
# Estimated error distribution less 'spread out' for full model
# so likelihood terms (the density values at the residuals) will be greater

par(mfrow=c(2,2))
par(ps=15)

# Reduced model and residuals
plot(x,y, pch=16)
abline(lm.reduced)
for(i in 1:10){
	lines(c(x[i],x[i]), c(y[i],fitted(lm.reduced)[i]), lty=2)
}
e<-seq(from=-4,to=4,length=100)

# Estimated error distribution for reduced model
plot(e,dnorm(e,0,reduced.sigmasq.mle),type="l",main="reduced",ylab="density",ylim=c(0,0.8), xlab="residual")
points(residuals(lm.reduced),rep(0,10),pch=4)
for(i in 1:10){
	lines(residuals(lm.reduced), dnorm(residuals(lm.reduced),0,reduced.sigmasq.mle),type="h",lty=2)
}

# Full model and residuals
plot(x,y, pch=16)
abline(lm.full)
for(i in 1:10){
	lines(c(x[i],x[i]), c(y[i],fitted(lm.full)[i]), lty=2)
}
full.sigmasq.mle<-mean(residuals(lm.full)^2)

# Estimated error distribution for full model
plot(e,dnorm(e,0,full.sigmasq.mle),type="l",main="full",ylab="density",ylim=c(0,0.8), xlab="residual")
points(residuals(lm.full),rep(0,10),pch=4)
for(i in 1:10){
	lines(residuals(lm.full), dnorm(residuals(lm.full),0,full.sigmasq.mle),type="h",lty=2)
}


#########################
##### Bootstrapping #####
#########################

# The bootstrap hypothesis test 

# Choose how many test-statistics to sample, and set up a placeholder
N<-200
boot.test.stats<-rep(0,N)

for(i in 1:N){
	
	# Generate new data from the reduced model
	y.new<-unlist(simulate(lm.reduced))
	
	# Fit both models to the new data
	lm.full.new<-lm(y.new~x)
	lm.reduced.new<-lm(y.new~1)
	
	# Calculate the test statistic for the new models
	boot.test.stats[i]<- -2*(logLik(lm.reduced.new) - logLik(lm.full.new))
	par(mfrow=c(1,2))
	plot(x,y.new,ylim=c(-1,4), pch=16)
	abline(lm.reduced.new, col="red",lty=2, lwd=2)
	abline(lm.full.new,col="blue", lwd=2)
	
	hist(boot.test.stats[1:i],prob=T,xlim=c(0,15),main="")
	points(boot.test.stats[i],0,pch=4,col="blue", cex=2)
	abline(v=obs.test.statistic, lwd=2)
	print(c(logLik(lm.reduced.new) , logLik(lm.full.new), -2*(logLik(lm.reduced.new) - logLik(lm.full.new))))
	if(i<10){readline()}
}

# Plot test statistics and chis-sq approximation
par(mfrow=c(1,1))
hist(boot.test.stats,prob=T)
curve(dchisq(x,df=1),from=0,to=20,add=T)
points(obs.test.statistic,0,pch=4,cex=4,col="red")
abline(v = obs.test.statistic, col=2)
# Get the p-value - how many simulated test statstics are larger than the observed one?
# How does the p-value compare with the anova test p-value in step 3?
mean(boot.test.stats>obs.test.statistic)

anova(lm.reduced, lm.full)
#############################################################
# Bootstrapping for confidence intervals for the full model #
#############################################################

# (Computationally, better to do this simultaneously with the hypothesis test)

# Data will be generated from this model
summary(lm.full)

# Choose how many test-statistics to sample, and set up a placeholder
N<-200
boot.parameters<-matrix(0,N,3)

for(i in 1:N){
	
	# Generate new data from the full model
	y.new<-unlist(simulate(lm.full))
	
	# Fit full model to the new data
	lm.full.new<-lm(y.new~x)
	
	# Extract parameter estimates for full model
	boot.parameters[i,1:2]<-lm.full.new$coefficients
	boot.parameters[i,3]<-summary(lm.full.new)$sigma
	par(mfrow=c(2,2))
	plot(x,y.new, ylim=c(-2,5),pch=16)
	abline(lm.full,lty=2,col="red",lwd=2)
	abline(lm.full.new,col="blue")
	
	hist(boot.parameters[1:i,1],prob=T,main="alpha", xlim=c(-2,2))
	abline(v=-.169,col="red",lwd=2)
	points(boot.parameters[i,1],0,pch=4,col="blue")
	
	hist(boot.parameters[1:i,2],prob=T,main="beta", xlim=c(-0.05,0.6))
	abline(v=.255,col="red",lwd=2)
	points(boot.parameters[i,2],0,pch=4,col="blue")
	
	hist(boot.parameters[1:i,3],prob=T,main="sigma", xlim=c(0,2))
	abline(v=0.809,col="red",lwd=2)
	points(boot.parameters[i,3],0,pch=4,col="blue")
	if(i<20){readline()}
	print(i)
	}
	
# Calculate sample confidence intervals
apply(boot.parameters,2,quantile,probs=c(0.025,0.975))

# Alternative CIs for intercept and slope
(boot.m<-apply(boot.parameters[,1:2],2,mean))
(boot.se<-apply(boot.parameters[,1:2],2,sd))
matrix(c(boot.m-2*boot.se,boot.m+2*boot.se),2,2, byrow=T)

# Compare with usual confidence intervals
confint(lm.full)

# boot.m and boot.se are very similar to estimates and standard errors in 
summary(lm.full)$coefficients
# Bootstrap interval is narrower because correct interval uses +/-
qt(0.975,n-2)
# estimated standard errors, rather than +/-2 estimated standard errors. 
# Boostrap CI hasn't accounted for uncertainty in sigma.

# Interval for sigma^2 based on chi-sq distribution
(sigma.CI<-c( ((n-2)/qchisq(0.975,n-2))^0.5*summary(lm.full)$sigma, ((n-2)/qchisq(0.025,n-2))^0.5*summary(lm.full)$sigma))

#########################
# Section 4.2.3 Example #
#########################

library(lme4)
attach(ergoStool) 

(fm.full<-lmer(effort~Type + (1|Subject),ergoStool,REML=F))
logLik(fm.full)

(fm.reduced<-lmer(effort~1 + (1|Subject),ergoStool,REML=F))
logLik(fm.reduced)

# GLRT with chi-squared approximation (not recommended!)
# We need to use ordinary maximum likelihood and not REML here
# The REML likelihood has the term |X^T V X| which is not invariant 
# to reparameterising the fixed effects. To see this, try

(fm1a<-lmer(effort~Type + (1|Subject),ergoStool))
(fm1b<-lmer(effort~Type + (1|Subject), contrasts=list(Type=contr.sum) ,ergoStool))
logLik(fm1a, REML=T)
logLik(fm1b, REML=T)


obs.test.stat<- - 2*(logLik(fm.reduced)-logLik(fm.full)) 
attributes(obs.test.stat)<-NULL

obs.test.stat
1-pchisq(obs.test.stat,3)

# Quick way to do the GLRT
anova(fm.reduced,fm.full)

# If you try to do the GLRT with REML, R will refit the models with ML:
fm.full.reml<-lmer(effort~Type + (1|Subject),ergoStool)
fm.reduced.reml<-lmer(effort~1 + (1|Subject),ergoStool)
anova(fm.full.reml,fm.reduced.reml )

# The anova command on a single command will give you 
# an F-statistic corresponding to an F-test, but no p-value
# For unballanced datasets, cannot work out exact distribution
# of F under null hypothesis.
anova(fm.full)

# Using bootstrapping
N<-100
boot.test.stats<-rep(0,N)

for(i in 1:N){
  new.y<-unlist(simulate(fm.reduced))
  fm.reduced.new<-lmer(new.y~1 + (1|Subject), REML=F)
  fm.full.new<-lmer(new.y~Type + (1|Subject), REML=F)
  boot.test.stats[i]<- -2*(logLik(fm.reduced.new)-logLik(fm.full.new))
  par(mfrow=c(1,2))
  plot(as.numeric(Type),new.y,ylim=c(5,20),pch=16,xaxp=c(1,4,3))
  title(paste("iteration:",i))
  hist(boot.test.stats[1:i],xlim=c(1,40),main="GLRT")
  points(boot.test.stats[i],0,pch=4,col="red", cex=2)
  abline(v=obs.test.stat, lwd=2)
  if(i<10){readline()}
}
# Plot test statistics and chis-sq approximation
par(mfrow=c(1,1))
hist(boot.test.stats,prob=T,xlim=c(0,40),ylim=c(0,0.25))
curve(dchisq(x,df=3),from=0,to=40,add=T,lwd=2)
points(obs.test.stat,0,pch=4,col="red", cex=2)

# Get the p-value - how many simulated test statstics are larger than the observed one?
# How does the p-value compare with the anova test p-value in step 3?
mean(boot.test.stats>obs.test.stat)

#### Bootstrap intervals
fm1<-lmer(effort~Type-1 + (1|Subject),ergoStool)
N<-100 # Needs to be larger, but will take longer to run
fixed.effects<-matrix(0,N,4)
variance.components<-matrix(0,N,2)

for(i in 1:N){
  new.y<- unlist(simulate(fm.full))
  new.fm1 <- lmer(new.y~Type-1 + (1|Subject))
  fixed.effects[i,] <- fixef(new.fm1)
  vc <- VarCorr(new.fm1)
  variance.components[i,]<-c(unlist( lapply(vc, diag) ),attr(vc,"sc")^2)
  
  par(mfrow=c(2,4))
  plot(as.numeric(Type),new.y,ylim=c(5,20), pch=16)
  title(paste("iteration:",i))
  hist(fixed.effects[1:i,1],main="")
  title(expression(mu))
  points(fixed.effects[i,1],0,pch=4,col="red", cex=2)
  
  hist(fixed.effects[1:i,2],main="")
  title(expression(beta[2]))
  points(fixed.effects[i,2],0,pch=4,col="red", cex=2)
  
  hist(fixed.effects[1:i,3],main="")
  title(expression(beta[3]))
  points(fixed.effects[i,3],0,pch=4,col="red", cex=2)
  
  hist(fixed.effects[1:i,4],main="")
  title(expression(beta[4]))
  points(fixed.effects[i,4],0,pch=4,col="red", cex=2)
  
  hist(variance.components[1:i,1]^0.5,main="")
  title(expression(sigma[b]))
  points(variance.components[i,1]^0.5,0,pch=4,col="red", cex=2)
  
  hist(variance.components[1:i,2]^0.5,main="")
  title(expression(sigma))
  points(variance.components[i,2]^0.5,0,pch=4,col="red", cex=2)
  # Element [i,1] is random effect variance
  # Element [i,2] is residual variance
  if(i<10){readline()}}


# Extract bootstrap confidence intervals

apply(fixed.effects,2,quantile,probs=c(0.025,0.975))

apply(variance.components^0.5,2,quantile,probs=c(0.025,0.975))

#################
# Section 4.2.5 #
#################

(fm1<-lmer(score~Machine-1+(1|Worker/Machine),data=Machines))
(fm2<-lmer(score~Machine-1+(1|Worker),data=Machines))

# Test to see if the worker-machine interaction random effects
# are significant
anova(fm1, fm2)

# Can see this in the interaction plot
attach(Machines)
interaction.plot(Machine,Worker,score)


# Compare the more complex model (3.13) against (3.4)
(fm3<-lmer(score~Machine-1+(Machine-1|Worker),data=Machines))
anova(fm3, fm1)

