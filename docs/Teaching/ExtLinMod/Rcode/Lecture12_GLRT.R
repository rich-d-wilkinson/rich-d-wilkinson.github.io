
# 1. Generate some data. 
set.seed(1)
n<-10 # the sample size
x<-seq(from=1,to=10,length=n)
y<-0.2*x+rnorm(n) # use this for full model to be true
#y<-rnorm(10) # use this for reduced model to be true
plot(x,y)

# 2. Fit the two linear models
lm.full<-lm(y~x)
lm.reduced<-lm(y~1)


# 3. Usual F-test to test H_0: lm.reduced is the true model
anova(lm.reduced,lm.full)

# 4. Compute maximum likelihood estimates for parameters in reduced model
(mu.mle<-mean(y))
(sigmasq.mle<-mean((y-mu.mle)^2))

# compare estimates with those in the fitted model - will they be the same?
lm.reduced$coefficients[1]
summary(lm.reduced)$sigma ^2
summary(lm.reduced)$sigma ^2*9/10


# 5. Evaluate the log-likelihood for the full model
sum(log(dnorm(y,mu.mle,sigmasq.mle^0.5)))

# direct way to get the log likelihood in R. 
logLik(lm.reduced)

# Compute log likelihood for full model
logLik(lm.full)

# 6. Perform the GLRT with chi-sq approximation
(obs.test.statistic<- -2*(logLik(lm.reduced) - logLik(lm.full)))
attributes(obs.test.statistic)<-NULL # get rid of the 'clutter'!

# Calculate the p-value. How does it compare with the p-value in step 3?
1-pchisq(obs.test.statistic,1)


#########################
##### Bootstrapping #####
#########################

# 7. The bootstrap hypothesis test 

# Choose how many test-statistics to sample, and set up a placeholder
N<-10
boot.test.stats<-rep(0,N)
wait<-"yes"
for(i in 1:N){
	
	# Generate new data from the reduced model
	y.new<-unlist(simulate(lm.reduced))
	
	# Fit both models to the new data
	lm.full.new<-lm(y.new~x)
	lm.reduced.new<-lm(y.new~1)
	
	# Calculate the test statistic for the new models
	boot.test.stats[i]<- -2*(logLik(lm.reduced.new) - logLik(lm.full.new))
	par(mfrow=c(1,2))
	plot(x,y.new)
	abline(lm.reduced.new)
	abline(lm.full.new,col="red")
	
	hist(boot.test.stats[1:i],prob=T,xlim=c(0,15),main="")
	points(boot.test.stats[i],0,pch=4,col="red")
	abline(v=obs.test.statistic)
	print(c(logLik(lm.reduced.new) , logLik(lm.full.new)))
	if(wait!="r"){wait<-readline()}
}

# Plot test statistics and chis-sq approximation
hist(boot.test.stats,prob=T)
curve(dchisq(x,df=1),from=0,to=20,add=T)
points(obs.test.statistic,0,pch=4,col="red")

# Get the p-value - how many simulated test statstics are larger than the observed one?
# How does the p-value compare with the anova test p-value in step 3?
mean(boot.test.stats>obs.test.statistic)



# 8. Bootstrapping for confidence intervals for the full model
# (Computationally, better to do this simultaneously with the hypothesis test)

# Choose how many test-statistics to sample, and set up a placeholder
N<-100
boot.parameters<-matrix(0,N,3)
wait<-"no"
for(i in 1:N){
	
	# Generate new data from the full model
	y.new<-unlist(simulate(lm.full))
	
	# Fit full model to the new data
	lm.full.new<-lm(y.new~x)
	
	# Extract parameter estimates for full model
	boot.parameters[i,1:2]<-lm.full.new$coefficients
	boot.parameters[i,3]<-summary(lm.full.new)$sigma
	par(mfrow=c(2,2))
	plot(x,y.new)
	abline(lm.full)
	hist(boot.parameters[1:i,1],prob=T,main="alpha")
	points(boot.parameters[i,1],0,pch=4,col="red")
	hist(boot.parameters[1:i,2],prob=T,main="beta")
	points(boot.parameters[i,2],0,pch=4,col="red")
	hist(boot.parameters[1:i,3],prob=T,main="sigma")
	points(boot.parameters[i,3],0,pch=4,col="red")
	if(wait!="r"){wait<-readline()}
	}
	
# 8.1 Calculate sample confidence intervals
apply(boot.parameters,2,quantile,probs=c(0.025,0.975))

# 8.2 Alternative CIs for intercept and slope
(boot.m<-apply(boot.parameters[,1:2],2,mean))
(boot.se<-apply(boot.parameters[,1:2],2,sd))
matrix(c(boot.m-2*boot.se,boot.m+2*boot.se),2,2)

# 8.3 Compare with usual confidence intervals
confint(lm.full)

# boot.m and boot.se are very similar to estimates and standard errors in 
summary(lm.full)$coefficients
# Bootstrap interval is narrower because correct interval uses +/-
qt(0.975,n-2)
# estimated standard errors, rather than +/-2 estimated standard errors. 
# Boostrap CI hasn't accounted for uncertainty in sigma.

# 8.4 Interval for sigma^2 based on chi-sq distribution
(sigma.CI<-c( ((n-2)/qchisq(0.975,n-2))^0.5*summary(lm.full)$sigma, ((n-2)/qchisq(0.025,n-2))^0.5*summary(lm.full)$sigma))