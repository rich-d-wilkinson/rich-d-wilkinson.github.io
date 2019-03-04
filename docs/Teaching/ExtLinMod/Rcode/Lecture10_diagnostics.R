library(lme4)
load('MAS473.RData')
attach(ergoStool)
library(ggplot2)
ergoStool

# Start with a plot of the data, to look for outliers
qplot(Type, effort, facets=~Subject,data=ergoStool)

fm1<-lmer(effort~Type -1 + (1|Subject),ergoStool)

#### 3.5.6 assessing the fitted model

# Expected values of random effects
ranef(fm1)

plot(ranef(fm1)$Subject, pch=16)

# Check assumption that random effects are normally distributed
qqnorm(unlist(ranef(fm1)))
abline(0,1.332) # Reference gradient is estimated standard deviation of random effect

# For background, R uses these quantiles for the QQ plot
# Do ?ppoints for more details
(qu<-(1:9 - 3/8)/(9 + (1-3/8)-3/8))
z<-qnorm(qu)
plot(z,sort(unlist(ranef(fm1)$Subject)) )

# Fitted values
fitted(fm1)

# Residuals (not standardaised)
resid(fm1)

# Plot level 1 fitted values against residuals
qplot(fitted(fm1),resid(fm1))
qplot(fitted(fm1), resid(fm1), facets=~Type, data=ergoStool) # Residual plots by stool type
qplot(fitted(fm1), resid(fm1), facets=~Subject, data=ergoStool) # Residual plots by subject


# Check assumption that errors are normally distributed
qqnorm(resid(fm1))
#qqline(resid(fm1))
abline(0,1.1003)

# Plot level 0 fitted values against residuals
(fitted.level0<-fm1@pp$X %*% fixef(fm1))
resid.level0<-effort-fitted.level0
qplot(fitted.level0,resid.level0)

# Can jitter the points to make easier to see
qplot(fitted.level0,resid.level0)+geom_jitter() 


# Assess general fit of model
plot(fitted(fm1),ergoStool$effort)
abline(0,1)


################################
# Section 4.1.1 Oxides example #
################################

# Plot the data first to check for outlying lots, wafers within lots, sites within wafers

attach(Oxide)
qplot(Wafer, Thickness, facets=~Lot)



fm1<-lmer(Thickness~1+(1|Lot/Wafer),data=Oxide)
summary(fm1)

# Estimated random effects
ranef(fm1)

# Check model assumptions

qqnorm(unlist(ranef(fm1)$`Wafer:Lot`),ylab="Wafer-within-lot-level random effects")
abline(0,5.9891)
qqline(unlist(ranef(fm1)$`Wafer:Lot`))

qqnorm(unlist(ranef(fm1)$Lot),ylab="Lot-level random effects")
abline(0,11.3967)
# Doesn't look so good, but sample size is small

qqnorm(resid(fm1))
abline(0, 3.5453)
qqline(resid(fm1))

plot(fitted(fm1),residuals(fm1))

attach(Oxide)
qplot(fitted(fm1), residuals(fm1), facets=~Lot)

# Why isn't the following plot helpful?
qplot(fitted(fm1),residuals(fm1), facets=~Wafer)

# Assess overall model fit
plot(fitted(fm1),Thickness)
abline(0,1)

# Three levels of fitted values
# In matrix notation (Section 2.6), can extract Z with

(Zt<-fm1@pp$Zt) # this is transpose of Z

# Convert into normal matrix form, and take transpose
Z<-as.matrix(t(Zt))

# Put the predicted random effects into single column vector
b.hat<-as.matrix(rbind(ranef(fm1)$`Wafer:Lot`,ranef(fm1)$Lot ))

# Inner most level is 'level 2'
# \hat{y_{ijk}} = \hat{beta} + \hat{b_i} + \hat{b_{ij}}

# Calculate manually
fitted.level2<-as.matrix(fm1@pp$X %*% fixef(fm1) + Z%*%b.hat)
# inner most level also given by the fitted command
fitted(fm1)
# Check they are the same
cbind(fitted.level2, fitted(fm1))

# Level 1
# \hat{y_{ijk}} = \hat{beta} + \hat{b_i} 

fitted.level1<-as.matrix(fm1@pp$X %*% fixef(fm1) + Z[,25:32]%*%b.hat[25:32,1])
resid.level1<-Thickness-fitted.level1
plot(fitted.level1, resid.level1)

# Level 0 (not very informative, as only a single fixed effect parameter)
# \hat{y_{ijk}} = \hat{beta} 
(fitted.level0<-as.matrix(fm1@pp$X %*% fixef(fm1) ))

