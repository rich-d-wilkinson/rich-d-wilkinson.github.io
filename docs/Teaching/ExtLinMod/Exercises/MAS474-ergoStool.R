library(lme4)
library(lattice)

load('MAS473.RData')
attach(ergoStool)
ergoStool

#### Section 2.5.2 inspecting the data
str(ergoStool)
head(ergoStool)
tail(ergoStool)
xtabs(~ Type + Subject, ergoStool)
xtabs(~ Type, ergoStool)

matrix(by(effort,list(Type,Subject),mean),4,9)


#### Plotting the data
xyplot(effort~Type|Subject,data=ergoStool)
library(ggplot2)
qplot(Type, effort, facets= ~Subject, data = ergoStool)

dotplot(reorder(Subject,effort) ~ effort, data=ergoStool, groups = Type, type = c("p"),pch=Type, par.settings=simpleTheme(pch=Type), xlab = "Effort to arise (Borg scale)", auto.key = list(columns=4))
qplot(effort, reorder(Subject, effort), color = Type, shape=Type, xlab='Effort to arise (Borg scale)')


# the argument reorder(Subject,effort) arranges the groups ("Subject") in order of increasing mean effort. 

plot(Type,effort, xlab="Type", ylab="effort")
qplot(Type, effort, geom='boxplot')

plot(reorder(Subject,effort),effort, xlab="Subject", ylab="effort") 
qplot(reorder(Subject, effort), effort, geom='boxplot', xlab='Subject')



## What?

#### 2.5.3 Fitting the model

(fm1<-lmer(effort~Type -1 + (1|Subject),ergoStool))

summary(fm1)


# Compare fixed effects estimates with

by(effort,Type,mean)
mean(effort[Type=="T1"])


# Standard errors for the fixed effects estimates
(1.775/9+1.211/9)^0.5

# Compare with ordinary linear model
lm1<-lm(effort~Type-1+Subject,contrasts=list(Subject=contr.sum))
summary(lm1)

# Correlation between fixed effects estimates
1.775/(1.775+1.211)

# Task 7

(fm2<-lmer(effort~Type + (1|Subject), contrasts=list(Type=contr.sum), data=ergoStool))

# Different maximised REML values:
logLik(fm1,REML=T)
logLik(fm2,REML=T)

# Ordinary maximised likelihood value unchanged:
logLik(fm1,REML=F)
logLik(fm2,REML=F)




#### Section 4.1 assessing the fitted model

# Expected values of random effects
ranef(fm1)

# Check assumption that random effects are normally distributed
qqnorm(unlist(ranef(fm1)))
abline(0,1.332)

# Plot level 1 fitted values against residuals
fitted(fm1)
resid(fm1)
plot(fitted(fm1),resid(fm1))
xyplot(resid(fm1)~fitted(fm1)|Type) # Residual plots by group

# Check assumption that errors are normally distributed
qqnorm(resid(fm1))
abline(0,1.1003)

# Plot level 0 fitted values against residuals
fitted.level0<-fm1@pp$X %*% fixef(fm1)
resid.level0<-effort-fitted.level0
plot(fitted.level0,resid.level0)

# Can jitter the points to make easier to see
plot(fitted.level0+rnorm(36,0,0.025),resid.level0)


# Assess general fit of model
plot(fitted(fm1),ergoStool$effort)
abline(0,1)

