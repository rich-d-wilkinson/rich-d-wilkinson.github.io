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
library(ggplot2)
qplot(Type, effort, facets= ~Subject, data = ergoStool)
qplot(effort, reorder(Subject, effort), color = Type, shape=Type, xlab='Effort to arise (Borg scale)')


# the argument reorder(Subject,effort) arranges the groups ("Subject") in order of increasing mean effort. 
qplot(Type, effort, geom='boxplot')
qplot(reorder(Subject, effort), effort, geom='boxplot', xlab='Subject')


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

