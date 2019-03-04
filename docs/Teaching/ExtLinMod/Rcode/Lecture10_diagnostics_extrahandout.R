## ---- messages=F---------------------------------------------------------
set.seed(1)
library(ggplot2)
library(lme4)

ngroup <- 100
groupsize <- 200
group <- gl(ngroup,groupsize)
b_i <- rnorm(ngroup) # random effects
epsilon <- rnorm(ngroup*groupsize)
y <- 1+b_i[group]+epsilon
data = data.frame(y, group = group)

## ---- cache=T------------------------------------------------------------
fm1<- lmer(y~1+(1|group), data = data)

## ------------------------------------------------------------------------
qplot(fitted(fm1),resid(fm1))
qqnorm(resid(fm1), main='QQ plot for level 1 residuals')
qqline(resid(fm1))

## ------------------------------------------------------------------------
fitted.level0<-fm1@pp$X %*% fixef(fm1)
resid.level0<-y-fitted.level0

qqnorm(resid.level0,main='QQ plot for level 0 residuals')
qqline(resid.level0)

## ------------------------------------------------------------------------
# Assess general fit of model
plot(fitted(fm1),data$y)
abline(0,1)

## ---- cache=TRUE---------------------------------------------------------
b_i <- rnorm(ngroup) # random effects
epsilon <- rt(ngroup*groupsize, df=3)
y <- 1+b_i[group]+epsilon
data = data.frame(y, group = group)
fm1<- lmer(y~1+(1|group), data = data)
qqnorm(resid(fm1), main='QQ plot for level 1 residuals')
qqline(resid(fm1))
fitted.level0<-fm1@pp$X %*% fixef(fm1)
resid.level0<-y-fitted.level0
qqnorm(resid.level0, main='QQ plot for level 0 residuals')
qqline(resid.level0)

## ---- cache=TRUE---------------------------------------------------------
b_i <- rt(ngroup, df=2) # random effects
epsilon <- rnorm(ngroup*groupsize)
y <- 1+b_i[group]+epsilon
data = data.frame(y, group = group)
fm1<- lmer(y~1+(1|group), data = data)
qqnorm(resid(fm1), main='QQ plot for level 1 residuals')
qqline(resid(fm1))
fitted.level0<-fm1@pp$X %*% fixef(fm1)
resid.level0<-y-fitted.level0
qqnorm(resid.level0, main='QQ plot for level 0 residuals')
qqline(resid.level0)

