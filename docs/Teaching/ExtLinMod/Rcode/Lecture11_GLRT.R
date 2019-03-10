library(lme4)
library(lattice)
load('MAS473.RData')


###################################
##
## Compare fixed effects structures
##
###################################


attach(ergoStool)
ergoStool




fm.full<-lmer(effort~Type + (1|Subject),ergoStool,REML=F)
logLik(fm.full)

fm.reduced<-lmer(effort~1 + (1|Subject),ergoStool,REML=F)
logLik(fm.reduced)
print('NOTE: we must use maximum likelihood, not REML, i.e. setREML=FALSE')



#GLRT with chi-squared approximation
(obs.test.stat<- - 2*(logLik(fm.reduced)-logLik(fm.full)) )
attributes(obs.test.stat)<-NULL

obs.test.stat
1-pchisq(obs.test.stat,3)

# Quick way to do the GLRT
anova(fm.reduced,fm.full)

# Note that lmer package will not report p-values!
anova(fm.full)

# Using bootstrapping
N<-100
boot.test.stats<-rep(0,N)
#wait<-"yes"
for(i in 1:N){
  print(i)
  new.y<-unlist(simulate(fm.reduced))
  fm.reduced.new<-lmer(new.y~1 + (1|Subject), REML=F)
  fm.full.new<-lmer(new.y~Type + (1|Subject), REML=F)
  boot.test.stats[i]<- -2*(logLik(fm.reduced.new)-logLik(fm.full.new))
  # par(mfrow=c(1,2))
#   plot(Type,new.y,ylim=c(5,20))
#   title(paste("iteration:",i))
#   hist(boot.test.stats[1:i],xlim=c(1,40),main="GLRT")
#   points(boot.test.stats[i],0,pch=4,col="red")
#   abline(v=obs.test.stat)
#   if(wait!="r"){wait<-readline()}
}
# Plot test statistics and chis-sq approximation
hist(boot.test.stats,prob=T,xlim=c(0,40),ylim=c(0,0.25))
curve(dchisq(x,df=3),from=0,to=40,add=T)
points(obs.test.stat,0,pch=4,col="red", lwd=3)
abline(v=obs.test.stat, col=2, lwd=2)

# Get the p-value - how many simulated test statstics are larger than the observed one?
# How does the p-value compare with the anova test p-value in step 3?
mean(boot.test.stats>obs.test.stat)

#### Bootstrap intervals

fm1<-lmer(effort~Type-1 + (1|Subject),ergoStool,contrasts=list(Type=contr.sum))
N<-100

summary(fm1)
# gives standard errors
8.556 +2*0.5760
8.556 -2*0.5760




fixed.effects<-matrix(0,N,4)
variance.components<-matrix(0,N,2)

for(i in 1:N){
  print(i)
  new.y<- unlist(simulate(fm1))
  new.fm1 <- lmer(new.y~Type-1 + (1|Subject),contrasts=list(Type=contr.sum))
  fixed.effects[i,] <- fixef(new.fm1)
  vc <- VarCorr(new.fm1)
  variance.components[i,]<-c(unlist( lapply(vc, diag) ),attr(vc,"sc")^2)
  # Element [i,1] is random effect variance
  # Element [i,2] is residual variance
  
  }

par(mfrow=c(2,4))
plot(Type,new.y,ylim=c(5,20))
title(paste("iteration:",i))

hist(fixed.effects[1:i,1],main="")
title(expression(beta[1]))

hist(fixed.effects[1:i,2],main="")
title(expression(beta[2]))

hist(fixed.effects[1:i,3],main="")
title(expression(beta[3]))

hist(fixed.effects[1:i,4],main="")
title(expression(beta[4]))

hist(variance.components[1:i,1]^0.5,main="")
title(expression(sigma[b]))

hist(variance.components[1:i,2]^0.5,main="")
title(expression(sigma))

# Extract bootstrap confidence intervals

apply(fixed.effects,2,quantile,probs=c(0.025,0.975))

apply(variance.components^0.5,2,quantile,probs=c(0.025,0.975))


##########################################################################################
#
# Comparing random effects structures
#
########################################################################################


fm.full <- lmer(score~Machine + (1|Worker/Machine), data=Machines, REML=FALSE)
fm.reduced <- lmer(score~Machine + (1|Worker), data=Machines, REML=FALSE)
anova(fm.reduced, fm.full)

N<-100
boot.test.stats<-rep(0,N)
attach(Machines)
#wait<-"yes"
for(i in 1:N){
  print(i)
  new.score<-unlist(simulate(fm.reduced))
  fm.reduced.new<-lmer(new.score~ Machine + (1|Worker), REML=F)
  fm.full.new<-lmer(new.score~Machine + (1|Worker/Machine), REML=F)
  boot.test.stats[i]<- -2*(logLik(fm.reduced.new)-logLik(fm.full.new))
}

par(mfrow=c(1,1))
hist(boot.test.stats,prob=T,xlim=c(0,40),ylim=c(0,0.25))
curve(dchisq(x,df=3),from=0,to=40,add=T)
points(obs.test.stat,0,pch=4,col="red")
abline(v=obs.test.stat, col=2, lwd=2)

# Get the p-value - how many simulated test statstics are larger than the observed one?
# How does the p-value compare with the anova test p-value in step 3?
mean(boot.test.stats>obs.test.stat)
anova(fm.reduced, fm.full)

####################
#
fm.full <- lmer(score~Machine-1 + (Machine-1|Worker), data=Machines, REML=FALSE)
fm.reduced <- lmer(score~Machine-1 + (1|Worker), data=Machines, REML=FALSE)
anova(fm.full, fm.reduced)





###################################
#
# We can also compare score~Machine-1 + (Machine-1|Worker)
# to score~Machine + (1|Worker/Machine)
# because this second model is nested within the first model
#
# The second model has a covariance matrix specified by sigma1 and sigma2.
# The first model has a covariance matrix specified by three variance terms plus three correlations, so is more general.
#
###################################

fm.full2 <- lmer(score~Machine -1+ (1|Worker/Machine), data=Machines, REML=FALSE)
anova(fm.full2, fm.full)

(obs.test.stat<- - 2*(logLik(fm.full2)-logLik(fm.full)) )
attributes(obs.test.stat)<-NULL



N<-100
boot.test.stats<-rep(0,N)
attach(Machines)
#wait<-"yes"
for(i in 1:N){
  print(i)
  new.score<-unlist(simulate(fm.full2))
  fm.full2.new<-lmer(new.score~Machine -1+ (1|Worker/Machine), REML=FALSE)
  fm.full.new<-lmer(score~Machine-1 + (Machine-1|Worker), REML=FALSE)
  boot.test.stats[i]<- -2*(logLik(fm.full2.new)-logLik(fm.full.new))
}


hist(boot.test.stats,prob=T,xlim=c(0,40),ylim=c(0,0.25))
curve(dchisq(x,df=3),from=0,to=40,add=T)
points(obs.test.stat,0,pch=4,col="red")
abline(v=obs.test.stat, col=2, lwd=2)

# Get the p-value - how many simulated test statstics are larger than the observed one?
# How does the p-value compare with the anova test p-value in step 3?
mean(boot.test.stats>obs.test.stat)


