library(lme4)
load('MAS473.RData')


#######################################
####### Chapter 4 Section 3.2 #########
attach(Machines)


head(Machines)

library(ggplot2)
qplot(Machine, score, geom='boxplot')
(q<-qplot(Worker, score, color=Machine))
q+stat_summary(fun.y = mean, geom = "line", aes(group = Machine, lty=Machine))



(g <- qplot(Machine, score, color=Worker))
g + stat_summary(fun.y = mean, geom = "line", aes(group = Worker, lty=Worker))



# But for once, these types of plots are easier using the built in interaction.plot command
interaction.plot(Machine,Worker,score)

# Plot suggests significant machine effects, with C better than B better than A
# Suggests significant worker random effects, with worker 3 the best
# Also suggests interaction, looking at worker 6

(fm1<-lmer(score~Machine-1+(1|Worker/Machine),data=Machines))
summary(fm1)

# Now consider the following
(fm3<-lmer(score~Machine-1+(Machine-1|Worker),data=Machines))
summary(fm3)

# and compare
ranef(fm1)
ranef(fm3)
# Note that model fm1, we can consider a single random effect term as
# d_{ij} = b_i + b_{ij}
# Predicted values for this new term d_{ij} given by
matrix(unlist(ranef(fm1)$`Machine:Worker`),6,3) + matrix(unlist(ranef(fm1)$Worker),6,3)

# Compare these with the random effects in fm3:
ranef(fm3)



# Look again at interaction plot
interaction.plot(Machine,Worker,score)
# Can see why the MachineB variance is higher

# Can also see why error variance is small
machinecolours<-rep(c("red","blue","green"),each=18)
#plot(as.numeric(Worker), score, col=machinecolours, pch=16)
qplot(Worker, score, col=Machine, cex=16)



# Model fm1
# Y_{ijk} = beta_j + b_i + b_{ij} + epsilon_{ijk}
# with b_i~N(0,sigma^2_1), b_{ij}~N(0,sigma^2_2), epsilon_{ijk}~N(0,sigma^2)

# Var(Y_{ijk}) = sigma^2_1 + sigma^2_2 + sigma^2 for all i,j,k
# Covariance between observations on the same worker on *different* machines is sigma^2_1, for any pair of machines.

# Model fm3
# Y_{ijk} = beta_j + b_{ij} + epsilon_{ijk}
# with b_i1~N(0,sigma^2_1), b_{i2}~N(0,sigma^2_2),  b_{i3}~N(0,sigma^2_3), epsilon_{ijk}~N(0,sigma^2)
# so variance changes for each machine:
# Var(Y_{i1k}) = sigma^2_1 + sigma^2
# Var(Y_{i2k}) = sigma^2_2 + sigma^2
# Var(Y_{i3k}) = sigma^2_3 + sigma^2
# Covariance between observations on the same worker on *different* machines is allowed to vary for any pair of machines.

# Model fm1 equivalent to fm3, with the constraints that variance is constant for each machine, and covariance of worker effects on any pair of machines is fixed across all pairs. 

# We will do formal testing later, but what would we look for in 
# an exploratory data analysis to choose between the models?

# Calculate means for each combination of worker and machine
cellmeans<-matrix(by(score,list(Worker,Machine),mean),6,3)
# Note Worker is an ordered factor, ordered by increasing mean, so 1st row in cellmeans corresponds to worker 6,
# second row corresponds to worker 2 and so on

# Calculate variances of worker means for each machine
apply(cellmeans,2,var)
# variance for machine B is higher, suggesting unequal variance model fm3 might be suitable
# Calculate correlations of worker means between each machine
cor(cellmeans)
# unequal off-diagonal correlations would suggest trying fm3
# However, sample size (6) for Workers is small, so not strong evidence to support fm3.

# Note that these sample variances and correlations are very similar to the parameter estimates
summary(fm3)

