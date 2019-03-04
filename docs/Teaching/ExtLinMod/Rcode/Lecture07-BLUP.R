

###################
# BLUP for raildata
library(lme4)
library(magic)
load('MAS473.RData')
fm1<-lmer(travel~1+(1|Rail),raildata)
ranef(fm1)

V1<-matrix(615.31,3,3)
diag(V1)<-diag(V1)+16.17

# Estimated variance matrix for Y
V.Y<-adiag(V1,V1,V1,V1,V1,V1) 

# Estimated covariance between b and Y 
V.temp<-V.Y-diag(16.17,18)
(cov.b.Y<-V.temp[c(1,4,7,10,13,16),])


(Y<-matrix(raildata$travel,18,1))

# Estimated fixed effects
(beta.hat<-matrix(fixef(fm1),1,1))
X<-matrix(1,18,1)

# Predicted random effect
cov.b.Y %*% solve(V.Y)%*%(Y-X%*%beta.hat)

ranef(fm1)
