set.seed(1)

X1 <- rnorm(30,0,3)
X2 <- 0.6*X1 + rnorm(30,0,3)
Y = 0.3*X1 + 2*X2 + rnorm(30,0,1)

(fit <- lm(Y~X1+X2))
summary(fit)

data.full <- data.frame(Y=Y, X1=X1, X2=X2)
miss.rate <- 0.2
X1[which(as.logical(rbinom(30,1,miss.rate)))] <- NA
X2[which(as.logical(rbinom(30,1,miss.rate)))] <- NA


data.mis <- data.frame(Y=Y, X1=X1, X2=X2)
data.mis
(fit.completecase <- lm(Y~X1+X2, data.mis))

############################
# complete-case analysis

library(mice)
(data.cc <- cc(data.mis))


lm(Y~X1+X2, data=data.cc)
lm(Y~X1+X2, data=data.mis)
# R does CC by default

cov(data.mis[,2:3])
cov(data.mis[,2:3], use='complete.obs')

cov(data.cc[,2:3])


colMeans(data.mis)
colMeans(data.cc)
colMeans(data.mis, na.rm=T) # why different? - this is doing available- case

##
## Mean imputation

MeanImpute <- function(data){
  nc = dim(data)[2] # number of columns
  data.impute <- data
  
  for(ii in 1:nc){ # loop through columns
    data.impute[is.na(data[,ii]),ii] <- mean(data[,ii], na.rm=T)
  }
  return(data.impute)
} 

(data.meanImp <- MeanImpute(data.mis))
colMeans(data.meanImp)
colMeans(data.mis, na.rm=T) # no change -why

# mean imputation shrinks variances and standard error estimates
apply(data.mis, 2, function(x) var(x,na.rm=T))
apply(data.meanImp, 2, function(x) var(x))

cov(data.mis[,2:3], use='complete.obs')
cov(data.meanImp[,2:3])


summary(lm(Y~X1+X2, data=data.mis))
summary(lm(Y~X1+X2, data=data.meanImp))
# changes coefficients 

####
# Regression imputation
#
data.RegImp <- data.mis

fit <- lm(X1~Y+X2, data = data.mis)
data.RegImp[is.na(data.mis$X1)&!is.na(data.mis$X2),2] <- predict(fit, data.mis[is.na(data.mis$X1)&!is.na(data.mis$X2),])

fit <- lm(X2~Y+X1, data = data.mis)
data.RegImp[is.na(data.mis$X2)&!is.na(data.mis$X1),3] <- predict(fit, data.mis[is.na(data.mis$X2)&!is.na(data.mis$X1),])

fit = lm(X2~Y, data=data.mis)
data.RegImp[is.na(data.mis$X2)&is.na(data.mis$X1),3] <- predict(fit, data.mis[is.na(data.mis$X2)&is.na(data.mis$X1),])

fit = lm(X1~Y, data=data.mis)
data.RegImp[is.na(data.mis$X2)&is.na(data.mis$X1),2] <- predict(fit, data.mis[is.na(data.mis$X2)&is.na(data.mis$X1),])


data.RegImp

###
colMeans(data.RegImp)
colMeans(data.mis, na.rm=T) 

apply(data.mis, 2, function(x) var(x,na.rm=T))
apply(data.RegImp, 2, function(x) var(x))
# some shrinkage but not much


cov(data.mis[,2:3], use='complete.obs')
cov(data.RegImp[,2:3])
# largely maintained the covariance estimates

summary(lm(Y~X1+X2, data=data.mis))
summary(lm(Y~X1+X2, data=data.RegImp))

##########################################################
data.mis
M <- is.na(data.mis) # missing data pattern.
(M <- as.data.frame(M))


# Begin by filling in missing values by simulating at random from the data. 
# Hot-deck imputation
data.impute <- data.mis
data.impute$X1[M$X1] <- sample(data.mis$X1[!M$X1], size = sum(M$X1), replace=T)
data.impute$X2[M$X2] <- sample(data.mis$X2[!M$X2], size = sum(M$X2), replace=T)

# Check there are no incomplete cases left
sum(ici(data.impute))

# We now cycle around prediction steps
for(i in 1:10){
  # first fit a model to predict X1, using the data only in places where we observed X1.
  fit <- lm(X1~Y+X2, data.impute, subset = !M$X1)
  # Then impute using this regression model plus the added error term
  data.impute$X1[M$X1] <- predict(fit, newdata=data.frame(Y=data.impute$Y[M$X1], X2 = data.impute$X2[M$X1]))+
    rnorm(sum(M$X1), 0, sd=summary(fit)$sigma)

  # Now fit a model to predict X2, using the data only in places where we observed X2.
  fit <- lm(X2~Y+X1, data.impute, subset = !M$X2)
  data.impute$X2[M$X2] <- predict(fit, newdata=data.frame(Y=data.impute$Y[M$X2], 
                                  X1 = data.impute$X1[M$X2]))+ 
                            rnorm(sum(M$X2), 0, sd=summary(fit)$sigma)
}







