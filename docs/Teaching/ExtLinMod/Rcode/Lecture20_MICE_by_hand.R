
library(mice)
set.seed(1)

X1 <- rnorm(30,0,3)
X2 <- 0.6*X1 + rnorm(30,0,3)
Y = 0.3*X1 + 2*X2 + rnorm(30,0,1)

(fit <- lm(Y~X1+X2))
summary(fit)

miss.rate <- 0.2
X1[which(as.logical(rbinom(30,1,miss.rate)))] <- NA
X2[which(as.logical(rbinom(30,1,miss.rate)))] <- NA

data.mis <- data.frame(Y=Y, X1=X1, X2=X2)
(fit.completecase <- lm(Y~X1+X2, data.mis))


data.mis
M <- is.na(data.mis) # missing data pattern.
M <- as.data.frame(M)
# Begin by filling in missing values by simulating at random from the data.
data.impute <- data.mis
data.impute$X1[M$X1] <- sample(data.mis$X1[!M$X1], size = sum(M$X1), replace=T)
data.impute$X2[M$X2] <- sample(data.mis$X2[!M$X2], size = sum(M$X2), replace=T)

# Check there are no imcomplete cases left
sum(ici(data.impute))

# We now cycle around prediction steps
for(i in 1:10){
  # first fit a model to predict X1, using the data only in places where we observed X1.
  fit <- lm(X1~Y+X2, data.impute, subset = !M$X1)
  # Then impute using this regression model plus the added error term
  data.impute$X1[M$X1] <- predict(fit, newdata=data.frame(Y=data.impute$Y[M$X1], 
                                                          X2 = data.impute$X2[M$X1]))+
    rnorm(sum(M$X1), 0, sd=summary(fit)$sigma)
  
  # Now fit a model to predict X2, using the data only in places where we observed X2.
  fit <- lm(X2~Y+X1, data.impute, subset = !M$X2)
  data.impute$X2[M$X2] <- predict(fit, newdata=data.frame(Y=data.impute$Y[M$X2], 
                                                          X1 = data.impute$X1[M$X2]))+ 
    rnorm(sum(M$X2), 0, sd=summary(fit)$sigma)
  print(data.impute)
}

##################
#
# We can do the same thing with mice in a single command
#
##################

mice(data.mis, method='norm')
