#
# Let's start by making some toy data
#
set.seed(1)

X1 <- rnorm(30,0,3)
X2 <- 0.6*X1 + rnorm(30,0,3)
Y = 0.3*X1 + 2*X2 + rnorm(30,0,1)



(fit <- lm(Y~X1+X2))
summary(fit)

miss.rate <- 0.2
X1[which(as.logical(rbinom(30,1,miss.rate)))] <- NA

data.mis <- data.frame(Y=Y, X1=X1, X2=X2)
data.mis
(fit.completecase <- lm(Y~X1+X2, data.mis))

Impute.data <- function(data.mis){
  X1.imp = data.mis$X1
  X1.imp[is.na(X1)] <- sample(X1[!is.na(X1)], size = sum(is.na(X1)), replace=TRUE)
  # sample with replacement from missing values
  data.impute <- data.mis
  data.impute$X1 = X1.imp
  return(data.impute)
} 



m <- 3000
fit.imp <- list()
theta <- matrix(nr=m, nc=3)
V <- matrix(nr=m, nc=3) # we'll just store the diagonal terms in the variance-covariance matrices.
for(i in 1:m){
  data.imp <- Impute.data(data.mis)
  fit <- lm(Y~X1+X2, data = data.imp)
  theta[i,] <- coef(fit)
  V[i,] <- diag(vcov(fit))
}

head(theta)


(theta.hat <- colMeans(theta))
(var.theta <- colMeans(V) + (1+1/m)* diag(cov(theta)))
sqrt(var.theta)
(proportion_of_var_due_to_missingness <- (1+1/m)* diag(cov(theta))/var.theta)

####################################
# lets do the same thing using mice
####################################
library(mice)

m = 1000
data.mice <- mice(data.mis, m=m, method='sample', seed=1)
fit.mice <- with(data.mice, lm(Y~X1+X2))

(fit.mice.pool <- pool(fit.mice))
rowMeans(sapply(fit.mice$analyses,coef))
(theta.hat <- colMeans(theta))


round(summary(fit.mice.pool),2)
sqrt(var.theta)
var.theta <- rowMeans(sapply(lapply(fit.mice$analyses, vcov), diag)) + 
  (1+1/m)* diag(cov(t(sapply(fit.mice$analyses,coef))))#
sqrt(var.theta)

# check what the pool function does.

# compare to the se column in the mice summary

# compare this to the lambda column in the summary
(1+1/m)* diag(cov(t(sapply(fit.mice$analyses,coef))))/var.theta


