X <- matrix(c(1,-1,-1,1,2,2), nr=3, byrow=T)

var_u <- function(u){
var(X%*%u)
}

var_u(c(0,-1))

p=5
a <- matrix(c(1:p), nc=1)
sigma = 2*diag(rep(1,p))+a%*%t(a)
eigen(sigma)


## Q3
lambda <- c(4.22 , 2.38 , 1.88 , 1.11 , 0.91 , 0.82 , 0.58, 0.44 , 0.35 , 0.19 , 0.05 , 0.04 , 0.04)

plot(1:length(lambda), lambda/sum(lambda)*100, ylab="Propotion of variance explained", xlab='component number')

plot(1:length(lambda), 100*cumsum(lambda)/sum(lambda), ylab="Cumulative proportion of variance explained",
     xlab="Number of components")


which(cumsum(lambda)/sum(lambda)>0.9)

which(lambda>mean(lambda))

## Q4

R = matrix(c(1 , 0.5792 , 0.2414 , 0.5792 , 1 , 0.5816 , 0.2414 , 0.5816 , 1), nr=3, byrow=TRUE)
eigen(R)

z1 = rep(0,3)
z2= rep(1,3)

t(eigen(R)$vectors) %*%z1
t(eigen(R)$vectors) %*%z2

