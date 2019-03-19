
###################################################
#
# First rejection
# Sample from Beta(alpha, beta)
#
###################################################


alpha <-2
beta <- 3
(M <- (alpha-1)^(alpha-1)*(beta-1)^(beta-1)/(alpha+beta-2)^(alpha+beta-2))
#M<-1 # note the bound M doesn't need to be tight - it is just inefficient otherwise

n <- 10^4

Y<- runif(n) # U[0,1]
U <- runif(n)*M # U[0,M]

keep <- U< Y^(alpha-1)*(1-Y)^(beta-1)
X = Y[keep]

# Check this against the density
hist(X, probability = T)
x<-seq(0,1, 0.001)
y<-dbeta(x, alpha, beta)
lines(x,y, col=2)


plot(Y, U, xlab="x", ylab="Density", pch=3,cex=0.4)
lines(x,y*beta(alpha, beta), col=2, lwd=4)

polygon(c(x,rev(x)),c(rep(0,length(x)),rev(y))*beta(alpha, beta),col="red")

points(Y, U, pch=3,cex=0.4)

# note the bound on this box doesn't need to be tight

############################################
##
## Rejection sampling: Beta distribution 2 
##
############################################

n=10^6


alpha<-4
beta<-1.5

#Implement rejection sampling
Z<-runif(n=n, min=0, max=1)
Y<-Z^(1/alpha)

U <- runif(n=n)
keep <- (U <= (1-Y)^(beta-1))
X <- Y[keep]

# Check this against the pdf
hist(X, probability = T)
x<-seq(0,1, 0.001)
y<-dbeta(x, alpha, beta)
lines(x,y, col=2)



# show what is happening with a picture
y<-x^(alpha-1)*(1-x)^(beta-1)

g<-alpha*x^(alpha-1)/alpha

g.f<-function(x){
  alpha*x^(alpha-1)
}

plot(x,g, type="l", xlab="x", ylab="Density")
lines(x,y)
polygon(c(x,rev(x)),c(rep(0,length(x)),rev(y)),col="red")


Z<-runif(n=10^3, min=0, max=1)
Y<-Z^(1/alpha)
U1<-sapply(Y, function(x) runif(n=1, min=0, max=1/alpha*g.f(x)))
points(Y, U1,pch=3,cex=0.4)

