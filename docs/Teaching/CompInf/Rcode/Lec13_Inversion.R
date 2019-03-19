

#### Example 7.2.1 - inverse sampling for exponential distribution (continuous)
lambda<-2
u<-runif(100000)
x<- -lambda*log(1-u)
mean(x)

hist(x,prob=T)
z<-seq(from=qexp(0.001,1/lambda),to=qexp(0.999,1/lambda),length=100)
lines(z,dexp(z,1/lambda))

#### Example:  inverse sampling for binomial distribution  (discrete) - specific example Bin(4, 0.3)

# Suppose we have sampled these uniform values
example.u<-c(0.2,0.3,0.7,0.95,0.999)
# To invert then just use qbinom(u, n, p) --- the value x s.t. P(X <= x) > u and P(X <= x-1) < u
qbinom(example.u,4,0.3)

# Try a much larger sample
u<-runif(100000)
x<-qbinom(u,4,0.3) # Should be a sample from a Bin(4, 0.3)

# Check how many of the sample are 2 and compare with theoretical value
mean(x==2)
dbinom(2,4,0.3)

plot(0:4, table(x)/length(x), xlab='x', ylab="Observed frequency")
lines(0:4, dbinom(0:4,4,0.3), type='h', col=2) 

#### Example 7.2.3 - Normal Random Variables
u1<-runif(100000)
u2<-runif(100000)

x1<-(-2*log(u1))^0.5*cos(2*pi*u2)
x2<-(-2*log(u1))^0.5*sin(2*pi*u2)

# They look pretty like a normal
par(mfrow=c(2,1))
z<-seq(from=qnorm(0.001),to=qnorm(0.999),length=100)
hist(x1,prob=T)
lines(z,dnorm(z),col=2)
hist(x2,prob=T)
lines(z,dnorm(z),col=2)

# And they're independent
par(mfrow=c(1,1))
plot(x1,x2)


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

Y<- runif(n)
U <- runif(n)*M

keep <- U<= Y^(alpha-1)*(1-Y)^(beta-1)
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

n=10^4


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

