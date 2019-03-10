
#################################################################################
################### MC Example
#################################################################################

h<-function(x){
  (cos(50*x)+sin(20*x))^2
}


## Plot the function
x<-seq(0,1, 0.001)
y<-h(x)
plot(x,y, type="l", xlab="Function")


## Do the integral using Monte Carlo
N<-100000
X<-runif(n=N)
h_realisations <- h(X)
hist(h_realisations, xlab="Generated values of function")

I<-c()
se<-c()
i<-1
for(n in seq(1000,N, 1000)){
  
  I[i]<-mean(h_realisations[1:n])
  se[i]<-1/sqrt(n)*sd(h_realisations)
  i<-i+1
}

plot(seq(1000,N,1000), I, type="l", ylim=c(0.9, 1.10), xlab="Number of samples",
     main="Mean and standard errors")
lines(seq(1000,N,1000), I-2*se, type="l", lty=2)
lines(seq(1000,N,1000), I+2*se, type="l", lty=2)

########################################################################
########################################################################


f <- function(x) {if(x>-1&&x< 1) return(exp(-x^2))
  else return(0)}


# first lets use simpsons rule to estimate the integral
simpson<-function(n){
  x<-seq(from=-1,to=1,length=n)
  d<-x[2]-x[1]
  w<-c(1,4,rep(c(2,4),(n-3)/2),4,1)
  sum(w*sapply(x,f))*d/3 
}

simpson(100)

# now lets use g() = U[-1,1]

x<-runif(n,-1,1)
h<-2*exp(-x^2)
mean(h)

x.plt <- seq(-2,2,len=500)
plot(x.plt, lapply(x.plt, f), type='l', ylim=c(0,1.4), xlim=c(-2,2))
lines(x.plt, dunif(x.plt,-1,1), col=2)

#### 
# now lets use g = N(0,1/2)

x<-rnorm(n,0,0.5^0.5)
h<-pi^0.5*((x>=-1)&(x<=1))
mean(h) 
lines(x.plt, dnorm(x.plt, 0,sd=sqrt(1/2)), col=3)


# no try g(x) = U[0,1]
x<-runif(n,0,1)
h<-exp(-x^2)
mean(h)
lines(x.plt, dunif(x.plt, 0,1), col=4)

# finally try g(x) = N(0,0.09)
x<-rnorm(n,0,0.09^0.5)
h<-(pi*0.18)^0.5*((x>=-1)&(x<=1))*exp(4.56*x^2)
mean(h) 

lines(x.plt, dnorm(x.plt,0,0.09^0.5), col=5)


comparevariances<-function(n){
  x1<-runif(n,-1,1)
  h1<-2*exp(-x1^2)
  v1<-var(h1)
  
  x2<-rnorm(n,0,0.5^0.5)
  h2<-pi^0.5*((x2>=-1)&(x2<=1))
  v2<-var(h2)
  
  x3<-rnorm(n,0,0.09^0.5)
  h3<-(pi*0.18)^0.5*((x3>=-1)&(x3<=1))*exp(4.56*x3^2)
  v3<-var(h3)
  
  print(c(v1,v2,v3))
  
  ssize<-c(10,100,1000,10000,100000,1000000)
  par(mfrow=c(1,1))
  plot(ssize,4*(v3/ssize)^0.5,log="x",type="l",ylab="95% CI Width")
  lines(ssize,4*(v2/ssize)^0.5,lty=2)
  lines(ssize,4*(v1/ssize)^0.5,lty=4)
  legend("topright",c("uniform", "N(0,0.5)","N(0,0.09)"),lty=c(4,2,1))
  
  
} 


comparevariances(10000)
