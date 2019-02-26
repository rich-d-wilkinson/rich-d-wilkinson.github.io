


###########################################
## ECDF
#########################################

par(mfrow=c(1,1))
X <- c(0,1,1.5,2,4)
plot(ecdf(X), verticals=T)


### Cauchy example

set.seed(2)
X <- rcauchy(100)
for(i in 1:50){
  plot(ecdf(X[1:i]), verticals=T, xlim=c(-15,15))
  lines(seq(-15,15,0.1),pcauchy(seq(-15,15,0.1)), col=2)
  readline()
}


plot(ecdf(X), xlim=c(-15,15))
lines(seq(-15,15,0.01), pcauchy(seq(-15,15,0.01)), col=2, lwd=3)

#############################################
# Heart attack example
#############################################

boot<-function(n){
  r1<-rbinom(n,11037,0.00942) #aspirin
  r2<-rbinom(n,11034,0.0171) #placebo
  thetahat<-(r1/11037)/(r2/11034)
  return(thetahat)
}

# run boot function
thetasamples <- boot(10000)

hist(thetasamples,prob=T,col=0) 
quantile(thetasamples,c(0.025,0.975))

