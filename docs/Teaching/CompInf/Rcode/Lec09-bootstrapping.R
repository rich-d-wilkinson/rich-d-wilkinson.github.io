# 2019
#############################################
# Heart attack example
#############################################

(theta.obs <- (104/11037)/(189/11034))


boot<-function(n){
  r1<-rbinom(n,11037,0.00942) #aspirin
  r2<-rbinom(n,11034,0.0171) #placebo
  thetahat<-(r1/11037)/(r2/11034)
  return(thetahat)
}
# run boot function
thetasamples <- boot(10000)

hist(thetasamples,prob=T,col=0) 
(qtheta <- quantile(thetasamples,c(0.025,0.975)))

abline(v=qtheta, col=2, lty=2)
abline(v=theta.obs, col=4, lty=3, lwd=3)

# Or calculate with Gaussian approximation 
# estimate standard error
se = sd(thetasamples)
theta.obs+1.96*se
theta.obs-1.96*se

x<- seq(0,1,0.001)
lines(x, dnorm(x, mean(thetasamples), se), col=3, lty=3, lwd=3)
# Gaussian approximation looks good
abline(v =mean(thetasamples)+1.96*se, col=3, lty=3 , lwd=3)
abline(v =mean(thetasamples)-1.96*se, col=3, lty=3 , lwd=3)



#####################################
# Law Schools non-parametric

x<-c(576,635,558,578,666,580,555,661,651,605,653,575,545,572,594)
y<-c(3.39,3.3,2.81,3.03,3.44,3.07,3,3.43,3.36,3.13,3.12,2.74,2.76,2.88,2.96)    
plot(x,y, xlab='LSAT', ylab='GPA')

bootse<-function(n,x,y){
  nx<-length(x)
  bootcor<-rep(0,n)
  for(i in 1:n){
    index<-sample(c(1:nx),replace=T)
    bootcor[i]<-cor(x[index],y[index])
  }
  return(bootcor)  
}

cor.out <- bootse(10000, x, y)
hist(cor.out)

# 95% CI
(qcor <- quantile(cor.out,c(0.025,0.975)))
abline(v=qcor, col=2, lty=2)

# estimated standard error
(se = sd(cor.out))
# use Normal approximation to calculate 95% CI
cor(x,y) + 1.96*se
cor(x,y) - 1.96*se




##############################################
# Mice survival times


treatment<-c(94,197,16,38,99,141,23)
control<-c(52,104,146,10,50,31,40,27,46)

t.test(treatment,control,alternative="greater",var.equal=T)

boot.test<-function(x,y,n){
  alldata<-c(x,y)
  nx<-length(x)
  ny<-length(y)
  
  t.gen<-rep(0,n)
  for(i in 1:n){
    newx<-sample(alldata,nx,replace=T)
    newy<-sample(alldata,ny,replace=T)
    t.gen[i]<-t.test(newx,newy,alternative="greater",var.equal=T)$statistic
  }
  return(t.gen)  
}

t.gen <- boot.test(treatment,control,1000)

hist(t.gen,col=0,prob=T,ylim=c(0,0.4))
z<-seq(from=min(t.gen),to=max(t.gen),length=100)
lines(z,dt(z,14))
mean( t.gen >= t.test(treatment,control,alternative="greater",var.equal=T)$statistic)



##########################################################################
## LSAT EXAMPLE
## - example of model based resampling


x<-c(576,635,558,578,666,580,555,661,651,605,653,575,545,572,594)
y<-c(3.39,3.3,2.81,3.03,3.44,3.07,3,3.43,3.36,3.13,3.12,2.74,2.76,2.88,2.96)    

modelBasedBoot <- function(n){
  fit <- lm(y~x)
  resids <- resid(fit)  
  coefs <- coef(fit)
  fits.boot <- list()
  
  for(i in 1:n){
    yboot <- coefs[1]+coefs[2]*x+sample(resids, replace=TRUE)
    fits.boot[[i]] <- lm(yboot ~ x)
  }
  return(fits.boot)
}

bootstrap.fits <- modelBasedBoot(100)
bootstrap.coefs <- sapply(modelBasedBoot(10), function(x) coef(x))

fit <- lm(y~x)
apply(bootstrap.coefs,1, sd)
summary(fit)

plot(x,y)
abline(fit, lwd=3)

for(i in 1:50){
  abline(bootstrap.fits[[i]], col=2, lty=2, lwd=0.5)
}


