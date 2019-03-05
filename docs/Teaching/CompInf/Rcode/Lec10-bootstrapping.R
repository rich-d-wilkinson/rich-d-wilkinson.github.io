
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



####################################
##
## loess bootstrap

library(MASS)
nrow <- dim(mcycle)[1]
times <- seq(0,70)
plot(mcycle)
fit <- loess(accel~times, mcycle, span=1/3)
lines(times, predict(fit, times), lty=1)

plot(mcycle)

for(i in 1:100){
  index.boot <- sample(1:nrow, nrow, replace=TRUE)
  mcycle.boot <- mcycle[index.boot,]
  (fit.boot <- loess(accel~times, mcycle.boot, span=1/3))
  lines(times, predict(fit.boot, times), lty=2)
}

#I only did 100 samples here as otherwise the picture begins to look too busy. 
#To get the 95\% CI, we need do a larger number of runs.

times <- 35
pred35ms<-c()
for(i in 1:1000){
  index.boot <- sample(1:nrow, nrow, replace=TRUE)
  mcycle.boot <- mcycle[index.boot,]
  (fit.boot <- loess(accel~times, mcycle.boot, span=1/3))
  pred35ms[i] <- predict(fit.boot, times)
}

hist(pred35ms)
quantile(pred35ms, c(0.025, 0.975))
abline(v =quantile(pred35ms, c(0.025, 0.975)), col=2, lwd=2 )


##########################################################################

#CV
#Reset the plotting window
par(mfrow = c(1,1))

x<-seq(from=0,to=1,length=10)
z<-seq(from=0,to=1,length=100)


# Create some noisy observations
y<-100*(x+rnorm(10,0,0.1))

#Fit a polynomial of order 9 to the data
lm1 <-lm(y~poly(x,9,raw=T))
y2 <-model.matrix(~poly(z,9,raw=T))%*%lm1$coefficients

# Plot it
plot(x,y,ylim=c(min(y2),max(y2)))
lm1 <-lm(y~poly(x,9,raw=T))
y2 <-model.matrix(~poly(z,9,raw=T))%*%lm1$coefficients

lines(z,y2, lty=3, col=4, lwd=2)

lm2 <-lm(y~poly(x,2,raw=T))
y2 <-model.matrix(~poly(z,2,raw=T))%*%lm2$coefficients
lines(z,y2, lty=2, col=2, lwd=2)

deviance(lm2)
deviance(lm1)

fit0 <- lm(y~x)
deviance(lm(y~x))

# Show the real underlying function
#abline(a = 0, b = 100, col = 1, lty=1, lwd=2)
abline(fit0, col=1, lty=1, lwd=2)
legend("topleft", legend = c("Linear", "quadratic", "9th order"), lty=1:3, col=c(1,2,4))


l<-min(y2)
u<-max(y2)
y.pred <- c()

for(i in 1:10){
  readline()
  y.reduced<-y[-i]
  x.reduced<-x[-i]
  lm1<-lm(y.reduced~poly(x.reduced,8,raw=T))
  y2<-model.matrix(~poly(z,8,raw=T))%*%lm1$coefficients
  y.pred[i]<-model.matrix(~poly(x[i],8,raw=T))%*%lm1$coefficients
  plot(x.reduced,y.reduced,ylim=c(l,u),xlim=c(-0.1,1.1))
  lines(z,y2)
  points(x[i],y[i],pch=4,col=2)
  abline(lm(y.reduced~x.reduced),lty=4)
}



# Plot of predicted vs true
y5.pred <- c()
y1.pred <- c()

for(i in 1:10){
  y.reduced<-y[-i]
  x.reduced<-x[-i]
  lm5<-lm(y.reduced~poly(x.reduced,6,raw=T))
  y5.pred[i]<-model.matrix(~poly(x[i],6,raw=T))%*%lm5$coefficients
  
  lm1<-lm(y.reduced~poly(x.reduced,1,raw=T))
  y1.pred[i]<-model.matrix(~poly(x[i],1,raw=T))%*%lm1$coefficients
  
}

mean((y5.pred-y)^2)
mean((y1.pred-y)^2)


plot(y,y5.pred, col=2, xlab='true', ylab='predicted', pch=3)
points(y,y1.pred, pch=2)
abline(a=0,b=1)
legend("bottomright", c('Linear', '6th order', 'x=y (perfect)'), pch=c(2,3, NA), col=c(1,2,1), lty=c(NA, NA, 1))

################################
# 
# Lets now calculate the RMSE for polynomial of different order
#


RMSE <- c()
for(j in 1:8){
  y.pred <- c()
  for(i in 1:10){
    y.reduced<-y[-i]
    x.reduced<-x[-i]
    lm1<-lm(y.reduced~poly(x.reduced,j,raw=T))
    y.pred[i]<-model.matrix(~poly(x[i],j,raw=T))%*%lm1$coefficients
  }
  
  RMSE[j] <- sqrt(mean((y.pred-y)^2))
}

plot(1:8, RMSE, ylab='Root mean square prediction error', xlab = 'Degree of fitted polynomial')
plot(1:8, log(RMSE), ylab='Root mean square prediction error', xlab = 'Degree of fitted polynomial')


# Which model would you prefer?
