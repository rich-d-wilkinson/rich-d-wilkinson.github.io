
##########################################################################
#Cross validation
###########################################################################

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

