
# Monte Carlo z-test

n=10
sigma =1 
t.sample <- c() 
for(i in 1:999){
  temp <- rnorm(n=n, mean=0, sd=sigma)
  t.sample[i] <- mean(temp)
}
z.sample <- t.sample*sqrt(n)/sigma
z <- c(z.sample, 2)

(c<-quantile(z, 0.95))

sum(z>2)/1000

######################################
# Diet example

dietA<-c(233,291,312,250,246,197,268,224)
dietB<-c(185,263,246,224,212,188,250,148)

t.test(dietA,dietB)




r.test<-function(x,y,n){
	n1<-length(x)
	n2<-length(y)
	tobs<-t.test(x,y)$statistic
	t.gen<-rep(0,n)
	alldata<-c(x,y)
	for(i in 1:n){
	    sampledata<-sample(alldata)
	    t.gen[i]<-t.test(sampledata[1:n1],sampledata[(n1+1):(n1+n2)])$statistic
	}
	hist(t.gen,prob=T,col=0,ylim=c(0,0.5))
	z<-seq(from=-4,to=4,length=100)
	lines(z,dt(z,n1+n2-2))
	abline(v=tobs, col=2)
	mean(abs(t.gen)>=abs(tobs)) 
	
}

r.test(dietA,dietB,10000)


r.test2<-function(x,y,n){
	n1<-length(x)
	n2<-length(y)
	tobs2<-mean(x)-mean(y)
	t.gen2<-rep(0,n)
	alldata<-c(x,y)
	for(i in 1:n){
	    sampledata<-sample(alldata)
	    t.gen2[i]<-mean(sampledata[1:n1])-mean(sampledata[(n1+1):(n1+n2)])
	}
	mean(abs(t.gen2)>=abs(tobs2))
}

r.test2(dietA, dietB, 10000)

#######################################
#
# Outliers example
#
#######################################

# Paired data
# Suppose we are given paired data

A<-c(0.33,0.27,0.44,0.28,0.45,0.55,0.44,0.76,0.59,0.01)
B<-c(0.28,0.8,3.72,1.16,1,.63,1.14,.33,.26,.63)
data <- data.frame(A=A, B=B)

t.test(A,B, paired = TRUE)


r.paired<-function(data,n){
  t.vals <- c()
  for(ii in 1:n){
    data.sim <- t(apply(data,1,function(x) sample(x,replace=F)))
    t.vals[ii] <- t.test(data.sim[,1], data.sim[,2],paired=T)$statistic
  }
  print(t.vals)
  mean(abs(t.vals)>=abs(t.test(data[,1], data[,2],paired=T)$statistic))
}

r.paired(data,1000)
##############################################################################
#
# Analysis of varioance example
#
#
############################################################################

y<-c(-.1,-1.1,.74,-3.8,.94,-.3,.67,.86,1.19,-.25,.84,.04,.25,.99,.08,.98,.75,.53)
groups<-as.factor(c(rep("a",4),rep("b",5),rep("c",4),rep("d",5)))

lm0<-lm(y~1)
lm1<-lm(y~groups)
F.stat.obs<-anova(lm0,lm1)$F[2]

# Might be worried about normality for F-distribution
library(MASS)
e<-stdres(lm1)
par(mfrow = c(1,2))
hist(e)
qqnorm(e)
abline(0,1)
    
### Plot the data and the test statistic
par(mfrow=c(1,2))
stripchart(y~groups,vertical=T)
plot(F.stat.obs,0,pch=4,col=2,xlim=c(0,5),xlab="F statistic",yaxt="n",ylab="")

# one way ANOVA, assuming fixed patient effects
par(mfrow=c(1,2))

### Now resample once for illustration
# Assume H_0 true, and randomly re-allocate patients to groups
# Under H_0, dependent variable does not change
newgroups<-sample(groups)

# Plot new data
stripchart(y~newgroups,vertical=T)
new.lm1<-lm(y~newgroups)

# Calculate new F-statistic, plot and compare with old
F.stat.sample<-anova(lm0,new.lm1)$F[2]
plot(F.stat.sample,0,xlim=c(0,5),xlab="F statistic",yaxt="n",ylab="")
points(F.stat.obs,0,pch=4,col=2)


# Resample loads of times and start to build up a histogram of data
N<-200
F.stat.sample<-rep(0,N)
par(mfrow=c(1,2))
for(i in 1:N){
    newgroups<-sample(groups)
    stripchart(y~newgroups,vertical=T)
    new.lm1<-lm(y~newgroups)
    F.stat.sample[i]<-anova(lm0,new.lm1)$F[2]
    
    #Compare distribution of sampled F-statistics under randomisation
    # with theoretical F_{3,14} distribution under H_0
    hist(F.stat.sample[1:i],prob=T,xlab="F statistic",xlim=c(0,13),main="")
    z<-seq(from=0,to=12,length=100)
    lines(z,df(z,3,14),col=4)
    points(F.stat.obs,0,pch=4,col=2)   
    points(F.stat.sample[i],0) 
}

mean(F.stat.sample>=F.stat.obs)  

##############################################################################

A<-c(130,119,119,168,130)
B<-c(154,115,169,137,186)



r.test2<-function(x,y,n){
  n1<-length(x)
  n2<-length(y)
  tobs2<-mean(x)-mean(y)
  t.gen2<-rep(0,n)
  alldata<-c(x,y)
  for(i in 1:n){
    sampledata<-sample(alldata)
    t.gen2[i]<-mean(sampledata[1:n1])-mean(sampledata[(n1+1):(n1+n2)])
  }
  mean(abs(t.gen2)>=abs(tobs2))
}


# Make sure have sources rtest2()
# Test hypothesis that difference in means is k for different k - when do we reject?
k<-c(-20,-10,0,10,20,30,40,50,60)
for(i in 1:length(k)){
   print(c(k[i],r.test2(A+k[i],B,1000)))
}

#############################################################################

r.test.plot1<-function(n){
    x<-c(233,291,312,250,246,197,268,224)
    y<-c(185,263,246,224,212,188,250,148)
    
    z<-c(x,y)
    groups<-as.factor(c(rep("dietA",8),rep("dietB",8)))

tobs<-t.test(x,y)$statistic
t.gen<-rep(0,n)
alldata<-c(x,y)
for(i in 1:n){
readline()
    newgroups<-sample(groups)
    t.gen[i]<-t.test(z[newgroups=="dietA"],z[newgroups=="dietB"])$statistic
    par(mfrow=c(1,2))
    stripchart(z~newgroups,vertical=T)
    plot(t.gen[i],0,xlab="Test statistic",xlim=c(-4,4),main="")
    points(tobs,0,pch=4,col=2)
    
}

mean(abs(t.gen)>=abs(tobs)) 
}


####

r.test.plot2<-function(n){
    x<-c(233,291,312,250,246,197,268,224)
    y<-c(185,263,246,224,212,188,250,148)
    
    z<-c(x,y)
    groups<-as.factor(c(rep("dietA",8),rep("dietB",8)))

tobs<-t.test(x,y)$statistic
t.gen<-rep(0,n)
alldata<-c(x,y)
for(i in 1:n){

    newgroups<-sample(groups)
    t.gen[i]<-t.test(z[newgroups=="dietA"],z[newgroups=="dietB"])$statistic
    par(mfrow=c(1,2))
    stripchart(z~newgroups,vertical=T)
    hist(t.gen[1:i],prob=T,xlab="Test statistic",xlim=c(-4,4),main="",ylim=c(0,0.5))
    zx<-seq(from=-4,to=4,length=100)
    lines(zx,dt(zx,14),col=4)
    points(tobs,0,pch=4,col=2)
    points(t.gen[i],0)
    
}

mean(abs(t.gen)>=abs(tobs)) 
}
#######################
# Fisher's one sample randomisation test

x<-c(10.61,9.46,7.02,11.68,9.58,11.96,11.28,7.73,6.42,8.85)
y<-x-10
Tobs <- abs(mean(y))

y*sample(c(-1,1),10,replace=T)

N<-100
y.mean<-rep(0,N)
for(i in 1:100){
    par(mfrow=c(2,1))
    new.y<-y*sample(c(-1,1),10,replace=T)
    plot(new.y,rep(0,10),xlim=c(-4,4),ylim=c(-0.1,.1),ylab="",yaxt="n",col=c(1:10))
    lines(c(0,0),c(-0.1,.1))
    y.mean[i]<-abs(mean(new.y))
    hist(y.mean[1:i],prob=T,main="")
    points(abs(mean(x-10)),0,pch=4,col=2)
    points(y.mean[i],0)
    readline()
}
mean(y.mean > abs(mean(x-10)))


Tsim <- replicate(10^3, abs(mean(y*sample(c(-1,1),10,replace=T))))
mean(Tsim>Tobs)

t.test(x, mu=10)

#######################
r.test.mac<-function(n){
    x<-c(233,291,312,250,246,197,268,224)
    y<-c(185,263,246,224,212,188,250,148)
    
    z<-c(x,y)
    groups<-as.factor(c(rep("dietA",8),rep("dietB",8)))

tobs<-t.test(x,y)$statistic
t.gen<-rep(0,n)
alldata<-c(x,y)
par(mfrow=c(1,2))
for(i in 1:n){

    newgroups<-sample(groups)
    t.gen[i]<-t.test(z[newgroups=="dietA"],z[newgroups=="dietB"])$statistic
    
    stripchart(z~newgroups,vertical=T)
    hist(t.gen[1:i],prob=T,xlab="Test statistic",xlim=c(-4,4),main="",ylim=c(0,0.5))
    zx<-seq(from=-4,to=4,length=100)
    lines(zx,dt(zx,14),col=4)
    points(tobs,0,pch=4,col=2)
     points(t.gen[i],0)
    readline()
}

mean(abs(t.gen)>=abs(tobs)) 
}

