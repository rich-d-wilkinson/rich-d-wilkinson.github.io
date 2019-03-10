

#####################################################################################
# Normal example
#####################################################################################

n=20
t.sample <- replicate(999,
                mean(rnorm(n)))

#t.sample <- c() 

#for(i in 1:999){
 # temp <- rnorm(n=n, mean=0, sd=1)
  #t.sample[i] <- mean(temp)
#}
z.sample <- t.sample*sqrt(n)/1
z <- c(z.sample, 2) # recall Zobs = 2


sum(z>2)/1000
1-pnorm(2)

hist(z.sample, probability=T)
curve(dnorm, -4,4,add=TRUE)


#####################################################################################
# Exams example
#####################################################################################

exams<-matrix(c(3,8,4,8,5,3,4,0),nrow=2,ncol=4)
chisq.test(exams)

examsreduced<-matrix(c(3,8,4,8,9,3),nrow=2,ncol=3)
chisq.test(examsreduced)

#################################################


prob<-apply(exams,2,sum)/sum(exams)
expec<-c(prob*16,prob*19)
tobs<-sum( (c(exams[1,],exams[2,])-expec)^2/expec)

tgen  <- replicate(1000,
                   {
                     obs<-c(rmultinom(1,16,prob),rmultinom(1,19,prob))
                     sum((obs-expec)^2/expec)
                   }
)

par(mfrow=c(1,1))
hist(tgen,prob=T,ylim=c(0,0.25),col=0,xlim=c(0,30))
lines(seq(from=0,to=25,length=100),dchisq(seq(from=0,to=25,length=100),3))

sum(tgen>tobs)
sum(tgen>tobs)/length(tgen)




#####################################################################################
### Testing for spatial randomness ####
#####################################################################################


# Read in the observed data
x<-matrix(c(0.64183627, 0.37724091, 0.6238449, 0.27059961, 0.52431424, 0.07702952, 0.40041743, 
    0.60477365, 0.43723551, 0.77048353, 0.3708506, 0.0084252787, 0.24240117, 
    0.73177749, 0.48169126, 0.58705295, 0.70793248, 0.052239959, 0.36798376, 
    0.35522715, 0.28808044, 0.46695876, 0.83134469, 0.49769194, 0.7685349, 0.35695953,
    0.50752601, 0.33841213, 0.63149462, 0.45216901, 0.60864363, 0.6488713, 0.48915889,
    0.69615915, 1., 0.61816525, 0.50339415, 0.70128727, 0.51788491, 0.26905284, 
    0.60979083, 0.50484525, 0.44681789, 0.77901966, 0.90460022, 0.39411187, 0.5280443,
    0.33311544, 0.65149943, 0., 0.90152143, 0.56449707, 0.74448572, 0.68764737, 
    0.56885903, 0.28586782, 0.28512201, 0.34925324, 0.42320853, 0.85028059, 0.409249,
    0., 0.64011021, 0.62635965, 0.59147053, 0.66682285, 0.70156484, 0.37372079, 
    0.46225028, 0.1756221, 0.25491359, 0.63169307, 0.58757665, 0.5359571, 0.60643316,
    0.27393282, 0.40908448, 0.41437412, 0.73234747, 0.55014021, 0.89311939, 0.66129336,
    0.29782663, 0.74479574, 0.96833546, 0.75142247, 0.70796747, 1., 0.51019547, 
    0.49757075, 0.68018413, 0.49468926, 0.9932952, 0.74938044, 0.60818068, 0.88398524,
    0.69053409, 0.26836528, 0.83488132, 0.51063196)
, nrow = 50, ncol = 2)

plot(x)
##################################################
## A visual illustration of our test statistic
 
par(mfrow=c(1,1))
n<-length(x[,1])
mx1<-matrix(x[,1],n,n,byrow=F)
mx2<-matrix(x[,2],n,n,byrow=F)
distances<-((mx1-t(mx1))^2+(mx2-t(mx2))^2)^0.5
distances<-distances+diag(100,n,n)

plot(x)

#indexes = apply(distances,1, which.min)

for(i in 1:50){
  index<-which.min(distances[i,])  #  alternatively: which(distances[i,]==min(distances[i,]))
  lines(c(x[i,1],x[index,1]),c(x[i,2],x[index,2]))        
}

1/mean(apply(distances,2,min))


##################################################
# A function which, given observations x find the test stat 

nnsum<-function(x){
    n<-length(x[,1])
    mx1<-matrix(x[,1],n,n,byrow=F)
    mx2<-matrix(x[,2],n,n,byrow=F)
    distances<-((mx1-t(mx1))^2+(mx2-t(mx2))^2)^0.5
    distances<-distances+diag(100,n,n) # so that we don't pick the distance
    # to ourselves
    return(mean(apply(distances,2,min)))
}

##################################################
# A function which perform the Monte Carlo Test
# Arguments x (observed data) and n number of MC reps
mctest<-function(x,n) {
    	tobs<-1/nnsum(x)
    	test.stat<-vector("numeric",n)
	for(i in 1:n){
    		rx<-matrix(runif(100),nrow=50,ncol=2)
    		test.stat[i] <- 1/    nnsum(rx)
	} 
	sum(test.stat>tobs)
}
##################################################
# Illustration of four randomly sampled datasets under H0
#x11()
par(mfrow=c(2,2))
for(j in 1:4){

	rx<-matrix(runif(100),nrow=50,ncol=2)

	n<-length(rx[,1])
	mx1<-matrix(rx[,1],n,n,byrow=F)
	mx2<-matrix(rx[,2],n,n,byrow=F)
    	distances<-((mx1-t(mx1))^2+(mx2-t(mx2))^2)^0.5
    	distances<-distances+diag(100,n,n)

	plot(rx)


	for(i in 1:50){
		index<-which(distances[i,]==min(distances[i,]))
	lines(c(rx[i,1],rx[index,1]),c(rx[i,2],rx[index,2]))        
	}
}

