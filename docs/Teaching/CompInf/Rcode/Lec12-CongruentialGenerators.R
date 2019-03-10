
##########################################################
# Random number generation
#############################################################


rnorm(10)
rnorm(10) # Different numbers

# Now setting the seed
set.seed(1)
rnorm(10)

set.seed(1)
rnorm(10) # The same numbers


M<-2^11
a=51
c=1


cong<-function(x0,n){
  tmp<-c()
  tmp[1]<-x0
  for(i in 1:n){
    tmp[i+1]<-a*tmp[i]%%M 
  }
  return(tmp)
}


out<-cong(1,1000)/M
n<-length(out)
data<-cbind(out[-n], out[-1])
plot(data,xlab="U_i", ylab="U_{i+1}", main="M=2^11, a=51, c=1")



