
presidents <- read.table(file="PresidentData.txt",nrows=21, sep=',', header=TRUE)

i<-1
WinnerParty <- c()
for(pres in presidents$Winner) {
  if(pres  %in% c('Obama', 'Clinton', 'Carter', 'Johnson', 'Kennedy','Truman','Roosevelt')){
    WinnerParty[i] <- 'D'
  }
  else WinnerParty[i] <- 'R'
    i <- i+1
}
presidents$WinnerParty <- WinnerParty

i<-1
IncumbentParty <- c()
for(incumb in presidents$Incumbent) {
  if(incumb  %in% c('Obama', 'Clinton', 'Carter', 'Johnson', 'Kennedy','Truman','Roosevelt')){
    IncumbentParty[i] <- 'D'
  }
  else IncumbentParty[i] <- 'R'
  i <- i+1
}
presidents$IncumbentParty <- IncumbentParty



attach(presidents)
IncumbentRunning <- c()
for(i in 1:dim(presidents)[1]){
  if(as.character(Incumbent[i]) %in% c(as.character(Winner[i]), as.character(RunnerUp[i]))){
    IncumbentRunning[i] <- 'Y'
  }else IncumbentRunning[i] <- 'N'
}
detach(presidents)
presidents$IncumbentRunning <- IncumbentRunning


theta<- c()
#theta[1] = P(change of party when no incumbent)
#theta[2] = P(incumbent wins)

nIncumbentWin
nIncumbentLose
n

presIncumbent <- presidents[IncumbentRunning=='Y',]
presNoIncumbent <- presidents[IncumbentRunning=='N',]

nIncumbentWin <- sum(presIncumbent$WinnerParty==presIncumbent$IncumbentParty)
nIncumbentLose <- sum(presIncumbent$WinnerParty!=presIncumbent$IncumbentParty)

nChangeParty <- sum(presNoIncumbent$WinnerParty != presNoIncumbent$IncumbentParty)
nSameParty <- sum(presNoIncumbent$WinnerParty == presNoIncumbent$IncumbentParty)

likelihood <-function(theta){
  # theta[1] = P(change party | no incumbent)
  #theta[2] <- P(incumbent wins)
  theta[1]^nChangeParty*(1-theta[1])^nSameParty*theta[2]^nIncumbentWin*(1-theta[2])^nIncumbentLose
  }

Loglike <-function(theta){
  nChangeParty *log(theta[1])+
    nSameParty * log(1-theta[1])+
    nIncumbentWin*log(theta[2])+
    nIncumbentLose*log(1-theta[2])
}

# Do by hand
dLdtheta1 <- nChangeParty/theta[1]-nSameParty/(1-theta[1])
nLoglike_trans <-function(p){
  # let p be logit(theta)
  
  theta = exp(p)/(1+exp(p))
  -Loglike(theta)
  }


optim(runif(2), Loglike)

dx<-dy<-0.01
theta.grid <- expand.grid(seq(0+dx,1-dx,dx), seq(0+dy,1-dy,dy))
ll <- apply(theta.grid,1, Loglike)
contour(seq(0+dx,1-dx,dx), seq(0+dy,1-dy,dy), matrix(ll, nr=1/dx-1), nlevels=100, xlab='P(change party | no incumbent)', ylab='P(incumbent wins)')

