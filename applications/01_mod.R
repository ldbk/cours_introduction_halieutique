#devtools::install_github("https://github.com/haddonm/MQMF/")
library(MQMF)


#test modèle de Scheaffer
B1=10+0.2*10*(1-10/100)-1

#fonction
B2<-function(B1,C1,r,K){
  B2<-B1+r*B1*(1-B1/K)-C1
  return(B2)
}

B2(10,1,0.2,100)

B2(10,1,2.8,100)

B2(10,0,seq(0,1,.1),100)

#Code to produce Figure 3.1. Note the two one-line functions
surprod <- function(Nt,r,K) return((r*Nt)*(1-(Nt/K)))
densdep <- function(Nt,K) return((1-(Nt/K)))
r <- 1.2; K <- 1000.0; Nt <- seq(10,1000,10)
par(mfrow=c(2,1),mai=c(0.4,0.4,0.05,0.05),oma=c(0.0,0,0.0,0.0))
par(cex=0.75, mgp=c(1.35,0.35,0), font.axis=7,font=7,font.lab=7)
plot1(Nt,surprod(Nt,r,K),xlab="Population Nt",defpar=FALSE,
      ylab="Production")
plot1(Nt,densdep(Nt,K),xlab="Population Nt",defpar=FALSE,
      ylab="Density-Dependence")


#Code for Figure 3.2. Try varying the value of rv from 0.5-2.8
yrs <- 100; rv=2.8;  Kv <- 1000.0; Nz=100; catch=0.0; p=1.0
ans <- discretelogistic(r=rv,K=Kv,N0=Nz,Ct=catch,Yrs=yrs,p=p)
avcatch <- mean(ans[(yrs-50):yrs,"nt"],na.rm=TRUE) #used in text
label <- paste0("r=",rv," K=",Kv," Ct=",catch, " N0=",Nz," p=",p=p)
plot(ans, main=label, cex=0.9, font=7) #Schaefer dynamics  #


#run discrete logistic dynamics for 600 years
yrs=600
ans <- discretelogistic(r=2.55,K=1000.0,N0=100,Ct=0.0,Yrs=yrs)
plot(ans)
tail(ans,30)
plot(ans[500:600,])

#run discretelogistic and search for repeated values of Nt
yrs <- 600
ans <- discretelogistic(r=2.55,K=1000.0,N0=100,Ct=0.0,Yrs=yrs)
avt <- round(apply(ans[(yrs-100):(yrs-1),2:3],1,mean),2)
count <- table(avt)
count[count > 1] # with r=2.55 you should find an 8-cycle limit

#searches for unique solutions given an r value  see Table 3.2
testseq <- seq(1.9,2.59,0.01)
nseq <- length(testseq)
result <- matrix(0,nrow=nseq,ncol=2,
                 dimnames=list(testseq,c("r","Unique")))
yrs <- 600
for (i in 1:nseq) {  # i = 31
  rval <- testseq[i]
  ans <- discretelogistic(r=rval,K=1000.0,N0=100,Ct=0.0,Yrs=yrs)
  ans <- ans[-yrs,] # remove last year, see str(ans) for why
  ans[,"nt1"] <- round(ans[,"nt1"],3) #try hashing this out
  result[i,] <- c(rval,length(unique(tail(ans[,"nt1"],100))))
}

#Effect of catches on stability properties of discretelogistic
yrs=50; Kval=1000.0
nocatch <- discretelogistic(r=2.56,K=Kval,N0=500,Ct=0,Yrs=yrs)
catch50 <- discretelogistic(r=2.56,K=Kval,N0=500,Ct=50,Yrs=yrs)
catch200 <- discretelogistic(r=2.56,K=Kval,N0=500,Ct=200,Yrs=yrs)
catch300 <- discretelogistic(r=2.56,K=Kval,N0=500,Ct=300,Yrs=yrs)

#Effect of different catches on n-cyclic behaviour Fig3.4
plottime <- function(x,ylab) {
  yrs <- nrow(x)
  plot1(x[,"year"],x[,"nt"],ylab=ylab,defpar=FALSE)
  avB <- round(mean(x[(yrs-40):yrs,"nt"],na.rm=TRUE),3)
  mtext(avB,side=1,outer=F,line=-1.1,font=7,cex=1.0)
} # end of plottime
#the oma argument is used to adjust the space around the graph
par(mfrow=c(2,2),mai=c(0.25,0.4,0.05,0.05),oma=c(1.0,0,0.25,0))
par(cex=0.75, mgp=c(1.35,0.35,0), font.axis=7,font=7,font.lab=7)
plottime(nocatch,"Catch = 0")
plottime(catch50,"Catch = 50")
plottime(catch200,"Catch = 200")
plottime(catch300,"Catch = 300")
mtext("years",side=1,outer=TRUE,line=-0.2,font=7,cex=1.0)

#ajustement

Schaefer<-function(par,data,verbose=FALSE){
  r<-par[["r"]]
  K<-par[["K"]]
  Binit<-par[["Binit"]]
  q<-par[["q"]]
  year <- data$Year
  C <- data$Catch
  I <- data$Index
  n <- length(year)
  B <- numeric(n)
  B[1] <- Binit
  for(i in 1:(n-1))
  {
    B[i+1] <- max(B[i] + r*B[i]*(1-B[i]/K) - C[i], 1)
  }
  Ifit <- q * B

  res <- log(I) - log(Ifit)
  RSS <- sum(res^2)

  pars <- c(r=r, K=K, Binit=Binit, q=q)
  refpts <- c(HRmsy=0.5*r, Bmsy=0.5*K, MSY=0.25*r*K)

  if(verbose)
    list(B=B, HR=C/B, Ifit=Ifit, res=res, pars=pars, refpts=refpts, RSS=RSS)
  else
    RSS
}




Schaefer <- function(par, data, verbose=FALSE)
{
  r <- exp(par[["logr"]])
  K <- exp(par[["logK"]])
  Binit <- exp(par[["logBinit"]])
  q <- exp(par[["logq"]])
  year <- data$Year
  C <- data$Catch
  I <- data$Index
  n <- length(year)
  B <- numeric(n)
  B[1] <- Binit
  for(i in 1:(n-1))
  {
    B[i+1] <- max(B[i] + r*B[i]*(1-B[i]/K) - C[i], 1)
  }
  Ifit <- q * B

  res <- log(I) - log(Ifit)
  RSS <- sum(res^2)

  pars <- c(r=r, K=K, Binit=Binit, q=q)
  refpts <- c(HRmsy=0.5*r, Bmsy=0.5*K, MSY=0.25*r*K)

  if(verbose)
    list(B=B, HR=C/B, Ifit=Ifit, res=res, pars=pars, refpts=refpts, RSS=RSS)
  else
    RSS
}

################################################################################
## South Atlantic albacore

albacore <- read.table("./applications/albacore.dat", header=TRUE)
init <- c(r=0.5, K=200, Binit=100, q=0.5)

Schaefer(par=init, albacore)
optim(init, Schaefer, data=albacore)
est <- optim(init, Schaefer, data=albacore)$par
fit <- Schaefer(est, albacore, verbose=TRUE)

par(mfrow=c(2,2))

plot(albacore$Year, fit$Ifit, ylim=c(0,90), yaxs="i", type="l", lwd=4,
     col="gray", xlab="Year", ylab="Biomass index",
     main="Albacore: Fit to data")
points(Index~Year, albacore)

plot(albacore$Year, fit$B, type="l", ylim=c(0,300), yaxs="i", lwd=2,
     xlab="Year", ylab="Biomass and catch", main="Albacore: Biomass and catch")
points(Catch~Year, albacore, type="h", lwd=6)

plot(albacore$Year, fit$HR, ylim=c(0,0.35), yaxs="i", type="l",
     lwd=2, xlab="Year", ylab="Harvest rate", main="Albacore: Harvest rate")

fit$pars
fit$refpts





#https://haddonm.github.io/URMQMF/surplus-production-models.html#management-advice
################################################################################
## South Atlantic albacore

albacore <- read.table("albacore.dat", header=TRUE)
init <- c(logr=log(0.5), logK=log(200), logBinit=log(100), logq=log(0.5))

Schaefer(par=init, albacore)
optim(init, Schaefer, data=albacore)
est <- optim(init, Schaefer, data=albacore)$par
fit <- Schaefer(est, albacore, verbose=TRUE)

par(mfrow=c(2,2))

plot(albacore$Year, fit$Ifit, ylim=c(0,90), yaxs="i", type="l", lwd=4,
     col="gray", xlab="Year", ylab="Biomass index",
     main="Albacore: Fit to data")
points(Index~Year, albacore)

plot(albacore$Year, fit$B, type="l", ylim=c(0,300), yaxs="i", lwd=2,
     xlab="Year", ylab="Biomass and catch", main="Albacore: Biomass and catch")
points(Catch~Year, albacore, type="h", lwd=6)

plot(albacore$Year, fit$HR, ylim=c(0,0.35), yaxs="i", type="l",
     lwd=2, xlab="Year", ylab="Harvest rate", main="Albacore: Harvest rate")

fit$pars
fit$refpts

################################################################################
## Georges Bank winter flounder

flounder <- read.table("flounder.dat", header=TRUE)

K.init <- 8 * mean(flounder$Catch)
B.init <- 0.5 * K.init
q.init <- flounder$Index[1] / B.init
init <- c(logr=log(0.5), logK=log(K.init),
          logBinit=log(B.init), logq=log(q.init))

Schaefer(par=init, flounder)
optim(init, Schaefer, data=flounder)
optim(init, Schaefer, data=flounder, method="Nelder-Mead",
      control=list(maxit=1e5, reltol=1e-10))
nlminb(init, Schaefer, data=flounder, control=list(eval.max=1e4, iter.max=1e4))
est <- nlminb(init, Schaefer, data=flounder,
              control=list(eval.max=1e4, iter.max=1e4))$par
fit <- Schaefer(est, flounder, verbose=TRUE)

par(mfrow=c(2,2))

plot(flounder$Year, fit$Ifit, ylim=c(0,8), yaxs="i", lwd=4, col="gray",
     type="l", xlab="Year", ylab="Biomass index", main="Flounder: Fit to data")
points(Index~Year, flounder)

plot(flounder$Year, fit$B, type="l", ylim=c(0,15), yaxs="i", lwd=2,
     xlab="Year", ylab="Biomass and catch", main="Flounder: Biomass and catch")
points(Catch~Year, flounder, type="h", lwd=6)

plot(flounder$Year, fit$HR, ylim=c(0,0.6), yaxs="i", type="l",
     lwd=2, xlab="Year", ylab="Harvest rate", main="Flounder: Harvest rate")

t(fit$pars)
t(fit$refpts)
