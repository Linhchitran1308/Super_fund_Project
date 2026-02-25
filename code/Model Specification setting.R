rm(list=ls())
require(zoo)
require(ggplot2)
require(rugarch)
require(MASS)
require(moments)
source("MRFzR_FunctionsBlock1.R")

source("LoadFundData.R")

par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(P, main="Level", xlab="", ylab="")
plot(r, main="returns", xlab="", ylab="")

#####################################
# VaR model settings                #
#####################################

T  <- length(r)                        # full sample, nobs
N  <- 1250                             # we shorten the sample to the last 5 years
r  <- tail(r,N)                        # log-returns
R  <- coredata(r)
 
p  <- 0.05                             # tolerance level 
H  <- 1                                # horizon

##########################################################
# Historical simulation                                  #
##########################################################

# Value at Risk - VaR
R0       <- sort(R)              
N0       <- floor(N*p)                 # p-th quantile                    
VaR_HS  <- R0[N0]                      # compare to: quantile(R,p)
# Expected shortfall - ES
ES_HS  <- mean(R0[1:N0])               # mean value below VaR

# Results vs function from the MRFzR_FunctionsBlock1.R file
c(VaR_HS,ES_HS)
VaRhist(R,p)

# histogram and VaR
dev.off()
hist(R, main="Historical simulation", xlab="return", ylab="nobs", nclass = 500, col="black")
abline(v=c(VaR_HS,ES_HS) ,lwd=2, col=c("blue","red"),lty=c(2,2))
legend("left", c("VaR for hist. sim.","ES for hist. sim."), lty=c(2,2), col=c("blue","red"), bty="n")

# VaR confidence interval
d1   <- density(R)
x    <- d1$x                              
temp <- tail(which(x<VaR_HS),1)
fhat <- d1$y[temp]

seVaRHS <- sqrt(p*(1-p)/(N*fhat^2))                      # quantile std. deviation
CI      <- c(VaR_HS-1.96*seVaRHS, VaR_HS+1.96*seVaRHS)   # 95% confidence interval 
CI

#######################################################
# Normal distribution                                 #
#######################################################

# Analytical formula
m   <- mean(R) 
s   <- sd(R)

VaR_N     <- qnorm(p)*s + m            # VaR dla stopy zwrotu
ES_N      <- m - s*dnorm(qnorm(p))/p   # Expected shortfall 

# Results vs function from the MRFzR_FunctionsBlock1.R file
c(VaR_N,ES_N)
VaRnorm(R,p)

# density plot
d2  <- density(rnorm(100000, mean=m,sd=s))
plot(d2, main="Normal distribution", xlab="return", ylab="pdf", col="black")
abline(v=c(VaR_N,ES_N) ,lwd=2, col=c("blue","red"),lty=c(2,2))
legend("left", c("VaR for normal dist.","ES for normal dist."), lty=c(2,2), col=c("blue","red"), bty="n")

# And what if we don't know analytical formula?
# Monte Carlo simulations for normal distribution  

M      <- 10000     # number of Monte Carlo simulations
m      <- mean(R) 
s      <- sd(R)

draws0    <- rnorm(M, mean=m, sd=s)
draws1    <- sort(draws0)                 # soreted obs.
M0        <- floor(M*p)                   # p-th quantile                   
VaR_N0    <- draws1[M0]                   # compare to: quantile(draws1,p)
ES_N0     <- mean(draws1[1:M0])           

# Comparison of simulated and analytical VaR/ES
c(VaR_N,VaR_N0)
c(ES_N,ES_N0)

# numerical intergation (quick and accurate)
# VaR_t <- m + s*qdist("std",shape=v1)
qf      <- function(x) qdist("norm", p=x)
VaR_N1  <- m + s*qf(p)
ES_N1  <- m + s*(1/p * integrate(qf, 0, p)$value) 

c(VaR_N,VaR_N1)
c(ES_N,ES_N1)


################################################
# Fat tails - t-Student distribution           #
################################################

# t-Student
# standarization
m  <- mean(R)
s  <- sd(R)
R0 <- (R-m)/s

# estymating degrees of freedom

# method of moments
K     <- kurtosis(R0)
v0    <- 4+6/(K-3)         

# maximum likelihood
dt    <- fitdistr(R0, "t", m = 0, start = list(s=sqrt((v0-2)/v0), df=v0), lower=c(0.001,3))
v1    <- dt$estimate[[2]]  # ML method


VaR_t <- m + s*qt(p,v1) * sqrt((v1-2)/v1)           
qf       <- function(x) qt(x, v1)
ES_t     <- m + s*sqrt((v1-2)/v1)*(1/p * integrate(qf, 0, p)$value)   

# Alternative with rugarch
# require(rugarch)
# temp  <- fitdist("std", R0)
# v1    <- temp$pars["shape"]
# VaR_t <- m + s*qdist("std",shape=v1)
# qf    <- function(x) qdist("std", p=x, shape=v1)
# ES_t  <- m + s*(1/p * integrate(qf, 0, p)$value) 

draws <- rt(100000,v1)*sqrt((v1-2)/v1)*s + m
d3    <- density(draws)
plot(d3, main="t-Student", xlab="return", ylab="pdf", xlim=c(min(R),max(R)),col="black")
abline(v=c(VaR_t,ES_t) ,lwd=2, col=c("blue","red"),lty=c(2,2))
legend("left", c("VaR for t-Stud.","ES for t-Stud."), lty=c(2,2), col=c("blue","red"), bty="n")

# Results vs function from the MRFzR_FunctionsBlock1.R file
c(VaR_t,ES_t)
VaRt(R,p,v1)


# Cornish-Fisher expansion    #
###############################

require(moments)
m   = mean(R)
s   = sd(R)
S   = skewness(R)
K   = kurtosis(R)
p =0.05
PsiP    <- qnorm(p)
VaR_CF  <- m + s*(PsiP + (PsiP^2-1)/6*S + (PsiP^3-3*PsiP)/24*(K-3) - (2*PsiP^3-5*PsiP)/36*(S^2) ) 

# Results vs function from the MRFzR_FunctionsBlock1.R file
VaR_CF
VaRCF(R,p)
# How to compute Expected Shartfall?

1
###########################
# Summary plot and trable #
###########################

plot(d1, main="Comparison of VaR/ES calculations", xlab="retuns", ylab="pdf", xlim=c(min(R),max(R)),col="black", lwd=2)
lines(d2,col="red", lwd=2)
lines(d3,col="green", lwd=2)
abline(v=c(VaR_HS,VaR_N,VaR_t,VaR_CF) ,lwd=2, col=c("black", "red","green","blue"),lty=c(2,2,2,2))
abline(v=c(ES_HS,ES_N,ES_t) ,lwd=2, col=c("black", "red","green"),lty=c(3,3,3))
legend("left", c("HS","Normal","t-Stud.","CF"), lty=c(1,1,1,1), col=c("black", "red","green","blue"), bty="n")



tb              <- matrix(c(VaR_HS,VaR_N,VaR_t,VaR_CF, ES_HS,ES_N,ES_t,NaN),4,2)
colnames(tb)    <- c("VaR","ES")
rownames(tb)    <- c("HS", "Norm", "t-Stud", "CF")
round(-100*tb,2)

# VaR chart for various p  #
############################
tb = matrix(NA,100,5)
for(p in seq(0.001,0.1,0.001)){
  tb[p*1000,] <- c(p,VaRhist(R,p)$VaR, VaRnorm(R,p)$VaR, VaRt(R,p,v1)$VaR, VaRCF(R,p)$VaR)
}

plot(tb[,1],tb[,2], main="Tolerance level and VaR", 
     xlab="tolerance level (p)", ylab="VaR", 
     ylim=c(min(tb[,2:5]),max(tb[,2:5])), xlim=c(0,0.1),col="black", lwd=2, type="l")
lines(tb[,1],tb[,3],col="red", lwd=2, type="l")
lines(tb[,1],tb[,4],col="green", lwd=2, type="l")
lines(tb[,1],tb[,5],col="orange", lwd=2, type="l")
legend("right", c("HS", "Norm", "t-Stud", "CF"), lty=c(1,1,1,1), col=c("black", "red","green","orange"), bty="n")
abline(v=c(0.05) ,lwd=2, lty=2)

tb1 <- tb[tb[,1]<=0.05,]
colMeans(tb1)
c(ES_HS,ES_N,ES_t)



qfCF       <- function(x){
  PsiP     <- qnorm(x)
  VaRCF    <- (PsiP + (PsiP^2-1)/6*S + (PsiP^3-3*PsiP)/24*(K-3) - (2*PsiP^3-5*PsiP)/36*(S^2) ) 
  return(VaRCF)
} 
ES_CF     <- m + s*(1/p * integrate(qfCF, 0, p)$value)



