rm(list=ls())
require(zoo)
require(ggplot2)
require(rugarch)
#require(rje)
#require(lubridate)
#require(dplyr)
require(bizdays)
require(xts)
#require(XML)
#require(tools)
#require(timeDate)
#require(RQuantLib)

source("LoadFundData.R")

ret = r


par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(P, main="Level", xlab="", ylab="")
plot(r, main="returns", xlab="", ylab="")

#####################################
# VaR model settings                #
#####################################

N  <- 1260                             # we shorten the sample to the last 5 years
r  <- tail(r,N)                        # log-returns
R  <- coredata(r)
p  <- 0.025                             # tolerance level

#######################################################################
#  Computing VaR, VEV i MRM according to Annex II of PRIIP regulation #
#  Cornish-Fisher for p=2.5% and 3-year horizon                       #
#######################################################################

Nyear      <- 252               # nobs in year
RHP        <- 3                 # RHP - recommended holding period
N          <- RHP * Nyear       # RHP in days

# Stats
M0 <- length(R)
M1 <- sum(R)/M0          
R0 <- R - M1             
M2 <- sum(R0^2)/M0
M3 <- sum(R0^3)/M0
M4 <- sum(R0^4)/M0

sigma <- sqrt(M2)       

mu1 <- M3/(sigma^3)     ## skewness
mu2 <- M4/(sigma^4) -3  ## kurtosis

# VaR, VEV i MRM (in RHP)
VaR <-   sigma*sqrt(N)*(-1.96+0.474*mu1/sqrt(N) - 0.0687*mu2/N + 0.146*(mu1^2)/N) - 0.5*(sigma^2)*N
VEV <-   (sqrt(3.842 - 2*VaR) - 1.96)/ sqrt(RHP)

# VEV interpretation: compare to annualized std
sigA = sigma*sqrt(Nyear)

MRM  <- as.numeric(cut(VEV, breaks = c(-Inf, 0.005, 0.05, 0.12, 0.2, 0.3, 0.8, Inf), labels = 1:7, right = FALSE))
MRM

#####################################################################################
# Scenarios according to Annex Iv of the regulation                                 #
#####################################################################################
V <- 1000   #  initial value of a portfolio

# setting RHP in line with the regulation 
if(RHP<=1){ Horyz <- RHP } else if(RHP<3){
    Horyz <- c(1,RHP)
} else {
    Horyz <- c(1,ceiling(RHP/2),RHP)
}
# Horizon i days
HoryzN   <- Horyz*Nyear

# Results
Sc           <- matrix(NA,nrow=length(Horyz),ncol=5)
colnames(Sc) <- c("Years", "Unfavourable", "Moderate", "Favourable", "Stressed")
Sc[,1]       <- Horyz

for(h in 1:length(HoryzN)){
    H0         <- HoryzN[h]
    RHP0       <- Horyz[h]
    q = c(0.1,0.5,0.9)  # quantiles for scenarios 1,2 and 3
    for(ScN in 1:3){    # scenarios 1,2 and 3
        PsiU    <- qnorm(q[ScN])
        temp    <- M1*H0 + sigma*sqrt(H0)*(PsiU + (PsiU^2-1)/6*mu1/sqrt(H0) + (PsiU^3-3*PsiU)/24*mu2/H0 - (2*PsiU^3-5*PsiU)/36*(mu1^2)/H0) - 0.5*(sigma^2)*H0
        Sc[h,ScN+1]  <- exp(temp)*V
    }

    # Stressed scenario
    sigma_sstar <- SigStressed(r,RHP0)
    if(RHP0<=1){PsiU <- qnorm(0.01)} else {PsiU <- qnorm(0.05)}
    temp        <- sigma_sstar*sqrt(H0)*(PsiU + (PsiU^2-1)/6*mu1/sqrt(H0) + (PsiU^3-3*PsiU)/24*mu2/H0 - (2*PsiU^3-5*PsiU)/36*(mu1^2)/H0) - 0.5*(sigma_sstar^2)*H0
    Sc[h,5]     <- exp(temp)*V
}
  # Table
Sc <- data.frame(round(Sc,3))
Sc



#####################################################
#  stressed sigma                                   #
#  T0       - horizon in years                      #
#  ret0     - returns                               #
#####################################################

SigStressed <- function(ret,RHP){

  Ret <- coredata(ret)

  if(RHP<=1){w <- 21}
  if(RHP >1){w <- 63}
  # window length from point 10(c) of Annex IV
  w <- w+1

  # rolling std (point 10c, Annex IV)
  # interwals ends (from w to Tobs)
  seq <- sapply(w:(length(Ret)), function(t){
    retA       <- Ret[(t-w+1):t]                           # w+1 obsevations
    retA       <- retA- mean(retA)                         # de-meaning 
    Mw         <- length(retA)                             
    sigma_w    <- sqrt(sum(retA^2)/Mw)
    return(sigma_w)}
  )
  sigma_sw <- zoo(seq, order.by = index(ret)[w:(length(ret))])

  # values for scenario in line with the formula from point 11 of the annex from regulation 
  if(RHP<=1){
    sigma_sstar <- as.numeric(quantile(sigma_sw, 0.99))
    PsiU        <- qnorm(0.01)
  } else {
    sigma_sstar <- as.numeric(quantile(sigma_sw, 0.90))
    PsiU        <- qnorm(0.05)
  }
  return(sigma_sstar)
}



