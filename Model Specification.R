rm(list=ls())
require(zoo)
require(ggplot2)
require(rugarch)
require(MASS)
require(moments)
require(forecast)
library("forecast")
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

############################
# Time dependences         #
############################

dev.off()

# No autocorrelatin of returns
Acf(R)
Pacf(R)

# but volatility clustering 
Acf(R*R)
Pacf(R*R)

# and leverage effect?
Ccf(R*R, R)

########################################################################
# Volatility clustering and VaR / ES                                   #
########################################################################

#  For simplicity: we assume t-Student distribution with v=5 
q    <- qdist("std", p=p, shape=5)
qf   <- function(x) qdist("std", p=x, shape=5)

# Alternative, uncomment for normal distribution
# q    <- qdist("norm",p)
# qf   <- function(x) qdist("norm", p=x)

# 1. Unconditional distribution
#######################################################

mu    <- mean(r)
sigma <- sd(r)

# VaR i ES (as numeric integral)
VaR_const   <- mu + sigma*q
ES_const <- mu + sigma*(1/p * integrate(qf, 0, p)$value)   

# plots
par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(merge( r, mu + 2*sigma, mu - 2*sigma), plot.type="single", col=c(1,2,2),
     main="log-returns vs +-2sd", ylab="" )
plot(merge(r, VaR_const, ES_const), plot.type="single", col=c(1,2,3), main=paste(100*p,"% VaR and ES from unconditional variance model",  sep=""), ylab="", xlab="" )
legend("bottomright", c("VaR", "ES"), lty=1, col=2:3)

# 2. Moving average (MA) 
###############################################################

# window length
w_length = 20             

# rolling mean and sd
MAmean  <- rollapply(r, width=w_length, mean, by=1, align="right")
MAstd   <- rollapply(r, width=w_length, sd,   by=1, align="right")
# shifting so that we have distribution from t into t+1
MAmean <- stats::lag(MAmean, -1)
MAstd  <- stats::lag(MAstd,  -1)
# VaR i ES 
MAvar <- MAmean + MAstd*q 
MAes  <- MAmean + MAstd*(1/p * integrate(qf, 0, p)$value)

# plots
par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(merge( r, MAmean + 2*MAstd, MAmean - 2*MAstd), plot.type="single", col=c(1,2,2),
     main="log-returns vs +-2sd", ylab="" )
plot(merge( r, MAvar, MAes), plot.type="single", col=c(1,2,3), main=paste(100*p,"% VaR and ES from MA model",  sep=""), ylab="" )
legend("bottomright", c("VaR", "ES"), lty=1, col=2:3)

# VaR i ES for T+1
mT      <- mean(coredata(tail(r,w_length)))
sT      <- sd(coredata(tail(r,w_length)))
VaR_MA  <- mT + sT*q
ES_MA   <- mT + sT*(1/p * integrate(qf, 0, p)$value)

# 3. Calibrated EWMA (RiskMetrics) 
###############################################################

lambda      <- 0.94                          # smoothing parameter
EWMAsig2    <- rep(0,N)
EWMAsig2[1] <- var(R)                        # starting point
for (t in 2:N){             
  EWMAsig2[t] = lambda * EWMAsig2[t-1]  + (1-lambda) * R[t-1]^2
}
EWMAstd      <- EWMAsig2^0.5
# VaR i ES 
EWMAvar <- EWMAstd*q 
EWMAes  <- EWMAstd*(1/p * integrate(qf, 0, p)$value)

# plots
par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(merge( r, -2*EWMAstd, 2*EWMAstd), plot.type="single", col=c(1,2,2),
     main="log-returns vs EWMA +-2sd", ylab="" )
plot(merge( r, EWMAvar, EWMAes), plot.type="single", col=c(1,2,3), main=paste(100*p,"% VaR and ES from calibrated EWMA",  sep=""), ylab="" )
legend("bottomright", c("VaR", "ES"), lty=1, col=2:3)

# VaR i ES for T+1
mT  <- 0
s2T <- lambda * EWMAsig2[N]  + (1-lambda) * R[N]^2
sT  <- sqrt(s2T)
VaR_EWMA = mT + sT*q
ES_EWMA  = mT + sT*(1/p * integrate(qf, 0, p)$value)

# Results vs function from the MRFzR_FunctionsBlock1.R file
c(VaR_EWMA,ES_EWMA)
VaREWMA(R,p, v=5,lambda=0.94)

# EWMA in rugarch - introduction to GARCH models #
##################################################

lambda = 0.94

EWMAspec <- ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
                       variance.model=list(model="iGARCH", garchOrder=c(1,1)),
                       fixed.pars=list(alpha1=1-lambda, omega=0, shape=5), distribution.model = "std")
#                      fixed.pars=list(alpha1=1-lambda, omega=0),         distribution.model = "norm")

# fitting to the data
EWMAfit <- ugarchfit(data = r, spec = EWMAspec) 
# our, calibrated params
round(coef(EWMAfit),3)

# comparison
plot(merge(sigma(EWMAfit),EWMAstd), plot.type="single", col=c(1,2))
tail(merge(sigma(EWMAfit),EWMAstd))
plot(EWMAfit, which=1)

# 4. Estimating lambda in EWMA: model IGARCH  
###########################################################

IGARCHspec <- ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
                       variance.model=list(model="iGARCH", garchOrder=c(1,1)),
                       fixed.pars=list(omega=0,shape=5), start.pars = list(alpha1=0.06), distribution.model = "std")
#                      fixed.pars=list(omega=0), start.pars = list(alpha1=0.06), distribution.model = "norm")

IGARCHfit <- ugarchfit(data = r, spec = IGARCHspec, solver="hybrid") # if it doesn't work try different solver, e.g. solver="hybrid"
round(coef(IGARCHfit),3)

# VaR i ES 
IGARCHvar <- sigma(IGARCHfit)*q 
IGARCHes  <- sigma(IGARCHfit)*(1/p * integrate(qf, 0, p)$value)

# plots
par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(IGARCHfit, which=1)
plot(merge( r, IGARCHvar, IGARCHes), plot.type="single", col=c(1,2,3), main=paste(100*p,"% VaR and ES from estimated EWMA (IGARCH)",  sep=""), ylab="" )
legend("bottomright", c("VaR", "ES"), lty=1, col=2:3)

# VaR i ES for T+1
IGARCHfct   <- ugarchforecast(IGARCHfit,n.ahead = 1)
mT  <- as.numeric(fitted(IGARCHfct))
sT  <- as.numeric(sigma(IGARCHfct))

VaR_IGARCH  = mT + sT*q
ES_IGARCH   = mT + sT*(1/p * integrate(qf, 0, p)$value)

# Results vs function from the MRFzR_FunctionsBlock1.R file
c(VaR_IGARCH,ES_IGARCH)
VaRIGARCH(r,p, v=5,lambda=0.94)

# 5. Estymating all params: GARCH(1,1)  
###########################################################

GARCHspec <- ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=TRUE),
                         variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                         fixed.pars=list(shape=5), distribution.model = "std")
#                        distribution.model = "norm")

GARCHfit <- ugarchfit(data = r, spec = GARCHspec, solver="nlminb") # hybrid
round(coef(GARCHfit),5)

# VaR i ES 
GARCHvar   <- fitted(GARCHfit) + sigma(GARCHfit)*q 
GARCHes    <- fitted(GARCHfit) + sigma(GARCHfit)*(1/p * integrate(qf, 0, p)$value)

# plots
par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(GARCHfit, which=1) 
plot(merge( r, GARCHvar, GARCHes), plot.type="single", col=c(1,2,3), main=paste(100*p,"% VaR and ES from GARCH(1,1))",  sep=""), ylab="" )
legend("bottomright", c("VaR", "ES"), lty=1, col=2:3)

# VaR and ES for T+1
GARCHfct   <- ugarchforecast(GARCHfit,n.ahead = 1)
mT  <- as.numeric(fitted(GARCHfct))
sT  <- as.numeric(sigma(GARCHfct))

VaR_GARCH  = mT + sT*q
ES_GARCH   = mT + sT*(1/p * integrate(qf, 0, p)$value)

# Results vs function from the MRFzR_FunctionsBlock1.R file
c(VaR_GARCH,ES_GARCH)
VaRGARCH(r,p, v=5)

# 6. Leverage effect: model eGARCH(1,1)  
###########################################################

EGARCHspec <- ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=TRUE),
                        variance.model=list(model="eGARCH", garchOrder=c(1,1)),
                        fixed.pars=list(shape=5), distribution.model = "std")
#                        distribution.model = "norm")

EGARCHfit <- ugarchfit(data = r, spec = EGARCHspec, solver="nlminb")
round(coef(EGARCHfit),5)

# VaR and ES 
EGARCHvar   <- fitted(EGARCHfit) + sigma(EGARCHfit)*q 
EGARCHes    <- fitted(EGARCHfit) + sigma(EGARCHfit)*(1/p * integrate(qf, 0, p)$value)

# plots
par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(EGARCHfit, which=1)
plot(merge( r, EGARCHvar, EGARCHes), plot.type="single", col=c(1,2,3), main=paste(100*p,"% VaR and ES from EGARCH(1,1))",  sep=""), ylab="" )
legend("bottomright", c("VaR", "ES"), lty=1, col=2:3)

# VaR and ES for T+1
EGARCHfct   <- ugarchforecast(EGARCHfit,n.ahead = 1)
mT  <- as.numeric(fitted(EGARCHfct))
sT  <- as.numeric(sigma(EGARCHfct))

VaR_EGARCH  = mT + sT*q
ES_EGARCH   = mT + sT*(1/p * integrate(qf, 0, p)$value)

# Results vs function from the MRFzR_FunctionsBlock1.R file
c(VaR_EGARCH,ES_EGARCH)
VaRGARCH(r,p, 5, 1,1,TRUE)

#Model sGarch
#6. Leverage effect: model sGARCH(1,1)  
###########################################################

sGARCHspec <- ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=TRUE),
                         variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                         fixed.pars=list(shape=5), distribution.model = "std")
#                        distribution.model = "norm")

sGARCHfit <- ugarchfit(data = r, spec = sGARCHspec, solver="nlminb")
round(coef(sGARCHfit),5)

# VaR and ES 
sGARCHvar   <- fitted(sGARCHfit) + sigma(sGARCHfit)*q 
sGARCHes    <- fitted(sGARCHfit) + sigma(sGARCHfit)*(1/p * integrate(qf, 0, p)$value)

# plots
par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(sGARCHfit, which=1)
plot(merge( r, sGARCHvar, sGARCHes), plot.type="single", col=c(1,2,3), main=paste(100*p,"% VaR and ES from sGARCH(1,1))",  sep=""), ylab="" )
legend("bottomright", c("VaR", "ES"), lty=1, col=2:3)

# VaR and ES for T+1
sGARCHfct   <- ugarchforecast(sGARCHfit,n.ahead = 1)
mT  <- as.numeric(fitted(sGARCHfct))
sT  <- as.numeric(sigma(sGARCHfct))

VaR_sGARCH  = mT + sT*q
ES_sGARCH   = mT + sT*(1/p * integrate(qf, 0, p)$value)

# Results vs function from the MRFzR_FunctionsBlock1.R file
c(VaR_sGARCH,ES_sGARCH)
VaRGARCH(r,p, 5, 1,1,TRUE)

sGARCHspec_norm <- ugarchspec(
  mean.model=list(armaOrder=c(0,0), include.mean=TRUE),
  variance.model=list(model="sGARCH", garchOrder=c(1,1)),
  distribution.model="norm"
)
sGARCHfit_norm <- ugarchfit(data=r, spec=sGARCHspec_norm, solver="nlminb")

n <- length(r)

LL_std  <- likelihood(sGARCHfit)
LL_norm <- likelihood(sGARCHfit_norm)

LLn_std  <- LL_std / n
LLn_norm <- LL_norm / n

IC_std  <- infocriteria(sGARCHfit)       # AIC, BIC...
IC_norm <- infocriteria(sGARCHfit_norm)

c(LLn_std=LLn_std, LLn_norm=LLn_norm)
IC_std
IC_norm


# All together
##################################

dev.off()
labssd <- c("log-returns", "unconditional", "MA", "EWMA", "IGARCH", "GARCH", "eGARCH")
plot( merge(r,
            2*sigma,
            2*MAstd,
            2*EWMAstd,
            2*sigma(IGARCHfit),
            2*sigma(GARCHfit),
            2*sigma(EGARCHfit),
            -2*sigma,
            -2*MAstd,
            -2*EWMAstd,
            -2*sigma(IGARCHfit),
            -2*sigma(GARCHfit),
            -2*sigma(EGARCHfit)),
         facets=NULL,
         main="log-returns vs +-2sd",
     col= c(1:7, 2:7),
     plot.type="single",
     ylab="", xlab="")
legend("bottomright", labssd, lty=1, col=c(1:7, 2:7))


# VaR plots 
###################
labsvar <- c("log-returns", "staly rozklad", "MA", "EWMA", "IGARCH", "GARCH", "eGARCH")
plot(merge(r,
           VaR_const,
           MAvar,
           EWMAvar,
           IGARCHvar,
           GARCHvar,
           EGARCHvar),
     main=paste("log-returns vs ", 100*p,"% VaR", sep=""),
     col=1:7,
     plot.type="single",
     ylab="", xlab="")
legend("bottomright", labsvar, lty=1, col=1:7)

# all together
-100*c(VaR_const,VaR_MA, VaR_EWMA, VaR_IGARCH, VaR_GARCH, VaR_EGARCH)
-100*c(ES_const,ES_MA, ES_EWMA, ES_IGARCH, ES_GARCH, ES_EGARCH)

# using funkctions in MRFzR_FunkcjeBlok1.R
VaR_const1  <- VaRt(R,p,5)$VaR
VaR_MA1     <- VaRt(tail(R,w_length),p,5)$VaR
VaR_EWMA1   <- VaREWMA(r, p, v=5)$VaR  
VaR_IGARCH1 <- VaRIGARCH(r, p, v=5)$VaR
VaR_GARCH1  <- VaRGARCH(r, p, v=5)$VaR
VaR_EGARCH1  <- VaRGARCH(r, p, v=5, eGARCH=TRUE)$VaR


# comparison
-100*cbind(
c(VaR_const,VaR_MA, VaR_EWMA, VaR_IGARCH, VaR_GARCH, VaR_EGARCH),
c(VaR_const1,VaR_MA1, VaR_EWMA1, VaR_IGARCH1, VaR_GARCH1, VaR_EGARCH1))

