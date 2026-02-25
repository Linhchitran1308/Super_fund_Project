rm(list=ls())
require(zoo)
require(ggplot2)
require(rugarch)
require(forecast)
require(binom)
require(knitr)
require(car)
source("MRFzR_FunctionsBlock2.R")

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
backN    <- 250                        # backtest sample
w_length <- 1000                       # rolling window length



## Historical simulation
varHS  <- rollapply(r, width=w_length,
                   function(w){quantile(w,p)},
                   by=1, align="right")
varHS  <- stats::lag(varHS, -1)            # lags so that VaR(t+1) = quantile(t)
rr     <- r[index(varHS)]           # rr - realized returns
etaHS  <- tail(rr<=varHS, backN)    # VaR violations)
nHS    <- sum(etaHS); nHS           # number of exceedances
piHS   <- mean(etaHS); piHS         # share of exceedances

dev.off()
VaRplot(alpha=p, actual=tail(rr,backN), VaR=tail(varHS,backN))
title(main="Historical simulation: VaR violations")

# substitute for selected method (from Topic7a.R)
var = varHS
eta = etaHS
pi  = piHS
n1  = nHS

###############################
#    1. Regression test       #
###############################
c(pi,p)
Acf(as.numeric(coredata(eta)))

# 1a. Unconditional coverage. For eta[t] = pi + eps[t] H0: pi = p
reg1a <- lm(eta~1)
summary(reg1a)
linearHypothesis(reg1a, hypothesis.matrix=1, rhs=p, test="Chisq") # H0: pi = p

# 1b. Independence. For eta[t] = pi + rho*eta[t-1] + eps[t] H0: rho = 0  
reg1b <- lm(eta[-1]~stats::lag(eta,-1))
summary(reg1b)
linearHypothesis(reg1b, hypothesis.matrix=c(0,1), rhs=0, test="Chisq")

# 1c. Conditional coverage. For eta[t] = pi + rho*eta[t-1] + eps[t]   H0: pi = p, rho = 0  
reg1c <- lm(eta[-1]~stats::lag(eta,-1))
summary(reg1c)
linearHypothesis(reg1c, hypothesis.matrix=diag(1,2), rhs=c(p,0), test="Chisq")

######################################################
# 2. Unconditional coverage, Binomial distribution   #
#    H0: n1 ~ Binom(n,p), 95% confidence interval    #
######################################################

# 95% confidence interval
kable(cbind(seq(0,20,2),pbinom(seq(0,20,2),backN,p)),col.names = c("n1","cdf"),digits=3)

LB = qbinom(0.025,backN,p)
UB = qbinom(0.975,backN,p)
print(paste0("LB=",LB,"; UB=",UB,", n. violations =",n1))

# Kupiec test
n  <- length(eta); n
n1 <- sum(eta);    n1
n0 <- n - n1;      n0
pi <- n1/n;        pi

# test stat.
LR_uc   = (p/pi)^n1 *((1-p)/(1-pi))^n0
stat_uc = -2*log(LR_uc)
prob_uc = 1 - pchisq(q=stat_uc,df=1, lower.tail=TRUE)
paste0("Kupiec test stat.: ", round(stat_uc,3),"; p-value:", round(prob_uc,3))

# with rugarch package
temp <- VaRTest(alpha=p, actual=rr, VaR=var)
paste0("Kupiec test stat.:", round(temp$uc.LRstat,3),"; p-value:", round(temp$uc.LRp,3))


###########################################################
# 3. Christofersen independence test                      #
#    H0: p(eta[t]=1|eta[t-1]=1) = p(eta[t]=1|eta[t-1]=0)  #
###########################################################
Acf(as.numeric(coredata(tail(eta,backN))))  
## Ljung-Box test, H0: independence
Box.test(as.numeric(coredata(tail(eta, backN))), 1, type="Ljung-Box")


# Christofersen independence test
eta1 <- coredata(eta[-length(eta)])
eta0 <- coredata(eta[-1])
n00 = sum(!eta1 & !eta0)  # no exceedance after no exceedance
n01 = sum(!eta1 &  eta0)  # exceedance after no exceedance
n10 = sum( eta1 & !eta0)  # no exceedance after exceedance
n11 = sum( eta1 &  eta0)  # exceedance after exceedance

#n0  = n00 + n10
#n1  = n01 + n11

pi0 = n01 / (n00+n01) # prob. of exceedance after no exceedance
pi1 = n11 / (n10+n11) # prob. of exceedance after exceedance 
pi  = (n01+n11) / (n00+n01+n10+n11)

c(pi0,pi1,pi)

#LR_ind   = pi^n1*(1-pi)^n0 / (pi0^n01 * pi1^n11 * (1-pi0)^n00  * (1-pi1)^n10)
LR_ind   = (pi/pi0)^n01 * (pi/pi1)^n11 * ((1-pi)/(1-pi0))^n00 * ((1-pi)/(1-pi1))^n10
stat_ind = -2*log(LR_ind)
prob_ind = 1 - pchisq(stat_ind,1)

paste0("Ch1 test stat.: ", round(stat_ind,3),"; p-value:", round(prob_ind,3))


##############################################################
# 4. Christofersen conditional coverage test                 #
#    H0: p(eta[t]=1|eta[t-1]=1) = p(eta[t]=1|eta[t-1]=0) = p #
##############################################################

# Christofersen conditional coverage test
LR_cc   = p^n1*(1-p)^n0 / (pi0^n01 * pi1^n11 * (1-pi0)^n00  * (1-pi1)^n10)
stat_cc = -2*log(LR_cc)
prob_cc = 1 - pchisq(stat_cc,1)
prob_cc

paste0("Ch2 test stat.: ", round(stat_cc,3),"; p-value:", round(prob_cc,3))

# stat_cc = stat_uc + stat_ind
c(stat_uc+stat_ind,stat_cc)

# with rugarch package
temp <- VaRTest(alpha=p, actual=coredata(rr), VaR=coredata(var))
paste0("Ch2 test stat.: ", round(temp$cc.LRstat,3),"; p-value:", round(temp$cc.LRp,3))


#####################################################################
# 5. Tests based on distance between VaR violations                 #
#####################################################################

# calculating time between VaR violations
eta = as.numeric(rr<var)
Te  = length(eta)
tau = which(eta==1)                            # t_i, moments of VaR violation
d   = diff(tau)                                # time between VaR violations
Td  = length(d)

densityPlot(d)
mean(d)
sd(d)

# Simple regression -- just for illustration
# Regression d[i] = a + b*d[i-1] + e[i]
# For a well specified model
# H0: a = 1/p and b=0 

regD <- lm(d[-1]~d[-Td])
summary(regD)

linearHypothesis(regD, hypothesis.matrix=diag(1,2), rhs=c(1/p,0), test="Chisq")

#####################################################################
# Christoffersen and Pelletier (2004)                               #
# Weibull distribution for violations interval                      #
# b=1 constant hazard rate                                          #
# b>1 increasing hazard rate                                        # 
# b<1 decreasing hazard rate                                        #
# H0: b=1 --> exponential distribution                              #
# https://erm.ncsu.edu/wp-content/uploads/sites/41/migrated-files/Christoffersen-Pelletier-Backtesting.pdf
#####################################################################

# if(eta[1]  == 0){d = c(tau[1], d)}            # adjusting for first observation
# if(eta[Te] == 0){d = c(d, Te-tail(tau,1))}    # adjusting for last  observation
# for interested in details see: rugarch:::VaRDurTest and rugarch:::.likDurationW

# Hazard rate in Weibul distribution (probability of event as a function of time from the last event)
# h(t)=a/b(t/b)^(a-1) (if a=1 then h(t)=1/b)
# a - shape parameter
# b - scale parameter (expected duration between exceedances b*Gamma(1+1/a) )
# t - time 

weibull_hazard <- function(t, shape, scale) {
  (shape/scale) * (t/scale)^(shape - 1)
}

t = seq(1,100,1)
a = 1        # shape
b = 1/p      # scale

plot(t,weibull_hazard(t,a,b), type = "l", xlab = "time from exceeddance", ylab = "hazard rate" )

# the same in ggplot
tab = data.frame(time=t, W1=weibull_hazard(t,0.8,1/p), W2=weibull_hazard(t,1,1/p), W3=weibull_hazard(t,1.1,1/p))
require(reshape2)
tb <- melt(tab, id.var = 'time')
colnames(tb)[2]="Alpha"

ggplot(tb,  aes(x = time, y = value, colour = Alpha)) +
  geom_line() +
  theme_light()+
  scale_color_discrete(labels = c("0.8", "1.0", "1.1"))+
  labs(title="Weibull hazard rate", y="hazard rate", x="time from exceedance", caption="")

# Fit of Weibull distribution
library(MASS)
fit <- fitdistr(d, "weibull", method = "BFGS")
fit

# hazazad rate for our data
plot(t,weibull_hazard(t,fit$estimate[1],fit$estimate[2]), type = "l", xlab = "time from exceeddance", ylab = "hazard rate" )


# Log-Likelihood test for a=1
# Comparison of logLik from Weibull and Exponential

fitExp <- fitdistr(d, "exponential")
fitExp

fit$loglik
fitExp$loglik

LR = 2 * (fit$loglik - fitExp$loglik); LR
LRp = 1 - pchisq(LR, 1)              ; LRp  

paste0("Ch-P test stat.: ", round(LR,3),"; p-value:", round(LRp,3))


# this kind of test is ready in rugarch package
# There as some differences (discrete Weibull + accounting for first and last observation)
# To see differences type
# rugarch:::VaRDurTest 
# rugarch:::.likDurationW

temp <- VaRDurTest(alpha=p, actual=coredata(rr), VaR=coredata(var), conf.level = 0.95)
temp
paste0("Ch-P test stat.: ", round(2*(temp$uLL-temp$rLL),3),"; p-value:", round(temp$LRp,3))




##############################
# Parametric models          #
##############################

MAmean  <- rollapply(r, width=w_length, mean, by=1, align="right")
MAstd   <- rollapply(r, width=w_length, sd,   by=1, align="right")
#var     <- MAmean + MAstd*qdist("norm", p);    
var     <- MAmean + MAstd*qdist("std", p, shape=5); 
var     <- stats::lag(var, -1)

dev.off()
VaRplot(alpha=p, actual=rr, VaR=var)
title(main="Parametric model: VaR violations")

temp <- VaRTest(alpha=p, actual=coredata(rr), VaR=coredata(var))
temp1 <- VaRDurTest(alpha=p, actual=coredata(rr), VaR=coredata(var), conf.level = 0.95)
paste0("Kupiec test stat.: ", temp$uc.LRstat,"; p-value:", temp$uc.LRp)
paste0("Ch2 test stat.: ", temp$cc.LRstat,"; p-value:", temp$cc.LRp)
paste0("Ch-P test stat.: ", round(2*(temp1$uLL-temp1$rLL),3),"; p-value:", round(temp1$LRp,3))


#################
# EWMA model    #
#################

lambda <- 0.94
EWMAspec <- ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
                        variance.model=list(model="iGARCH", garchOrder=c(1,1)),
                        fixed.pars=list(alpha1=1-lambda, omega=0), distribution.model = "norm")
                        #fixed.pars=list(alpha1=1-lambda, omega=0,shape=5), distribution.model = "std")

## rolling VaR
varEWMA <- rollapply(r, width=w_length,
                        function(w) {
                            # here we don't estimate but calibrate
                            # forecast on the basis of specification
                            frc <- ugarchforecast(EWMAspec,data=w, n.ahead=1)
                            quantile(frc, p)
                        },
                        by=1, align="right")
varEWMA <- stats::lag(varEWMA, -1)
VaRplot(alpha=p, actual=rr, VaR=varEWMA)
title(main="EWMA: VaR violations")

temp <- VaRTest(alpha=p, actual=coredata(rr), VaR=coredata(varEWMA))
temp1 <- VaRDurTest(alpha=p, actual=coredata(rr), VaR=coredata(var), conf.level = 0.95)
paste0("Kupiec test stat.: ", temp$uc.LRstat,"; p-value:", temp$uc.LRp)
paste0("Ch2 test stat.: ", temp$cc.LRstat,"; p-value:", temp$cc.LRp)
paste0("Ch-P test stat.: ", round(2*(temp1$uLL-temp1$rLL),3),"; p-value:", round(temp1$LRp,3))


##################
# GARCH model    #
##################

GARCHspec <- ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=TRUE),
                         variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                         distribution.model = "norm")
                         # fixed.pars=list(shape=5), distribution.model = "std")

varGARCH <- rollapply(r, width=w_length,
                        function(w) {
                            fit <- ugarchfit(data=w, spec=GARCHspec, solver="hybrid")
                            frc <- ugarchforecast(fit, n.ahead=1)
                            quantile(frc, p)
                        },
                        by=1, align="right")
varGARCH <- stats::lag(varGARCH, -1)

VaRplot(alpha=p, actual=rr, VaR=varGARCH)
title(main="GARCH: VaR violation")

temp <- VaRTest(alpha=p, actual=coredata(rr), VaR=coredata(varGARCH))
temp1 <- VaRDurTest(alpha=p, actual=coredata(rr), VaR=coredata(var), conf.level = 0.95)
paste0("Kupiec test stat: ", temp$uc.LRstat,"; p-value:", temp$uc.LRp)
paste0("Ch2 test stat: ", temp$cc.LRstat,"; p-value:", temp$cc.LRp)
paste0("Ch-P test stat.: ", round(2*(temp1$uLL-temp1$rLL),3),"; p-value:", round(temp1$LRp,3))

#############################################
# Backtesing with ugarchroll from   rugarch #
#############################################

step = 10; # how often we reestimate the GARCH model

varGARCHRoll <- ugarchroll(spec=GARCHspec, data=r, refit.every=step, forecast.length=backN, refit.window="moving", window.size=w_length, calculate.VaR=TRUE, VaR.alpha=c(0.01, 0.025, 0.05))
report(varGARCHRoll, VaR.alpha = 0.05)
report(varGARCHRoll, VaR.alpha = 0.01)


#########################
# Appendix: tests power #
#########################

backN    <- 250                   
alpha    <- 0.05                  # significance level (donot confuse with tolerance level)

###############################
# Traffic lights method       #
###############################

nexceed <- list(0:4,5:9,10:backN)
pval    <- 0.01*c(1,2,3,5)

tab1     <- matrix(0, nrow=4, ncol=3)
colnames(tab1) <- c("n1<=4", "4<n<=9", "9<n"); rownames(tab1) <- paste0("p=",pval)
for(j in 1:4) {
  for(i in 1:3) {
    tab1[j,i] <- round(sum(dbinom(nexceed[[i]], backN, pval[j])), 3)
  }
}
tab1

##########################
# Binomial distribution  #
##########################

# No-rejection region for H0 (norejInt function)
norejInt(p=0.05, n=backN, alpha=alpha)
pbinom(18,backN,0.05) - pbinom(3,backN,0.05)

pval <- 0.01*c(1,2,3,5)
nval <- c(250, 500, 1000, 1500)

tab2 <- matrix(NA,nrow=length(pval), ncol=2*length(nval))
colnames(tab2)  <- paste0("n=",c(rbind(nval, nval)))
rownames(tab2)  <- paste0("p=",pval)

for(i in 1:length(nval)){
  for(j in 1:length(pval)){
    tab2[j,c(2*i-1, 2*i)] <- norejInt(pval[j], nval[i], alpha)
  }}
tab2

# Test power
p      = 0.01   # tolerance level for H0
pTrue  = 0.02   # prob. of VaR violation from DGP

temp <- norejInt(p=p, n=backN, alpha=alpha)  # no-rejection region 
temp[1]:temp[2]
moc = 1-sum(dbinom(temp[1]:temp[2], backN, pTrue))
paste0("The power for p=",p,", when the true prob. of exceedance is :",pTrue," amounts to:", round(moc,3))

#####################
# Kupiec test power #
#####################

p      = 0.01   # tolerance level for H0
pTrue  = 0.02   # prob. of VaR violation from DGP

LR_uc <- function(n, n1, p) {
  n0 <- n - n1
  pi <- n1/n
  (p/pi)^n1 *((1-p)/(1-pi))^n0
}

# no-rejection region
nexceed = 0:backN
temp <- nexceed[(-2*log(LR_uc(backN,nexceed,p)) <= qchisq(1-alpha,1))]
temp
moc = 1-sum(dbinom(temp, backN, pTrue))
paste0("The power for p=",p,", when the true prob. of exceedance is :",pTrue," amounts to:", round(moc,3))

################################################
### Backtesting for ES - McNeil and Frey test  #
################################################

## Historical simulation
##########################################
varHS <- rollapply(r, width=w_length,
                   function(w){quantile(w,p)},
                   by=1, align="right")

esHS  <- rollapply(r, width=w_length,
                   function(w){mean(w[w < quantile(w,p)])},
                   by=1, align="right")

varHS <- stats::lag(varHS, -1)           # lags so that VaR(t+1) = quantile(t)
esHS  <- stats::lag(esHS, -1)

rr     <- r[index(varHS)]          # rr - realized returns
etaHS  <- tail(rr<=varHS, backN)   # VaR violations

# Kupiec and Cristoffersen test for VaR
VaRTest(alpha=p, actual=coredata(rr), VaR=coredata(varHS))

# McNeil and Frey tests
mean(varHS[etaHS])
mean(esHS[etaHS])
mean(rr[etaHS])
ESVaRplot(alpha=p, actual=rr, VaR=varHS, ES=esHS); title(main="Historical simulation")

temp      <- rr[etaHS]-esHS[etaHS]
stat_MF   <- mean(temp)/ (sd(temp)/sqrt(length(temp)))
prob_MF   <- 1- pnorm(abs(stat_MF))
prob_MF 

# The same in rugarch
ESTest(alpha=p, actual=coredata(rr), ES=coredata(esHS), VaR=coredata(varHS))

##########################
## Parametric models     #
##########################

MAmean  <- rollapply(r, width=w_length, mean, by=1, align="right")
MAstd   <- rollapply(r, width=w_length, sd,   by=1, align="right")

# Normal distribution
# VaR
varN <- MAmean + MAstd*qdist("norm", p);          
# VaR violations
etaN <- tail(rr<=varN, backN)
# ES
qf   <- function(x) qdist("norm", p=x)
esN  <- MAmean + MAstd*(1/p * integrate(qf, 0, p)$value)

varN <- stats::lag(varN, -1)           
esN  <- stats::lag(esN, -1)
etaN <- tail(rr<=varN, backN)  

ESVaRplot(alpha=p, actual=rr, VaR=varN, ES=esN); title(main="Normal distribution")
ESTest(alpha=p, actual=coredata(rr), ES=coredata(esN), VaR=coredata(varN))

# t-Student 
# VaR
varT <- MAmean + MAstd*qdist("std", p, shape=5)         
# VaR violations
etaT <- tail(rr<=varT, backN)
# ES
qf       <- function(x) qdist("std", p=x, shape=5)
esT  <- MAmean + MAstd*(1/p * integrate(qf, 0, p)$value)

varT <- stats::lag(varT, -1)           
esT  <- stats::lag(esT, -1)
etaT <- tail(rr<=varT, backN)  

ESVaRplot(alpha=p, actual=rr, VaR=varT, ES=esT); title(main="t-Student distribution")
ESTest(alpha=p, actual=coredata(rr), ES=coredata(esT), VaR=coredata(varT))

#################
## EWMA model   #
#################

lambda <- 0.94

EWMAspec <- ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
                       variance.model=list(model="iGARCH", garchOrder=c(1,1)),
                       fixed.pars=list(alpha1=1-lambda, omega=0), distribution.model = "norm")
#fixed.pars=list(alpha1=1-lambda, omega=0,shape=5), distribution.model = "std")

## quantile function
qf     <- function(x) qdist("norm", p=x)
#qf     <- function(x) qdist("std", p=x, shape=5)
eqt <-  integrate(qf, 0, p)$value

varesEWMA <- rollapply(r, width=w_length,
                       function(w) {
                         frc <- ugarchforecast(EWMAspec,data=w, n.ahead=1)
                         var <- quantile(frc, p)
                         sigma <- sigma(frc)
                         es <- sigma*eqt/p
                         return(c(var, es))
                       },
                       by=1, align="right")
varesEWMA <- stats::lag(varesEWMA, -1)
varEWMA <- varesEWMA[,1]
esEWMA  <- varesEWMA[,2]
etaEWMA <- (rr<varEWMA)

ESVaRplot(alpha=p, actual=rr, VaR=varEWMA, ES=esEWMA); title(main="Model EWMA")
ESTest(alpha=p, actual=coredata(rr), ES=coredata(esEWMA), VaR=coredata(varEWMA))

#######################
## GARCH  model #
#######################

GARCHspec <- ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=TRUE),
                        variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                        distribution.model = "norm")
#fixed.pars=list(shape=5), distribution.model = "std")

qf     <- function(x) qdist("norm", p=x)
#qf     <- function(x) qdist("std", p=x, shape=5)
eqt <-  integrate(qf, 0, p)$value

varesGARCH <- rollapply(r, width=w_length,
                        function(w) {
                          fit <- ugarchfit(data=w, spec=GARCHspec, solver="hybrid")
                          frc <- ugarchforecast(fit, n.ahead=1)
                          var <- quantile(frc, p)
                          sigma <- sigma(frc)
                          mu <- fitted(frc)
                          es <- mu + sigma*eqt/p
                          return(c(var, es))
                        },
                        by=1, align="right")
varesGARCH <- stats::lag(varesGARCH, -1)
varGARCH <- varesGARCH[,1]
esGARCH  <- varesGARCH[,2]
etaGARCH <- (rr<varGARCH)

ESVaRplot(alpha=p, actual=rr, VaR=varGARCH, ES=esGARCH); title(main="GARCH model")
ESTest(alpha=p, actual=coredata(rr), ES=coredata(esGARCH), VaR=coredata(varGARCH))


#############################################################
#        Backtesting entire distribution                    #
#                   Berkowitz test                          #
#############################################################

###########################
# Parametric model       #
###########################

MAmean  <- rollapply(r, width=w_length, mean, by=1, align="right"); MAmean=stats::lag(MAmean,-1)
MAstd   <- rollapply(r, width=w_length, sd,   by=1, align="right"); MAstd=stats::lag(MAstd,-1)

## rstar must be N(0,1))
rstar <- (rr-MAmean)/ MAstd
# PIT~U(0,1)
PIT   <- pdist("norm",rstar)
# PIT   <- pdist("std",rstar,shape=5) # t-Student
rstar <- qnorm(PIT)

regB <- lm(rstar[-1]~stats::lag(rstar,-1))
summary(regB)
linearHypothesis(regB, hypothesis.matrix=diag(1,2), rhs=c(0,0), test="Chisq")

# the same in rugarch (test LR instead of LM)
testDist <- BerkowitzTest(data = rstar, lags = 1, significance = 0.05)
testDist
testDist$LRp 
testDist$H0  
testDist$Decision  

##########
## GARCH #
##########

GARCHspec <- ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=TRUE),
                        variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                        fixed.pars=list(shape=5), distribution.model = "norm")
# fixed.pars=list(shape=5), distribution.model = "std")

temp <- rollapply(r, width=w_length,
                  function(w) {
                    fit <- ugarchfit(data=w, spec=GARCHspec, solver="hybrid")
                    frc <- ugarchforecast(fit, n.ahead=1)
                    sigma <- sigma(frc)
                    mu <- fitted(frc)
                    return(c(mu, sigma))
                  },
                  by=1, align="right")
temp <- stats::lag(temp, -1)
MAmean = temp[,1]; MAstd = temp[,2];
rstar <- (rr-MAmean)/ MAstd
# PIT~U(0,1)
PIT   <- pdist("norm",rstar)
# PIT   <- pdist("std",rstar,shape=5) # t-Student
rstar <- qnorm(PIT)

BerkowitzTest(data = rstar, lags = 1, significance = 0.05)


