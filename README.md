# Superfund Silver Fund (1068.N) — Risk Analysis
This project analyzes the Superfund Silver dataset using time-series and risk modeling techniques
**Data:** Weekly prices 2009–2026 (Stooq)
**Methods:** GARCH(1,1)-t, VaR & ES (6 methods), Backtesting
- How to run:
1. Load and prepare data
source("Dataexploration.R")

# 2. Descriptive statistics and distributional analysis
source("Model__Descriptive_Statistic.R")

# 3. Static VaR/ES (HS, Normal, t-Student, Cornish-Fisher)
source("Model_Specification_setting.R")

# 4. Dynamic volatility models (EWMA, IGARCH, GARCH, eGARCH)
source("Model_Specification.R")

# 5. Formal backtesting (Kupiec, Christoffersen, ES tests)
source("Backtesting.R")

# 6. PRIIP regulation scenarios
source("Exedance_forecast.R")

1. Data exploraion:
Handles data loading and preprocessing. Downloads or reads the 1068.N CSV file from Stooq, constructs the price series P as a zoo object, and computes daily log returns r = diff(log(P)). Includes an automated data quality check (minimum 2,000 observations, no gaps > 7 days). This script is sourced by all other scripts via source("LoadFundData.R").

2.Model__Descriptive_Statistic.R

Performs distributional analysis of the log return series over the last 5 years (N = 1,250 daily observations). Computes and annualises the four central moments (mean, variance, skewness, kurtosis). Runs three formal normality tests: D'Agostino (skewness), Anscombe-Glynn (kurtosis), and Jarque-Bera. Fits a Student-t distribution by method of moments and maximum likelihood, and produces the t-Student density comparison plot and QQ plot used in the report.

3. Model_Specification_setting.R
 
Estimates static VaR and ES under four methods: Historical Simulation, Normal, Student-t (v = 5), and Cornish-Fisher expansion. Computes confidence intervals for the HS VaR using kernel density estimation. Produces a summary density comparison plot and a VaR sensitivity chart across tolerance levels p ∈ [0.001, 0.10].

4.Model_Specification.R

Estimates dynamic volatility models and computes time-varying VaR and ES. Implements five volatility specifications in sequence: (1) unconditional constant variance, (2) Moving Average (MA), (3) calibrated EWMA (λ = 0.94), (4) estimated IGARCH, (5) GARCH(1,1) with Student-t innovations (v = 5), and (6) eGARCH for the leverage effect. Compares Normal vs Student-t GARCH using AIC, BIC and log-likelihood per observation. Produces conditional volatility plots and overlaid VaR charts across all models.

5.Backtesting.R
Implements formal VaR and ES backtesting using a rolling window of 1,000 observations over a 250-day backtest sample. For VaR backtesting, runs five tests:

Kupiec (1995) — unconditional coverage (proportion of failures)
Christoffersen (1998) — independence and conditional coverage
Christoffersen & Pelletier (2004) — Weibull-based duration test (hazard rate analysis)
Regression-based tests (unconditional, independence, conditional coverage via OLS)
McNeil & Frey — ES backtest based on standardised exceedance residuals

Also includes the Berkowitz test for full distributional backtesting (PIT-based), and a traffic-lights / binomial power analysis. Covers HS, parametric (Normal, t-Student), EWMA, and GARCH models, with exceedance plots for visual diagnostics.
Key outputs: Kupiec LR statistic, Christoffersen CC statistic, Ch-P duration test, ES test p-values, exceedance plots
