# Superfund Silver Fund (1068.N) — Risk Analysis
This project analyzes the Superfund Silver dataset using time-series and risk modeling techniques
**Data:** Weekly prices 2009–2026 (Stooq)
**Methods:** GARCH(1,1)-t, VaR & ES (6 methods), Backtesting
superfund-silver-risk-analysis/

How to run:
# 1. Load and prepare data
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
