# Superfund Silver Fund (1068.N) — Risk Analysis
This project analyzes the Superfund Silver dataset using time-series and risk modeling techniques
**Data:** Weekly prices 2009–2026 (Stooq)
**Methods:** GARCH(1,1)-t, VaR & ES (6 methods), Backtesting
superfund-silver-risk-analysis/
│
├── README.md
│
├── data/
│   └── 1068_n.csv                        # Weekly NAV prices from Stooq (2009–2026)
│
├── code/
│   ├── Dataexploration.R                 # Data loading and log-return construction
│   ├── Model__Descriptive_Statistic.R    # Descriptive stats, normality tests, t-Student fit
│   ├── Model_Specification_setting.R     # Static VaR/ES: HS, Normal, t-Student, CF
│   ├── Model_Specification.R             # Dynamic VaR/ES: MA, EWMA, IGARCH, GARCH, eGARCH
│   ├── Backtesting.R                     # VaR backtesting: Kupiec, Christoffersen, Ch-P tests
│   └── Exedance_forecast.R               # PRIIP regulation scenarios, Cornish-Fisher VEV/MRM
│
└── report/
    └── report_1068N_final.docx           # Full written report (Sections 1–5 + References)
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
