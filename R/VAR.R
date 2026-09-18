# Re-analysis of the Macroeconomic Analysis (VAR) project.
#
# The original script had two problems. First, the IIP cycle was shifted one
# month relative to the differenced variables: diff() keeps months 2-75 while
# window() kept months 1-74, so the VAR mixed up the timing. Second, the
# vars::causality() output was reported as three pairwise tests, but it is
# actually a block test of one variable against all the other equations.
# This version aligns the data properly and reports the tests for what they
# are.
#
# The sample stops in March 2020, because the next rows in the data file jump
# to June 2020 and the file is not continuous after that point.

library(readr)
library(lubridate)
library(ggplot2)
library(patchwork)

library(vars)
library(urca)
library(ARDL)
library(mFilter)

library(sandwich)
library(sessioninfo)

set.seed(123)

# --- 1. Data ----------------------------------------------------------------

data <- read_csv("data/filename.csv", show_col_types = FALSE)
data$Date <- dmy(data$Date)
data <- data[order(data$Date), ]
data <- data[data$Date <= "2020-03-01", ]
data$IIP_growth_yoy <- data$Actual_IIP
data$iip_cycle <- as.numeric(hpfilter(data$IIP_growth_yoy, freq = 14400)$cycle)

# The VAR uses changes for the two I(1) variables and the cycle in levels.
# The cycle is shifted back one month so that all three series are dated the
# same way.
var_data <- data.frame(
  d_exc_rate = diff(data$monthly_exc_rate),
  d_inflation = diff(data$Inflation),
  iip_cycle = data$iip_cycle[-1]
)

# --- 2. Plots ---------------------------------------------------------------

p_exc <- ggplot(data, aes(Date, monthly_exc_rate)) +
  geom_line(colour = "#0072B2") +
  labs(title = "Monthly exchange rate", y = "INR per USD", x = NULL) +
  theme_minimal(base_size = 11)

p_inf <- ggplot(data, aes(Date, Inflation)) +
  geom_line(colour = "#D55E00") +
  labs(title = "Inflation", y = "Percent", x = NULL) +
  theme_minimal(base_size = 11)

p_cyc <- ggplot(data, aes(Date, iip_cycle)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_line(colour = "#009E73") +
  labs(title = "IIP growth cycle (HP)", y = "Percentage points", x = "Date") +
  theme_minimal(base_size = 11)

ggsave("plots/time_series_original.png", p_exc / p_inf / p_cyc,
       width = 8, height = 8, dpi = 300)

# Levels and first differences of the two I(1) series
p1 <- ggplot(data, aes(Date, monthly_exc_rate)) +
  geom_line(colour = "#0072B2") +
  labs(title = "Exchange rate (level)", x = NULL, y = NULL) +
  theme_minimal(base_size = 10)

p2 <- ggplot(data.frame(Date = data$Date[-1], x = diff(data$monthly_exc_rate)),
             aes(Date, x)) +
  geom_line(colour = "#0072B2") +
  labs(title = "Exchange rate (change)", x = NULL, y = NULL) +
  theme_minimal(base_size = 10)

p3 <- ggplot(data, aes(Date, Inflation)) +
  geom_line(colour = "#D55E00") +
  labs(title = "Inflation (level)", x = NULL, y = NULL) +
  theme_minimal(base_size = 10)

p4 <- ggplot(data.frame(Date = data$Date[-1], x = diff(data$Inflation)),
             aes(Date, x)) +
  geom_line(colour = "#D55E00") +
  labs(title = "Inflation (change)", x = NULL, y = NULL) +
  theme_minimal(base_size = 10)

ggsave("plots/stationarity_transforms.png", (p1 + p2) / (p3 + p4),
       width = 9, height = 6, dpi = 300)

# --- 3. Unit root tests -----------------------------------------------------
# ADF with up to 6 lags chosen by AIC, and KPSS. Same deterministic terms for
# every series: an intercept, no trend.

summary(ur.df(data$monthly_exc_rate, type = "drift", lags = 6, selectlags = "AIC"))
summary(ur.kpss(data$monthly_exc_rate, type = "mu", lags = "short"))

summary(ur.df(data$Inflation, type = "drift", lags = 6, selectlags = "AIC"))
summary(ur.kpss(data$Inflation, type = "mu", lags = "short"))

summary(ur.df(data$iip_cycle, type = "drift", lags = 6, selectlags = "AIC"))
summary(ur.kpss(data$iip_cycle, type = "mu", lags = "short"))

# --- 4. Cointegration -------------------------------------------------------

ardl_model <- auto_ardl(monthly_exc_rate ~ Inflation + iip_cycle,
                        data = as.data.frame(data), max_order = 5, selection = "AIC")

bounds_f_test(ardl_model$best_model, case = 3, exact = TRUE, R = 20000)
bounds_t_test(uecm(ardl_model$best_model), case = 3, exact = TRUE, R = 20000)

# --- 5. VAR -----------------------------------------------------------------

VARselect(var_data, lag.max = 10, type = "const")
var_model <- VAR(var_data, p = 1, type = "const")
summary(var_model)
roots(var_model)

serial.test(var_model, lags.pt = 5, type = "PT.asymptotic")
arch.test(var_model, lags.multi = 5)
normality.test(var_model)

png("plots/var_residuals.png", width = 1200, height = 800, res = 150)
plot.ts(residuals(var_model), main = "VAR residuals")
dev.off()

# --- 6. Block Granger tests -------------------------------------------------
# Each call tests whether one variable helps predict all the other equations
# together. These are block tests, not pairwise tests.

causality(var_model, cause = "d_inflation")$Granger
causality(var_model, cause = "iip_cycle")$Granger
causality(var_model, cause = "d_exc_rate")$Granger

causality(var_model, cause = "d_inflation", vcov. = vcovHC)$Granger
causality(var_model, cause = "iip_cycle", vcov. = vcovHC)$Granger
causality(var_model, cause = "d_exc_rate", vcov. = vcovHC)$Granger

# --- 7. Robustness: drop March 2020 -----------------------------------------
# March 2020 is the first pandemic month and it drives the non-normal
# residuals in the full sample.

data_short <- data[data$Date <= "2020-02-01", ]
var_data_short <- data.frame(
  d_exc_rate = diff(data_short$monthly_exc_rate),
  d_inflation = diff(data_short$Inflation),
  iip_cycle = as.numeric(hpfilter(data_short$IIP_growth_yoy, freq = 14400)$cycle)[-1]
)

var_short <- VAR(var_data_short, p = 1, type = "const")
normality.test(var_short)
causality(var_short, cause = "d_inflation")$Granger

# --- 8. Package versions ----------------------------------------------------

writeLines(capture.output(session_info()), "session_info.txt")
