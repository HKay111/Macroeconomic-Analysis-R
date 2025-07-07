library(vars)
library(urca)
library(tidyverse)
library(dplyr)
library(mFilter)
library(tseries)
library(ARDL)
library(sandwich)
library(lmtest)
library(dynlm)
library(ggplot2)
library(patchwork)

data <- filename
data <- data[1:75,]
data$Date <- dmy(data$Date)
data <- data[order(data$Date), ]
hp_filtered <- hpfilter(data$Actual_IIP, freq = 14400)
data$Potential_IIP <- hp_filtered$trend
data$Output_Gap    <- hp_filtered$cycle

test_data <- data %>% select("monthly_exc_rate", "Output_Gap", "Inflation")
test_data <- na.omit(test_data)
str(test_data)

# Monthly Exchnage Rate
adf.test(test_data$monthly_exc_rate)
pp.test(test_data$monthly_exc_rate)
pp.test(diff(test_data$monthly_exc_rate))
kpss.test(test_data$monthly_exc_rate)
kpss.test(diff(test_data$monthly_exc_rate))

# Output Gap
adf.test(test_data$Output_Gap)
pp.test(test_data$Output_Gap)
kpss.test(test_data$Output_Gap)


# Inflation
adf.test(test_data$Inflation)
pp.test(test_data$Inflation)
pp.test(diff(test_data$Inflation))
kpss.test(test_data$Inflation, null = "Trend")
kpss.test(diff(test_data$Inflation), null = "Trend")

# ARDL
models <-auto_ardl(monthly_exc_rate ~ Output_Gap+Inflation, data = test_data, max_order = 5)
models$top_orders

# Best model is ARDL(2, 5, 2)

ardl_252 <- models$best_model
ardl_252$order
summary(ardl_252)

uecm_252 <- uecm(ardl_252)
summary(uecm_252)

recm_252 <- recm(uecm_252, case = 2) # Long Run
summary(recm_252)

recm_252_sr <- recm(uecm_252, case = 3) # Short Run
summary(recm_252_sr)

bounds_f_test(ardl_252, case = 3)

tbounds <- bounds_t_test(uecm_252, case = 3, alpha = 0.01)
tbounds
tbounds$tab


# 1. Transform I(1) variables by differencing
d_exc_rate <- diff(test_data$monthly_exc_rate)
d_Inflation <- diff(test_data$Inflation)

# 2. Align I(0) variables (if any)
aligned_Output_Gap <- window(test_data$Output_Gap, start = start(d_exc_rate), end = end(d_exc_rate))

# 3. Combine into one object for the package
var_data <- cbind(d_exc_rate, d_Inflation, aligned_Output_Gap)
colnames(var_data) <- c("d_exc_rate", "d_Inflation", "Output_Gap")


VARselect(var_data, lag.max = 10, type = "const")


# Estimate the VAR(1) model
# The type="const" includes a constant (intercept) in each equation, which is standard.
var_model <- VAR(var_data, p = 1, type = "const")

# Look at the summary output
summary(var_model)

# 1. Test for serial correlation in the residuals
# H0: no serial correlation. You want a high p-value.
serial.test(var_model, lags.pt = 5, type = "PT.asymptotic") # Using 5 lags for the test is a common choice

# 2. Test for heteroskedasticity in the residuals
# H0: no heteroskedasticity. You want a high p-value.
arch.test(var_model, lags.multi = 5)

# 3. Test for normality of the residuals
# H0: residuals are normally distributed. You want a high p-value.
normality.test(var_model)

# 4. Check for structural breaks in the residuals
plot(stability(var_model))

# Exchange Rate
a <- var_model$varresult$d_exc_rate
a_f <- formula(a)
a_d <- model.frame(a)
lm_d_exc_rate <- dynlm(a_f, data = a_d)

# Inflation
b <- var_model$varresult$d_Inflation
b_f <- formula(b)
b_d <- model.frame(b)
lm_d_Inflation <- dynlm(b_f, data = b_d)

# Output Gap
c <- var_model$varresult$Output_Gap
c_f <- formula(c)
c_d <- model.frame(c)
lm_Output_Gap <- dynlm(c_f, data=c_d)

# --- ROBUST COEFFICIENT TESTS ---
# Now this will work without error!
coeftest(lm_d_exc_rate, vcov. = vcovHC(lm_d_exc_rate, type = "HC1"))
coeftest(lm_d_Inflation, vcov. = vcovHC(lm_d_Inflation, type = "HC1"))
coeftest(lm_Output_Gap, vcov. = vcovHC(lm_Output_Gap, type = "HC1"))

# Test: Does d_Inflation Granger-cause d_exc_rate?
waldtest(lm_d_exc_rate, . ~ . - d_Inflation.l1 - d_Inflation.l2, vcov = vcovHC)

# Test: Does Output_Gap Granger-cause d_Inflation?
waldtest(lm_d_Inflation, . ~ . - Output_Gap.l1 - Output_Gap.l2, vcov = vcovHC)

# Test: Does d_exc_rate Granger-cause Output_Gap?
waldtest(lm_Output_Gap, . ~ . - d_exc_rate.l1 - d_exc_rate.l2, vcov = vcovHC)

# Create individual plots
p_exc_rate <- ggplot(data, aes(x = Date, y = monthly_exc_rate)) +
  geom_line(color = "blue") +
  labs(title = "Monthly Exchange Rate", y = "Rate", x = "") +
  theme_minimal()

p_inflation <- ggplot(data, aes(x = Date, y = Inflation)) +
  geom_line(color = "red") +
  labs(title = "Inflation", y = "Percent", x = "") +
  theme_minimal()

p_output_gap <- ggplot(data, aes(x = Date, y = Output_Gap)) +
  geom_line(color = "darkgreen") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Output Gap", y = "Cycle", x = "Date") +
  theme_minimal()

# Combine plots and save
combined_plot <- p_exc_rate / p_inflation / p_output_gap
ggsave("./plots/time_series_original.png", plot = combined_plot, width = 8, height = 6, bg = "white")


# Create a data frame for plotting
plot_df <- data.frame(
  Date = data$Date,
  `Exchange Rate (Level)` = data$monthly_exc_rate,
  `Exchange Rate (Diff)` = c(NA, diff(data$monthly_exc_rate)),
  `Inflation (Level)` = data$Inflation,
  `Inflation (Diff)` = c(NA, diff(data$Inflation))
) %>%
  pivot_longer(-Date, names_to = "Series", values_to = "Value")

# Plot using facets
ggplot(na.omit(plot_df), aes(x = Date, y = Value)) +
  geom_line(aes(color = Series)) +
  facet_wrap(~ Series, scales = "free_y", ncol = 1) +
  labs(title = "Effect of First-Differencing on I(1) Variables", x = "Date", y = "Value") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("./plots/stationarity_transforms.png", width = 8, height = 6, bg = "white")

# Get residuals from the VAR model
residuals_df <- as.data.frame(residuals(var_model))
residuals_df$Date <- tail(data$Date, nrow(residuals_df)) # Align dates

# Pivot for ggplot
residuals_long <- pivot_longer(residuals_df, -Date, names_to = "Equation", values_to = "Residual")

# Create the plot
ggplot(residuals_long, aes(x = Date, y = Residual)) +
  geom_line(alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ Equation, scales = "free_y", ncol = 1) +
  labs(title = "Residuals from VAR(1) Model by Equation",
       subtitle = "Note the changing volatility, confirming heteroskedasticity.",
       x = "Date", y = "Residual") +
  theme_minimal()

ggsave("./plots/var_residuals.png", width = 8, height = 6, bg = "white")

