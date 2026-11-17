# Swiss government bond yields: simple time-series lecture example
# Run from this folder: Rscript swiss_bond_yields_lecture.R
# The only add-on package is forecast. Everything else is base R.

library(forecast)

# 1. Load the five monthly Swiss Confederation yields (percent per year).
data <- read.csv("data/snb_swiss_government_bond_yields_monthly_clean.csv")
dates <- as.Date(data$date)
stopifnot(identical(dates, seq(dates[1], tail(dates, 1), by = "month")),
          !anyNA(data))

columns <- paste0("swiss_government_", c("2y", "5y", "10y", "20y", "30y"))
labels <- c("2 years", "5 years", "10 years", "20 years", "30 years")
yields <- ts(data[columns],
             start = c(as.integer(format(dates[1], "%Y")),
                       as.integer(format(dates[1], "%m"))),
             frequency = 12)
y10 <- yields[, "swiss_government_10y"]

# 2. Describe and plot all five maturities.
describe <- function(x) c(mean = mean(x), sd = sd(x),
                          min = min(x), max = max(x), latest = tail(x, 1))
description <- t(apply(yields, 2, describe))
rownames(description) <- labels
cat("Descriptive statistics (percent per year):\n")
print(round(description, 2))
cat("\nCorrelations among yield levels:\n")
print(round(cor(yields), 2))

dir.create("output", showWarnings = FALSE)
colors <- c("#D97706", "#176B91", "#344054", "#8E44AD", "#4B8B3B")
png("output/ts_01_five_yields.png", width = 1200, height = 600, res = 130)
matplot(as.numeric(time(yields)), yields, type = "l", lty = 1, lwd = 2,
        col = colors, main = "Swiss Confederation bond yields",
        xlab = "Year", ylab = "Yield (% per year)")
abline(h = 0, col = "grey75")
legend("topleft", labels, col = colors, lty = 1, lwd = 2, bty = "n")
dev.off()

# From here on, analyze only the 10-year yield.
# Exercise: repeat the analysis for another maturity.

# 3. Correlograms for the level and the monthly change.
# The AR model below still uses the undifferenced yield level.
change10 <- diff(y10)
png("output/ts_02_correlograms.png", width = 1100, height = 850, res = 130)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
Acf(y10, main = "10-year yield: ACF")
Pacf(y10, main = "10-year yield: PACF")
Acf(change10, main = "Monthly change: ACF")
Pacf(change10, main = "Monthly change: PACF")
dev.off()

# 4. Select an AR(p) model using BIC and observations through 2021.
# Each candidate includes a mean, with no differences or seasonal terms.
# Keeping 2022 onward out of model selection makes the forecast test honest.
training <- window(y10, end = c(2021, 12))
orders <- 1:6
candidates <- lapply(orders, function(p) {
  Arima(training, order = c(p, 0, 0), include.mean = TRUE)
})
bic <- sapply(candidates, BIC)
p <- orders[which.min(bic)]
cat("\nBIC by AR order (lower is better):\n")
print(data.frame(AR_order = orders, BIC = round(bic, 2)), row.names = FALSE)
cat("\nSelected model using data through 2021: AR(", p, ")\n", sep = "")
print(candidates[[which.min(bic)]])

# Refit that AR order on all available data.
fit <- Arima(y10, order = c(p, 0, 0), include.mean = TRUE)
cat("\nFull-sample fit:\n")
print(fit)

# 5. Check residuals. A large Ljung-Box p-value is reassuring.
png("output/ts_03_residual_diagnostics.png", width = 1100, height = 800, res = 130)
checkresiduals(fit)
dev.off()
lb <- Box.test(residuals(fit), lag = 12, type = "Ljung-Box", fitdf = p)
cat("Ljung-Box p-value at lag 12:", round(lb$p.value, 3), "\n")
if (lb$p.value < 0.05) {
  cat("Residual autocorrelation remains; interpret forecasts cautiously.\n")
}

# 6. Forecast the next 12 months from the latest observation.
next_year <- forecast(fit, h = 12, level = c(80, 95))
cat("\n12-month-ahead forecast:\n")
print(next_year)
png("output/ts_04_one_year_forecast.png", width = 1100, height = 600, res = 130)
plot(next_year, xlim = c(tail(time(y10), 1) - 5,
                         tail(time(y10), 1) + 1.1),
     main = "Swiss 10-year yield: next 12 months",
     xlab = "Year", ylab = "Yield (% per year)")
dev.off()

# 7. Recursive forecasts, starting with the information at December 2021.
# Refit the same AR order each month. Compare against the historical mean,
# also recalculated each month. Score only forecasts already realized.
horizons <- c(1, 6, 12)
first_origin <- match(as.Date("2021-12-01"), dates)
if (is.na(first_origin)) stop("December 2021 is missing")
errors <- data.frame(origin = character(), target = character(),
                     horizon = integer(), model_error = numeric(),
                     mean_error = numeric())

for (origin in first_origin:(length(y10) - 1)) {
  history <- ts(y10[1:origin], start = start(y10), frequency = 12)
  refit <- try(Arima(history, order = c(p, 0, 0), include.mean = TRUE),
               silent = TRUE)
  # Near the 2022 yield jump, the usual fitting method may fail; retry with ML.
  if (inherits(refit, "try-error")) {
    refit <- Arima(history, order = c(p, 0, 0),
                   include.mean = TRUE, method = "ML")
  }
  model_forecast <- as.numeric(forecast(refit, h = 12)$mean)
  mean_forecast <- as.numeric(meanf(history, h = 12)$mean)

  for (h in horizons) {
    target <- origin + h
    if (target > length(y10)) next
    errors <- rbind(errors, data.frame(
      origin = as.character(dates[origin]),
      target = as.character(dates[target]), horizon = h,
      model_error = y10[target] - model_forecast[h],
      mean_error = y10[target] - mean_forecast[h]))
  }
}

# RMSFE is root mean squared forecast error. Below 1 favors the AR model.
rmsfe <- function(x) sqrt(mean(x^2))
results <- do.call(rbind, lapply(horizons, function(h) {
  e <- errors[errors$horizon == h, ]
  model <- rmsfe(e$model_error)
  benchmark <- rmsfe(e$mean_error)
  data.frame(horizon_months = h, forecasts = nrow(e),
             model_RMSFE = model, mean_RMSFE = benchmark,
             relative_RMSFE = model / benchmark)
}))
cat("\nRecursive out-of-sample forecast performance:\n")
print(round(results, 3), row.names = FALSE)
write.csv(errors, "output/ts_recursive_forecast_errors.csv", row.names = FALSE)
write.csv(results, "output/ts_relative_rmsfe.csv", row.names = FALSE)

# Plot the errors against the month being forecast.
png("output/ts_05_recursive_forecast_errors.png", width = 1100, height = 850, res = 130)
par(mfrow = c(3, 1), mar = c(3.5, 4.5, 2.5, 1))
for (h in horizons) {
  e <- errors[errors$horizon == h, ]
  target_date <- as.Date(e$target)
  plot(target_date, e$model_error, type = "l", col = "#176B91", lwd = 2,
       ylim = range(c(e$model_error, e$mean_error)),
       main = paste(h, "month(s) ahead"), xlab = "", ylab = "Error (pp)")
  lines(target_date, e$mean_error, col = "#D97706", lwd = 1.5)
  abline(h = 0, col = "grey65", lty = 2)
  if (h == 1) legend("topleft", c("AR model", "Historical mean"),
                     col = c("#176B91", "#D97706"), lty = 1,
                     lwd = c(2, 1.5), bty = "n")
}
dev.off()

cat("\nTime-series charts and forecast-error tables saved in output/.\n")
