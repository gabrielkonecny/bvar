library("BVAR")

# Access a subset of the fred_qd dataset
data <- fred_qd[, c("CPIAUCSL", "UNRATE", "FEDFUNDS")]

# Transform it to be stationary
data <- fred_transform(data, codes = c(4, 4, 1))
#data <- fred_transform(data, codes = c(4, 4, 1))

# Estimate a BVAR using one lag, default settings and very few draws
x <- bvar(data, lags = 1, n_draw = 600L, n_burn = 100L, verbose = FALSE)

# Compute + store IRF with a longer horizon, no identification and thinning
irf(x) <- irf(x, bv_irf(horizon = 24L, identification = FALSE), n_thin = 5L)
u <- resid(x)

#############################################################################

library(readxl)

#alternatively without dependence on dplyr
instrument_data <- read_xlsx(file.path(".", "data", "ff4_instruments_shared.xlsx"), range="A2:B325", col_names = c("date", "MPI"))
instrument_data$date <- as.Date(paste0(instrument_data$date, "01"), format = "%Y%m%d")

# Create a vector using dates as names
instrument <- setNames(instrument_data$MPI, as.character(instrument_data$date))

# For identification, if IV is shorter than residuals, subset residuals.
y <- intersect_vectors_by_date(resid(x)[,3], instrument) #From 62b_proxy_var.R
y

iv_stats(y$residuals, y$instrument)
