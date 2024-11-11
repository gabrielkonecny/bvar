#Load functions ----
r_files <- list.files(path = "R", pattern = "\\.R$", full.names = TRUE)

for (file in r_files) {
  tryCatch({
    source(file)
    message("Successfully sourced: ", file)
  }, error = function(e) {
    message("Error in file: ", file)
  })
}

# Example 1 - dates provided in rownames ----

data <- readRDS(data, file = "./data/data_with_ebp.rds")
instrument <- readRDS(file = "./data/instrument_MAR21.rds")

# In case of correctly provided rownames:
data <- fred_transform(data, codes = c(1, 4, 4, 1, 1))
x <- bvar(data, lags = 12, n_draw = 1000L, n_burn = 500L, verbose = T)
irf(x) <- irf.bvar(x, bv_irf(horizon = 24, identification = TRUE, instrument = instrument), n_thin = 1L)
plot(irf(x)) # Here, only instrumented shock should be displayed,
# other shocks are not sensibly identified (= garbage)
plot(irf(x), vars_impulse = 1)

#check_iv
check_iv_results <- check_iv(data, instrument)
check_iv_results

# Example 2 - Rownames via interface ----
rownames(data) <- NULL
names(instrument) <- NULL
x <- bvar(data, lags = 12, n_draw = 1000L, n_burn = 500L, verbose = T,
          start_date = "1979-01-01", frequency = "month")

irf(x) <- irf.bvar(x, bv_irf(horizon = 24L, identification = TRUE,
                             instrument = instrument,
                             start_date = "1991-01-01",
                             frequency = "month"), n_thin = 1L)

plot(irf(x), vars_impulse = "GS1")


# Example 3 - No rownames ----
# All IRF will be garbage or calculation will throw error.
x <- bvar(data, lags = 12, n_draw = 1000L, n_burn = 500L, verbose = T)

irf(x) <- irf.bvar(x,
          bv_irf(horizon = 24L, identification = TRUE, instrument = instrument
                 ), n_thin = 1L)
plot(irf(x))




#plot(summary(irf(x), vars_impulse="GS1"))


# Appendix ----
# Data wrangling from example 1, not nice

# Data wrangling skipped, already saved as .rds
# data <- fred_md[, c("GS1", "INDPRO","CPIAUCSL", "UNRATE")]
#
# # Create a sequence of dates starting from 1959-01-01, with monthly frequency
# start_date <- as.Date("1959-01-01")
# num_rows <- nrow(data)  # Get the number of rows in the dataframe
# date_sequence <- seq(start_date, by = "month", length.out = num_rows)
#
# # Set rownames to be this date sequence
# rownames(data) <- date_sequence
#
# # Convert row names to Date type if they are not already
# data$Date <- as.Date(rownames(data))
#
# # 1. Read the EBP data
# EBP <- read.csv(file.path(".", "data", "ebp_csv.csv"), header = TRUE)
#
# EBP <- EBP[,c("date", "ebp")]
#
# # 2. Select the 1st and 3rd columns (assuming 3rd is date) and prepare the date as rownames
# EBP_date <- as.Date(paste0("01", EBP$date), format = "%d%b%Y")
#
# # Set the 'date' column as rownames and remove the date column
#
# EBP <- EBP[, "ebp"]
#
# names(EBP) <- EBP_date
#
# # 3. Assuming `data` is already loaded
# # Convert the rownames of `data` to Date format if not already done
# data$Date <- as.Date(rownames(data))
#
# # 4. Merge `data` and `EBP` by rownames (dates)
# merged_data <- merge(data, EBP, by = "row.names", all = TRUE)
#
# # 5. Set the merged dates back as rownames and remove the redundant 'Row.names' column
# rownames(merged_data) <- merged_data$Row.names
# merged_data <- merged_data[, -which(names(merged_data) == "Row.names")]
#
# # 6. Subset the data between 1979-01-01 and 2019-12-31
# subset_data <- merged_data[rownames(merged_data) >= as.Date("1979-01-01") &
#                              rownames(merged_data) <= as.Date("2019-01-01"), ]
#
# # 7. Rename the column "y" to "ebp"
# colnames(subset_data)[colnames(subset_data) == "y"] <- "EBP"
#
# # 8. Delete the 'Date' column if it exists
# subset_data$Date <- NULL
#
# head(subset_data)
# tail(subset_data)
#
# data <- subset_data
#
# saveRDS(data, file = "./data/data_with_ebp.rds")

# Instrument
# library(readxl)
#
# #Load instrument
# instrument_data <- read_xlsx(file.path(".", "data", "ff4_instruments_shared.xlsx"), range="A2:B325", col_names = c("date", "MPI"))
# instrument_data$date <- as.Date(paste0(instrument_data$date, "01"), format = "%Y%m%d")
#
# # Create a vector using dates as names
# instrument <- setNames(instrument_data$MPI, as.character(instrument_data$date))
# saveRDS(instrument, file = "./data/instrument_MAR21.rds")


