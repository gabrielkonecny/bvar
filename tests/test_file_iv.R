# List all .R files in the subdirectory
r_files <- list.files(path = "R", pattern = "\\.R$", full.names = TRUE)

# Filter only the files that start with a number
r_files <- r_files[grepl("^\\d", basename(r_files))]


for (file in r_files) {
  tryCatch({
    source(file)
    message("Successfully sourced: ", file)
  }, error = function(e) {
    message("Error in file: ", file)
    # You can also print the error message if you want
    # message("Error message: ", e$message)
  })
}

#library("BVAR")

# Access a subset of the fred_qd dataset
data <- fred_md[, c("GS1", "INDPRO","CPIAUCSL", "UNRATE")]

# Create a sequence of dates starting from 1959-01-01, with monthly frequency
start_date <- as.Date("1959-01-01")
num_rows <- nrow(data)  # Get the number of rows in the dataframe
date_sequence <- seq(start_date, by = "month", length.out = num_rows)

# Set rownames to be this date sequence
rownames(data) <- date_sequence

# Convert row names to Date type if they are not already
data$Date <- as.Date(rownames(data))

# 1. Read the EBP data
EBP <- read.csv(file.path(".", "data", "ebp_csv.csv"), header = TRUE)

EBP <- EBP[,c("date", "ebp")]

# 2. Select the 1st and 3rd columns (assuming 3rd is date) and prepare the date as rownames
EBP_date <- as.Date(paste0("01", EBP$date), format = "%d%b%Y")

# Set the 'date' column as rownames and remove the date column

EBP <- EBP[, "ebp"]

names(EBP) <- EBP_date

# 3. Assuming `data` is already loaded
# Convert the rownames of `data` to Date format if not already done
data$Date <- as.Date(rownames(data))

# 4. Merge `data` and `EBP` by rownames (dates)
merged_data <- merge(data, EBP, by = "row.names", all = TRUE)

# 5. Set the merged dates back as rownames and remove the redundant 'Row.names' column
rownames(merged_data) <- merged_data$Row.names
merged_data <- merged_data[, -which(names(merged_data) == "Row.names")]

# 6. Subset the data between 1979-01-01 and 2019-12-31
subset_data <- merged_data[rownames(merged_data) >= as.Date("1979-01-01") &
                             rownames(merged_data) <= as.Date("2019-01-01"), ]

# 7. Rename the column "y" to "ebp"
colnames(subset_data)[colnames(subset_data) == "y"] <- "ebp"

# 8. Delete the 'Date' column if it exists
subset_data$Date <- NULL

head(subset_data)
tail(subset_data)

data <- subset_data

head(data)
tail(data)

rownames(data) <- NULL
data


# Transform data
data <- fred_transform(data, codes = c(1, 4, 4, 1, 1))


# Estimate a BVAR
x <- bvar(data, lags = 12, n_draw = 1000L, n_burn = 500L, verbose = T, start_date = "1979-01-01", frequency = "month")
resid(x)[,1:5]

library(readxl)

#alternatively without dependence on dplyr
instrument_data <- read_xlsx(file.path(".", "data", "ff4_instruments_shared.xlsx"), range="A2:B325", col_names = c("date", "MPI"))
instrument_data$date <- as.Date(paste0(instrument_data$date, "01"), format = "%Y%m%d")

# Create a vector using dates as names
instrument <- setNames(instrument_data$MPI, as.character(instrument_data$date))
names(instrument) <- NULL
# For identification, if IV is shorter than residuals, subset residuals.
# y <- intersect_vectors_by_date(resid(x)[,1:5], instrument) #From 62b_proxy_var.R
# y

irf(x) <- irf.bvar(x, bv_irf(horizon = 24L, identification = TRUE,
          instrument = instrument, start_date = "1991-01-01",
          frequency = "month"), n_thin = 1L)

#u <- resid(x)
plot(irf(x))

#plot(summary(irf(x), vars_impulse="GS1"))

#############################################################################


