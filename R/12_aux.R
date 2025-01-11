
#' Check numeric scalar
#'
#' Check whether an object is bounded and coercible to a numeric value.
#'
#' @param x Numeric scalar.
#' @param min Numeric scalar. Minimum value of \emph{x}.
#' @param max Numeric scalar. Maximum value of \emph{x}.
#' @param fun Function to apply to \emph{x} before returning.
#' @param msg String fed to \code{\link[base]{stop}} if an error occurs.
#'
#' @return Returns \code{fun(x)}.
#'
#' @noRd
num_check <- function(
  x, min = 0, max = Inf,
  msg = "Please check the numeric parameters.",
  fun = as.numeric) {

  if(!is.numeric(x) || length(x) != 1 || x < min || x > max) {stop(msg)}

  return(fun(x))
}


#' @noRd
int_check <- function(
  x, min = 0L, max = Inf,
  msg = "Please check the integer parameters.") {

  num_check(x, min, max, msg, fun = as.integer)
}


#' Name hyperparameters
#'
#' Function to help name hyperparameters. Accounts for multiple occurences
#' of \emph{psi} by adding sequential numbers.
#'
#' @param x Character vector. Parameter names.
#' @param M Integer scalar. Number of columns in the data.
#'
#' @return Returns a character vector of adjusted parameter names.
#'
#' @noRd
name_pars <- function(x, M) {

  out <- Reduce(c, sapply(x, function(y) {
    if(y == "psi") {paste0(y, 1:M)} else {y}}))

  return(out)
}


#' Fill credible intervals
#'
#' Helper function to fill data, colours or similar things based on credible
#' intervals. These are used in \code{\link{plot.bvar_irf}} and
#' \code{\link{plot.bvar_fcast}}.
#'
#' Note that transparency may get appended to recycled HEX colours. Also note
#' that no, i.e. a length 0 central element is required when drawing polygons.
#'
#' @param x Scalar or vector. The central element.
#' @param y Scalar or vector. Value(s) to surround the central element with.
#' The first value is closest, values may get recycled.
#' @param P Odd integer scalar. Number of total bands.
#'
#' @return Returns a vector or matrix (if \emph{x} is a vector) of \emph{x},
#' surrounded by \emph{y}.
#'
#' @noRd
fill_ci <- function(x, y, P) {

  n_y <- if(P %% 2 == 0) {
    stop("No central position for x found.")
  } else {P %/% 2}

  fill <- rep(y, length.out = n_y)

  if(length(x) > 1) { # Matrix
    n_row <- length(x)
    return(cbind(t(rev(fill))[rep(1, n_row), ], x, t(fill)[rep(1, n_row), ]))
  } else { # Vector
    return(c(rev(fill), x, fill))
  }
}


#' @noRd
fill_ci_na <- function(x, P) {

  # Corner case when quantiles are missing (t_back or conditional forecasts)
  if(P == 2) {return(if(length(x > 1)) {cbind(x, NA)} else {c(x, NA)})}

  fill_ci(x = x, y = NA, P = P)
}


#' @noRd
fill_ci_col <- function(x, y, P) {

  # Apply transparency to HEX colours
  if(length(y) == 1 && is_hex(y, alpha = FALSE)) {
    y <- paste0(y, alpha_hex(P))
  }

  fill_ci(x = x, y = y, P = P)
}


#' Get a transparency HEX code
#'
#' @param P Integer scalar. Number of total bands.
#'
#' @return Returns a character vector of transparency codes.
#'
#' @importFrom grDevices rgb
#'
#' @noRd
alpha_hex <- function(P) {

  n_trans <- P %/% 2
  out <- switch(n_trans, # Handpicked with love
    "FF", c("FF", "80"), c("FF", "A8", "54"),
    c("FF", "BF", "80", "40"), c("FF", "CC", "99", "66", "33"))

  if(is.null(out)) { # Let rgb() sort it out otherwise
    out <- substr(rgb(1, 1, 1, seq(1, 0, length.out = n_trans)), 8, 10)
  }

  return(out)
}


#' Check valid HEX colour
#'
#' @param x Character scalar or vector. String(s) to check.
#' @param alpha Logical scalar. Whether the string may contain alpha values.
#'
#' @return Returns a logical scalar or vector.
#'
#' @noRd
is_hex <- function(x, alpha = FALSE) {

  if(alpha) return(grepl("^#[0-9a-fA-F]{6,8}$", x))

  return(grepl("^#[0-9a-fA-F]{6,6}$", x))
}


#' Get variable positions
#'
#' Helper functions to aid with variable selection, e.g. in
#' \code{\link{plot.bvar_irf}} and \code{\link{plot.bvar_fcast}}.
#'
#' @param vars Numeric or character vector of variables to subset to.
#' @param variables Character vector of all variable names. Required if
#' \emph{vars} is provided as character vector.
#' @param M Integer scalar. Count of all variables.
#'
#' @return Returns a numeric vector with the positions of desired variables.
#'
#' @noRd
pos_vars <- function(vars, variables, M) {

  if(is.null(vars) || length(vars) == 0L) {
    return(1:M) # Full set
  }
  if(is.numeric(vars)) {
    return(vapply(vars, int_check, # By position
      min = 1, max = M, msg = "Variable(s) not found.", integer(1)))
  }
  if(is.character(vars) && !is.null(variables)) {
    out <- do.call(c, lapply(vars, grep, variables)) # By name
    if(length(out) > 0) {return(out)}
  }

  stop("Variable(s) not found.")
}


#' Name dependent / explanatory variables
#'
#' @param variables Character vector of all variable names.
#' @param M Integer scalar. Count of all variables.
#' @param lags Integer scalar. Number of lags applied.
#'
#' @return Returns a character vector of variable names.
#'
#' @noRd
name_deps <- function(variables, M) {

  if(is.null(variables)) {
    variables <- paste0("var", seq(M))
  } else if(length(variables) != M) {
    stop("Vector with variables is incomplete.")
  }

  return(variables)
}


#' @noRd
name_expl <- function(variables, M, lags) {

  if(is.null(variables)) {
    variables <- name_deps(variables, M)
  }
  explanatories <- c("constant", paste0(rep(variables, lags), "-lag",
    rep(seq(lags), each = length(variables))))

  return(explanatories)
}


#' Compute log distribution function of Inverse Gamma
#'
#' @param x Numeric scalar. Draw of the IG-distributed variable.
#' @param shape Numeric scalar.
#' @param scale Numeric scalar.
#'
#' @return Returns the log Inverse Gamma distribution function.
#'
#' @noRd
p_log_ig <- function(x, shape, scale) {

  return(shape * log(scale) - (shape + 1) * log(x) - scale / x - lgamma(shape))
}


#' Check whether a package is installed
#'
#' @param package Character scalar.
#'
#' @noRd
has_package <- function(package) {

  if(!requireNamespace(package, quietly = TRUE)) {
    stop("Package \'", package, "\' required for this method.", call. = FALSE)
  }

  return(NULL)
}


#' Generate quantiles
#'
#' Check a vector of confidence bands and create quantiles from it.
#'
#' @param conf_bands Numeric vector of probabilities (\eqn{(0, 1)}).
#'
#' @return Returns a sorted, symmetric vector of quantiles.
#'
#' @noRd
quantile_check <- function(conf_bands) {

  conf_bands <- sapply(conf_bands, num_check,
    min = 0 + 1e-16, max = 1 - 1e-16, msg = "Confidence bands misspecified.")

  # Allow only returning the median
  if(length(conf_bands) == 1 && conf_bands == 0.5) {return(conf_bands)}

  # Sort and make sure we have no duplicates (thank mr float)
  quants <- sort(c(conf_bands, 0.5, (1 - conf_bands)))
  quants <- quants[!duplicated(round(quants, digits = 12L))]

  return(quants)
}

#' Set dates
#'
#' Changes rownames to a character string containing dates. Common format of
#' dates or other index is needed if the length of the \emph{instrument} and
#' reduced form residuals differ. In such case, the identification is based on
#' a common subset of residuals and the instrument.
#'
#' @param data Numeric vector, matrix or a data.frame. Note that observations
#' are expected to be ordered from earliest to latest, and variables in the
#' columns.
#' @param start Starting date of the data in form "YYYY-MM-DD".
#' @param frequency Frequency of the data: "year", "quarter", "month", "week" or "day".
#'
#' @return Returns data with rownames containing sequence of dates as character.
set_dates <- function(data, start, frequency) {

  # Check if the input is a data frame, matrix, or numeric vector
  if (!is.data.frame(data) && !is.matrix(data) && !(is.vector(data) && is.numeric(data))) {
    stop("Error: data must be either a data frame, a matrix, or a numeric vector.")
  }

  # Determine the number of elements to set row names for
  num_elements <- if (is.data.frame(data)) {
    nrow(data)
  } else if (is.matrix(data)) {
    nrow(data)
  } else {
    length(data)
  }

  # Set the row names or names based on the sequence of dates
  rownames(data) <- seq(as.Date(start), by = frequency, length.out = num_elements)

  return(data)
}


#' Check IV
#'
#' Intersects the \emph{residuals} and the \emph{instrument} based on an index
#' (e.g. dates) provided in rownames and names, respectively. Returns an error
#'  if the intersection is empty.
#'
#' @param residuals Data frame or numeric matrix. During estimation, draw specific residuals
#' are used and rownames are inherited from data.
#' @param instrument Numeric vector.
#' @param proxyvar character string. Variable for which instrument is provided -
#' Needs to correspond to one of the column names of the data used in
#' \emph{bvar}.
#' @param manual_matching If set to TRUE, user is not expected to specify common
#' index for the data and instrument. Instead the length of instrument needs to
#' match the length of residuals (= length(data inputted in bvar) - lags).
#' @return Returns a list of class \code{check_iv} with information on
#' values and indices corresponding to the intersection. Additionally the list
#' contains inputs provided to the function.
#' @noRd
check_iv <- function(residuals, instrument, manual_matching = FALSE,
                     proxyvar) {

residuals <- data.frame(residuals)

if(manual_matching == TRUE){
  dates_residuals <- NULL
  dates_instrument <- NULL
  common_dates <- NULL

  if(!nrow(residuals) == length(instrument)){
  stop("manual_matching is TRUE but the length of residuals does not match the length of instrument!
  Note that residuals are shorter than the data inputted to bvar, since first l observations of data are not part of residuals,
  where l corresponds to number of lags used.")
  }

  residuals_shortened <- residuals
  instrument_shortened <- instrument
} else{

  if (!proxyvar %in% colnames(residuals)) {
    stop(paste("The proxied variable", proxyvar, "was not found in the column names of residuals."))
  }
  # Get the dates or other indices
  dates_residuals <- rownames(residuals)
  dates_instrument <- names(instrument)

  # Find the intersection of dates (common dates)
  common_dates <- intersect(dates_residuals, dates_instrument)

  # If there are no matches and throw an error
  if (length(common_dates) == 0) {
    stop("Could not match instrument with the data; please provide common indices for both. See `bv_irf` for details.")
    }

  # Shorten both objects to the common dates
  residuals_shortened <- residuals[common_dates,]
  instrument_shortened <- instrument[common_dates]
}

  out <- structure(mget(c("residuals", "instrument",
                          "dates_residuals", "dates_instrument",
                          "residuals_shortened", "instrument_shortened",
                          "common_dates", "manual_matching")),
                       class = "check_iv")

  return(out)
}

# display_matrix <- function(mat) {
#   n <- nrow(mat)
#
#   if (n > 6) {
#     # Display the first 3 rows
#     print(head(mat, 3))
#
#     # Print dots to indicate omitted rows
#     cat("...\n...\n...\n")
#
#     # Display the last 3 rows
#     print(tail(mat, 3))
#   } else {
#     # If the matrix has 6 or fewer rows, display the entire matrix
#     print(mat)
#   }
# }

get_swapped_index <- function(residuals, proxyvar) {
  # Find the column index of the proxy variable
  j <- which(colnames(residuals) == proxyvar)

  # Create a sequence of column indices
  index <- seq_len(ncol(residuals))

  # Swap the first column index with the proxyvar column index
  index[1] <- index[j]
  index[j] <- 1

  return(index)
}


