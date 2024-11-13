
#' @export
print.bvar <- function(x, ...) {

  y <- x[["meta"]]

  cat("Bayesian VAR consisting of", y[["N"]], "observations,",
      y[["M"]], "variables and", y[["lags"]], "lags.")
  cat("\nTime spent calculating:", format(round(y[["timer"]], 2)))
  cat("\nHyperparameters:",
      paste(x[["priors"]][["hyper"]], collapse = ", "),
      "\nHyperparameter values after optimisation:",
      paste(round(x[["optim"]][["par"]], 5), collapse = ", "))
  cat("\nIterations (burnt / thinning): ", y[["n_draw"]], " (", y[["n_burn"]],
      " / ", y[["n_thin"]], ")", sep = "")
  cat("\nAccepted draws (rate): ", y[["accepted"]], " (",
      round(y[["accepted"]] / (y[["n_draw"]] - y[["n_burn"]]), 3),
      ")\n", sep = "")

  return(invisible(x))
}

#' @export
print.check_iv <- function(x, ...) {

  # Extract the start and end dates from the data and instrument date ranges
  start_date_data <- x$dates_residuals[1]
  end_date_data <- tail(x$dates_residuals, 1)
  start_date_instrument <- x$dates_instrument[1]
  end_date_instrument <- tail(x$dates_instrument, 1)

  # Get the number of matched dates and display the first and last matched observations
  num_matches <- length(x$common_dates)
  matched_display <- paste0(
    x$common_dates[1], " (value ", x$instrument_shortened[x$common_dates[1]],
    ") until ",
    x$common_dates[num_matches], " (value ", x$instrument_shortened[x$common_dates[num_matches]],
    ")"
  )

  # Print the summary information
  cat("Found residuals from", start_date_data, "until", end_date_data,
      "and instrument from", start_date_instrument, "until", end_date_instrument, "\n")
  cat("Matched", num_matches, "observations: from ", matched_display, "\n")

  return(invisible(x))
}


