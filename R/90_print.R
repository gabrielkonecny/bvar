
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
print.check_iv <- function(out, ...) {

  # Extract the start and end dates from the data and instrument date ranges
  start_date_data <- out$dates_residuals[1]
  end_date_data <- tail(out$dates_residuals, 1)
  start_date_instrument <- out$dates_instrument[1]
  end_date_instrument <- tail(out$dates_instrument, 1)

  if(out$manual_matching == TRUE){
    print(class(out$residuals_shortened))
    cat("Performing exact matching... Note that this disregards any information
        on the indexing provided in rownames.", "\n")
    cat("Head of cbind(residuals, instrument):", "\n")
    print(head(cbind(out$residuals_shortened, out$instrument_shortened)))
    # cat("Found residuals:", head(irf$residuals_shortened,1), "...", tail(out$residuals_shortened,1),
    #     "and instrument:",  head(out$instrument_shortened,1), "...",  tail(out$instrument_shortened,1))
  } else{

  # Get the number of matched dates and display the first and last matched observations
  num_matches <- length(out$common_dates)
  matched_display <- paste0(
    out$common_dates[1], " (value ", format(out$instrument_shortened[out$common_dates[1]], scientific = TRUE),
    ") until ",
    out$common_dates[num_matches], " (value ", format(out$instrument_shortened[out$common_dates[num_matches]], scientific = TRUE),
    ")"
  )

  # Print the summary information
  cat("Found residuals from", start_date_data, "until", end_date_data,
      "and instrument from", start_date_instrument, "until", end_date_instrument, "\n")
  cat("Matched", num_matches, "observations: from ", matched_display, "\n")

   }

  return(invisible(out))
}


