#' Impulse response draws
#'
#' Computes impulse responses using the posterior draws of the VAR coefficients
#' and VCOV-matrix obtained from \code{\link{draw_post}}.
#'
#' @param beta_comp Numeric matrix. Posterior draw of the VAR coefficients in
#' state space representation.
#' @param sigma Numeric matrix. Posterior draw of the VCOV-matrix of the
#' model.
#' @param sigma_chol Numeric matrix. Lower part of the Cholesky decomposition
#' of \emph{sigma}. Calculated as \code{t(chol(sigma))}.
#' @param M Integer scalar. Number of columns in \emph{Y}.
#' @param lags Integer scalar. Number of lags in the model.
#' @param horizon Integer scalar. Horizon for which impulse responses should be
#' computed. Note that the first period corresponds to impacts i.e.
#' contemporaneous effects.
#' @param identification Logical scalar. Whether or not the shocks used for
#' calculating the impulse should be identified. Defaults to \code{TRUE},
#' meaning identification will be performed recursively through a
#' Cholesky decomposition of the VCOV-matrix as long as \emph{sign_restr}
#' and \emph{zero_restr} are \code{NULL}. If set to \code{FALSE}, shocks will
#' be unidentified.
#' @param sign_restr  Numeric matrix. Elements inform about expected impacts
#' of certain shocks. Can be either \eqn{1}, \eqn{-1} or \eqn{0} depending
#' on whether a positive, a negative or no contemporaneous effect of a
#' certain shock is expected. Elements set to \eqn{NA} indicate that there are
#' no particular expectations for the contemporaneous effects.
#' @param zero Logical scalar. Whether to impose zero and sign restrictions,
#' following Arias et al. (2018).
#' @param sign_lim Integer scalar. Maximum number of rotational matrices to
#' draw and check for fitting sign restrictions.
#'
#' @return Returns a numeric array of impulse responses.
#'
#' @noRd
compute_irf <- function(
  beta_comp,
  sigma, sigma_chol,
  M, lags,
  horizon,
  identification,
  sign_restr, zero = FALSE, sign_lim = 10000,
  residuals = NULL, instrument = NULL, manual_matching = FALSE, proxyvar = NULL) {



  # Identification
  if(identification) {
    sigma_chol <- t(chol(sigma))
    if(is.null(sign_restr) & is.null(instrument)) {
      shock <- sigma_chol
    }
    if(!is.null(sign_restr) & is.null(instrument)){
    shock <- sign_restr(sigma_chol = sigma_chol,
    sign_restr = sign_restr, M = M, sign_lim = sign_lim, zero = zero)
    }
    if(is.null(sign_restr) & !is.null(instrument)){
          shock <- diag(M)

          # Find the column index of the proxy variable
          col_index <- which(colnames(residuals) == proxyvar)

          # If the proxy variable is not already the first column, reorder columns
          if (col_index != 1) {
            # Instrumented variable swaps place with first variable in system
            proxy_svar_ordering <- get_swapped_index(residuals, proxyvar)

            # For proxy svar the instrumented variable is now ordered first
            proxy_svar_output <- proxy_svar(
              residuals[, proxy_svar_ordering],instrument)

            # Reverse the ordering for output
            proxy_svar_output$impact <- proxy_svar_output$impact[proxy_svar_ordering, , drop = FALSE]
            shock[,col_index] <- proxy_svar_output$impact

            } else{


          # Assumed instrumented variable is ordered first
          proxy_svar_output <- proxy_svar(residuals,instrument)
          shock[,1] <- proxy_svar_output$impact
            }
    }
    if(!is.null(sign_restr) & !is.null(instrument)){
    stop("Sign restrictions and instrument cannot be used at the same time!")
    }
  } else {shock <- sigma}

  # Impulse responses
  irf_comp <- array(0, c(M * lags, horizon, M * lags))
  irf_comp[1:M, 1, 1:M] <- shock
  for(i in 2:horizon) {
    irf_comp[, i, ] <- beta_comp %*% irf_comp[, i - 1, ] # Could vectorise
  }
  irf_comp <- irf_comp[1:M, , 1:M]

  output <- list()
  output$irf_comp <- irf_comp


  output$iv_f_stat <- if(!is.null(instrument)){
    proxy_svar_output$f_stat} else{NULL}

  return(output)
}
