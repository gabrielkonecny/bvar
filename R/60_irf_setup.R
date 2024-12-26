
#' Impulse response settings and identification
#'
#' Provides settings for the computation of impulse responses to
#' \code{\link{bvar}}, \code{\link{irf.bvar}} or \code{\link{fevd.bvar}}. Allows
#' setting the horizon for which impulse responses should be computed, whether
#' or not forecast error variance decompositions (FEVDs) should be included
#' as well as if and what kind of identification should be used. See the Details
#' section for further information on identification. Identification can be
#' achieved via Cholesky decomposition, sign restrictions (Rubio-Ramirez,
#' Waggoner and Zha, 2010), combination of zero and sign restrictions (Arias,
#' Rubio-Ramirez and Waggoner, 2018) or using an external instrument
#' (Stock and Watson 2012, Mertens and Ravn 2013)
#'
#' Identification can be performed via Cholesky decomposition, sign
#' restrictions, zero and sign restrictions or using external instrument.
#' The algorithm for generating suitable sign restrictions follows
#' Rubio-Ramirez, Waggoner and Zha (2010), while the one for zero and sign
#' restrictions follows Arias, Rubio-Ramirez and Waggoner (2018).
#' Note the possiblity of finding no suitable zero/sign restrictions.
#' The proxy SVAR identification follows the implementation in
#' Miranda-Agrippino and Ricco (2021).
#'
#' @param horizon Integer scalar. The horizon for which impulse responses
#' (and FEVDs) should be computed. Note that the first period corresponds to
#' impacts i.e. contemporaneous effects.
#' @param fevd Logical scalar. Whether or not forecast error variance
#' decompositions should be calculated.
#' @param identification Logical scalar. Whether or not the shocks used for
#' calculating impulses should be identified. Defaults to \code{TRUE}, i.e.
#' identification via Cholesky decomposition of the VCOV-matrix unless
#' \emph{sign_restr} or \emph{instrument} is provided.
#' @param sign_restr Elements inform about expected impacts
#' of certain shocks. Can be either \eqn{1}, \eqn{-1} or \eqn{0} depending
#' on whether a positive, a negative or no contemporaneous effect of a
#' certain shock is expected. Elements set to \eqn{NA} indicate that there are
#' no particular expectations for the contemporaneous effects. The default
#' value is \code{NULL}. Note that in order to be fully identified at least
#' \eqn{M * (M - 1) / 2} restrictions have to be set and a maximum of
#' \eqn{M - j} zero restrictions can be imposed on the \eqn{j}'th column.
#' @param sign_lim Integer scalar. Maximum number of tries to find suitable
#' matrices to for fitting sign or zero and sign restrictions.
#' @param instrument Numeric vector. If provided, the identification is performed using proxy
#' SVAR. Since the residuals and instrument need not to have equal length in general,
#' user is expected to provide a common index (rownames) for the data and the
#' instrument. See examples and helper function \emph{set_dates}. This can be
#' avoided by setting \emph{manual_matching} to TRUE.
#' @param proxyvar character string. Variable for which instrument is provided -
#' Needs to correspond to one of the column names of the data used in
#' \emph{bvar}.
#' @param manual_matching If set to TRUE, user is not expected to specify common
#' index for the data and instrument. Instead the length of instrument needs to
#' match the length of residuals (= length(data inputted in bvar) - lags).
#'
#' @return Returns a named list of class \code{bv_irf} with options for
#' \code{\link{bvar}}, \code{\link{irf.bvar}} or \code{\link{fevd.bvar}}.
#'
#' @references
#'   Rubio-Ramirez, J. F. and Waggoner, D. F. and Zha, T. (2010) Structural
#'   Vector Autoregressions: Theory of Identification and Algorithms for
#'   Inference. \emph{The Review of Economic Studies}, \bold{77}, 665-696,
#'   doi:10.1111/j.1467-937X.2009.00578.x.
#'   Arias, J.E. and Rubio-Ramirez, J. F. and Waggoner, D. F. (2018)
#'   Inference Based on Structural Vector Autoregressions Identifiied with
#'   Sign and Zero Restrictions: Theory and Applications.
#'   \emph{Econometrica}, \bold{86}, 2, 685-720,
#'   doi:10.3982/ECTA14468.
#'   Miranda-Agrippino, S., & Ricco, G. (2021). The transmission of monetary
#'   policy shocks. \emph{American Economic Journal: Macroeconomics},
#'    \bold{13(3)}, 74-107, doi:10.1257/mac.20180124
#'    Stock, J. H., and Watson, M. W. (2012). Disentangling the Channels of the
#'    2007-2009 Recession (No. w18094) \emph{National Bureau of Economic
#'    Research}, doi:10.3386/w18094.
#'    Mertens, K., and Ravn, M. O. (2013). The dynamic effects of personal and
#'    corporate income tax changes in the United States. \emph{American
#'    economic review}, \bold{103(4)}, 1212-1247,
#'    doi:10.1257/aer.103.4.1212.
#'
#' @seealso \code{\link{irf.bvar}}; \code{\link{plot.bvar_irf}}
#'
#' @keywords BVAR irf fevd settings
#'
#' @export
#'
#' @examples
#' # Set impulse responses to a horizon of 20 time periods and enable FEVD
#' # (Identification is performed via Cholesky decomposition)
#' bv_irf(horizon = 20, fevd = TRUE)
#'
#' # Set up structural impulse responses using sign restrictions
#' signs <- matrix(c(1, NA, NA, -1, 1, -1, -1, 1, 1), nrow = 3)
#' bv_irf(sign_restr = signs)
#'
#' # Set up structural impulse responses using zero and sign restrictions
#' zero_signs <- matrix(c(1, 0, NA, -1, 1, 0, -1, 1, 1), nrow = 3)
#' bv_irf(sign_restr = zero_signs)
#'
#' # Set up structural impulse responses using external instrument with length
#' corresponding to reduced form residuals
#' bv_irf(instrument = istrument)
#'
#' # Prepare to estimate unidentified impulse responses
#' bv_irf(identification = FALSE)
bv_irf <- function(
  horizon = 12,
  fevd = FALSE,
  identification = TRUE,
  sign_restr = NULL,
  sign_lim = 1000,
  instrument = NULL,
  proxyvar = NULL,
  manual_matching = FALSE
  )  {

  # Input checks
  horizon <- int_check(horizon, min = 1, max = 1e6,
    msg = "Invalid value for horizon (outside of [1, 1e6]).")
  sign_lim <- int_check(sign_lim, min = 100, max = Inf,
    msg = "Invalid value for sign_lim (outside of [100, Inf]).")

  if(!is.logical(c(identification, fevd))){
    stop("Please provide fevd and identification as logical scalars.")
  }

  zero <- FALSE # Zero or sign restrictions
  if(identification) {
    if(!is.null(sign_restr)) {
      restr_len <- length(sign_restr)
      if(!is.numeric(sign_restr) && !all(sign_restr %in% c(-1, 0, NA, 1)) &&
        sqrt(restr_len) %% 1 != 0) {
        stop("Please provide sign_restr as a numeric square matrix ",
          "containing NAs, 1s and -1s (and 0s for zero restrictions).")
      }
      if(0 %in% sign_restr) {zero <- TRUE}
      if(is.vector(sign_restr)) {
        sign_restr <- matrix(sign_restr, nrow = sqrt(restr_len))
      }
      if(zero && any(colSums(sign_restr == 0, na.rm = TRUE) >
        rev(seq_len(sqrt(restr_len)) - 1))) {
        stop("Number of zero restrictions on at least one of the shocks is ",
          "too high. Please reduce or change the order of variables.")
      }

      if(sum(!is.na(sign_restr)) <
        (sqrt(restr_len) - 1) * sqrt(restr_len) / 2) {
        message("Number of restrictions implies an underidentified system.")
      }
    }
    # Cholesky

  }

  # IV

  if(!is.null(instrument)){
    if (is.numeric(instrument)) {
    } else {
      stop("Input must be a numeric vector. Multiple instruments are not
           supported.")
    }
  }

  if(!manual_matching && !is.null(instrument) &&
     is.null(names(instrument))){
    stop("No index detected for the object instrument. If manual_matching is
  FALSE, names(instrument) need to be specified. For identification to work
  correctly a common index for the data and the instrument must be specified.
  Alternatively, switch manual_matching to TRUE.")
  }

  if (!is.null(instrument) && is.null(proxyvar)) {
    stop("If 'instrument' is specified, 'proxyvar' must also be specified.")
  }

  if(!is.null(proxyvar) && !is.character(proxyvar)){
    stop("Proxyvar needs to be an object of type character.")
  }

  # Outputs
  out <- list("horizon" = horizon, "fevd" = fevd,
    "identification" = identification,
    "sign_restr" = sign_restr, "zero" = zero,
    "sign_lim" = sign_lim,
    "instrument" = instrument,
    "manual_matching" = manual_matching,
    "proxyvar" = proxyvar
  )

  class(out) <- "bv_irf"

  return(out)
}
