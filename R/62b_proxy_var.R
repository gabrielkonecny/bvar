#' Proxy SVAR
#'
#' Performs identification using external \emph{instrument}. At this stage the
#' \emph{instrument} and \emph{residuals} need to have equal number of
#' observations. Transformations such as intersecting the instrument and
#' residuals need to be done before calling this function.
#'
#' @param residuals Numeric matrix. During estimation, draw specific residuals
#' are used.
#' @param instrument Numeric vector containing the instrument. Possibly
#' shortened beforehand to match the length of residuals.
#'
#' @return Returns a list with the following elements:
#' \itemize{
#'   \item \code{impact} - Contemporaneous transmission coefficients: Bzero
#'   \item \code{Gamma}  - Estimated correlation between shock and instrument
#'   \item \code{f_stat} - First stage F statistic
#'   \item \code{lambda} - Reliability of instrument
#' }
#' @noRd

proxy_svar <- function(residuals, instrument){

  residuals <- as.matrix(residuals)

  t <- nrow(residuals)
  n <- ncol(residuals)
  m <- 1 #only m=1 instrument implemented at the moment

  # Create the Kronecker product
  eye_n <- diag(rep(1, n)) #Matrix(diag(rep(1, n)), sparse = TRUE)
  tempX <- cbind(1, instrument)
  kron_matrix <- kronecker(eye_n, tempX)

  # Solve the linear equation
  betaIV <- solve(t(kron_matrix) %*% kron_matrix) %*% t(kron_matrix) %*% as.vector(residuals)

  # Reshape betaIV
  betaIV <- matrix(betaIV, nrow = length(betaIV) / n, byrow = FALSE)

  # Transpose to get the final result
  betaIV <- as.matrix(t(betaIV))

  # F stat (regression on instruments of relevant innovations)
  proxyVar <- instrument #Using notation from Agrip Ricco for this
  tempU <- residuals[, 1:m] - tempX %*% betaIV[1:m, ]

  tempY <- tempX %*% betaIV[1:m, ] - matrix(rep(mean(residuals[, 1:m]), t), ncol = m, byrow = TRUE)
  k <- length(betaIV[1:m, ]) - 1

  f_stat <- ((t(tempY) %*% tempY) / k) / ((t(tempU) %*% tempU) / (t - k - 1))


  ################
  # Assuming m is defined and betaIV is already a numeric matrix [1:5, 1:2]

  # beta_11 and beta_21
  beta_11 <- betaIV[1:m, 2:(m+1)]
  beta_21 <- betaIV[(m + 1):nrow(betaIV), 2:(m+1)]

  # ratio of regression coefficients
  B21B11 <- beta_21 / beta_11

  # Covariance matrix
  SigmaU <- cov(residuals)

  # Identification
  Zeta <- (B21B11 * SigmaU[1:m, 1:m]) %*% t(B21B11) -
    (SigmaU[(m + 1):nrow(SigmaU), 1:m] %*% t(B21B11) + B21B11 %*% t(SigmaU[(m + 1):nrow(SigmaU), 1:m])) +
    SigmaU[(m + 1):nrow(SigmaU), (m + 1):ncol(SigmaU)]

  B12B12 <- t(SigmaU[(m + 1):nrow(SigmaU), 1:m] - B21B11 * SigmaU[1:m, 1:m]) %*% solve(Zeta) %*%
    (SigmaU[(m + 1):nrow(SigmaU), 1:m] - B21B11 * SigmaU[1:m, 1:m])

  B11B11 <- SigmaU[1:m, 1:m] - B12B12

  #Since we only have 1 instrument
  B11 <- sqrt(B11B11)   # beta_{11}
  B <- as.vector(B11) * c(1, B21B11)  # first column of B (u_t = B * e_t)

  ###########################
  #%realized shock sequences (Montiel-Olea, Stock and Watson)
  # Assuming T, residuals, and proxyVar are defined

  # Create tempX matrix
  tempX <- cbind(1, residuals)

  # Calculate e using matrix multiplication and solving the system of equations
  e <-  tempX %*% solve(t(tempX) %*% tempX) %*% t(tempX) %*% proxyVar

  # Standardize e to have unit variance
  e <- scale(e, center = TRUE, scale = apply(e, 2, sd)) # Now, e contains the unit variance shock series

  ##################
  #Again since m=1

  # Assuming proxyVar, residuals, and e are defined in R

  # Proportion of uncensored data
  D <- (proxyVar != 0)

  # Covariance matrix
  SigmaMU <- cov(cbind(proxyVar[D], residuals[D, ]), use = "pairwise")

  # Relevance
  A_ <- SigmaMU[1,2:ncol(SigmaMU)]
  Phi <- 1/(MASS::ginv(A_) %*% B)
  G <- 1/(t %*% MASS::ginv(sum(D)))
  Gamma <- (Phi*(1/G))

  # Calculate eSquare and zSquare
  eSquare <- e^2
  zSquare <- (proxyVar - rep(Gamma,t) * as.vector(e))^2

  # Calculate Lambda
  Lambda <- 1/((Gamma^2 * sum(eSquare[D]) + sum(zSquare[D])) / (Gamma^2 * sum(eSquare[D])))


  ####################
  iP <- 1 #instrument position
  # Load output structure

  impact = matrix(NA, n, m)      # contemporaneous transmission coefficients: Bzero
  impact[iP,] <- B[1:m]
  impact[(iP+1):n,] <- B[(m + 1):n]

  #normalize
  impact <- impact / impact[1]


  #Summary
  # impact = impact,             # contemporaneous transmission coefficients: Bzero
  # Gamma = Gamma,             # estimated correlation between shock and instrument
  # L = Lambda,                # reliability of instrument
  # e = e,                     # realized shocks series
  # fstat = diag(f_stat)       # F statistic of regression on instrument

  output <- list(impact, Gamma, f_stat, Lambda)
  names(output) <- c("impact", "Gamma", "f_stat", "Lambda")

  return(output)
}


