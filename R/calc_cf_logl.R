#' Calculate/Simulate model counterfactual in matrix one sector form
#'

#' \code{calc_cf_log_m_1} is a function that calculates the
#' model for one sector assumption for a given set of data and a
#' shock. This function gets sample data in matrix form
#' ( either arbitrary matrix or real world data)
#' , shocks, and parameters and returns changes of welfare C in two forms:

#' 1. C_hat_1 returns results based on X' computed based on psuedo inverse
#' 2. C_hat_2 returns result based on X' computed based on normalized Q.
#' At the end, C_hat_1 and C_hat_2 should convert to a same number with
#' with an infinitesimal difference


#' @details This script calculates welfare changes as a result of shock to one sector model:                                   #
#' Here we are using matrix representation of the mode.

#' using options() function we set digits=22 to make sure about the accuracy of the model
#' of  convergence.

#' The objective of using option function is to make sure
#' about the convergence of X' using psuedo inverse in comparison
#' with inverse of normalized one using Q matrix which is
#' Q=gamma*R as it is mentioned in footnote-page 10-task 1.
#'
#' @param data a list of model variables given as examples
#'   \describe{
#'   \
#'   \item{R}{matrix of country revenues with
#'            \code{J = length(country) columns}}
#'   \item{pi}{vector of import shares across locations for each
#'            destination combination, where the import share of
#'             origin \code{i} products used in destination
#'             \code{j}
#'   \item{gamma}{gamma with items gamma_j vector of intermediate cost shares of each country countries for each
#'                    . after receiving this vector as input
#'                    it is going to be turned into diagonal matrix J*J
#' @param shock a list containing a shock with values
#'   \describe{
#'     \item{T_hat}{technological shock}}
#'     \item{tau_hat}{tariff shocks that may affect both importer and exporter}
#' @return after computing all X,X',A,Ap,Y1,Y1p,Y2,Y2p,B,Bp,Z1,Z1p,Z2,Z2p
#' Z3,Z3p,Z4p,Y3, the function returns C_hat_1 and C_hat_2 based on psuedo and normalized X'
#' @export calc_cf_logl

calc_cf_logl <- function(J, R, pi_I, pi_F, gamma, T_hat, tau_hat_I,
                            tau_hat_F, epsilon, use_pseudo = FALSE) {
  # options(digits=22) this only changes the number of digits that R prints but
  # not precision. Without extra packages you are limited to ~15-16 digits of
  # precision

  diaggamma <- diag(gamma)
  IJ <- diag(J)
  Igammadiag <- IJ - diaggamma
  pi_Igamma <- crossprod(pi_I, Igammadiag)
  pi_Fgamma <- crossprod(pi_F, Igammadiag)

  # step 1 Matrix X calc
  X <- solve(IJ - pi_Igamma)   ##### matrix X def checked

  # step 2 Matrix A calc
  # getting the reciprocal of R vector (1/Ri)
  # R defined as JxS matrix, with S = 1 its a column not a row vector!
  R <- t(R)
  Rrecip <- 1/R
  A <- crossprod(Rrecip, R) * (pi_I %*% Igammadiag)

  # step 3 Matrix diagZ calc
  # matrix(diag(z))
  diagz <- diag(as.vector(epsilon * Rrecip * (R %*% t(pi_I %*% Igammadiag))))

  # step 4 Matrix A' calc
  Ap <- crossprod(Rrecip, R) * (pi_F %*% diaggamma)

  # step 5 Matrix Y2 calc
  Y2 <- crossprod(pi_F, diaggamma +
                    crossprod(Igammadiag, crossprod(t(X),
                                                    crossprod(pi_I, diaggamma)
        )))

  # step 6 Matrix Z' calc
  #####matrix diag(z')
  diagzp <- diag(as.vector(epsilon * Rrecip * (R %*% t(pi_F %*% diaggamma))))

  # step 7 Matrix Z2 calc
  ddzx <- diagz * Igammadiag * (X %*% t(pi_I))
  ddzxp <- diagzp * Igammadiag * (X %*% t(pi_I))
  zg <- diagz * diaggamma
  zIg <- diagz * Igammadiag
  zpg <- diagzp %*% diaggamma
  zpIg <- diagzp %*% Igammadiag

  #######matrix Z2 correct
  Z2 <- A - zg + (epsilon * A - zIg) %*% X %*% t(pi_I) %*% diaggamma

  # step 8 Matrix Z2' calc
  Z2P <- Ap + epsilon * (Ap %*% Y2) - zpg - zpIg %*% X %*% t(pi_I) %*% diaggamma

  # step 9 Matrix Z1 calc
  Z1 <- 1 / epsilon * (diagz - (epsilon * A - zIg) %*% X %*% t(pi_I))

  # step 10 Matrix Y1 calc
  Y1 <- -1 / epsilon *
    (t(pi_F) +
       crossprod(pi_F, crossprod(Igammadiag, crossprod(t(X), t(pi_I)))))

  # step 11 Matrix Z1' calc
  Z1p <- 1 / epsilon * (diagzp + (zpIg %*% X %*% t(pi_I))) + epsilon * Ap %*% Y1

  # step 12 Matrix B calc
  B <- epsilon * crossprod(Rrecip, R) *
    (pi_I %*% Igammadiag %*% diag(colSums(pi_I)))

  # step 13 Matrix B' calc
  Bp <- epsilon * crossprod(Rrecip, R) *
    (pi_F %*% diaggamma %*% diag(colSums(pi_F)))

  # step 14 Matrix Z3 calc
  Z3 <- (epsilon * A - zIg) %*% X

  # step 15 Matrix Y3 calc
  Y3 <- crossprod(t(pi_Fgamma), X)

  # step 16 Matrix Z3' calc
  Z3p <- epsilon * (Ap %*% Y3) - (zpIg) %*% X

  # step 17 Matrix Z4' calc
  Z4p <- epsilon * Ap

  #comparison of X' using psuedo inverse with normalized format
  # step 18 Matrix X' calc

  ######### matrix X'? using normalized def with Q
  gR <- R * gamma
  Q <- matrix(rep(gR, each = J), ncol = J, byrow = FALSE)

  if(use_pseudo) {
    ######### matrix X'? using psuedo inverse
    Xp <- MASS::ginv(IJ - Z2 - Z2P)
  } else {
    Xp <- solve(IJ - (Z2 + Z2P) - Q)
  }

  XZ1 <- Xp %*% (Z1 + Z1p)
  XZ3 <- Xp %*% (Z3 + Z3p)

  # return C_hat
  (XZ1 - (Y1 + Y2 %*% XZ1)) %*% T_hat +
      (Y2 %*% Xp - Xp) %*% diag(B %*% t(tau_hat_I)) +
      (Y2 %*% Xp - Xp) %*% diag(Bp %*% t(tau_hat_F)) +
      (XZ3 - (Y3 + (Y2 %*% XZ3))) %*% diag(crossprod(pi_I, tau_hat_I)) +
      ((Xp %*% Z4p) - (IJ + Y2 %*% Xp %*% Z4p)) %*%
      diag(crossprod(pi_F, tau_hat_F))

}
