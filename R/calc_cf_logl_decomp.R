#' This file is a function,calc_cf_log_t2  that represents task 2 updated
#'  version.
#' The basic idea behind this task is to rewrite our matrix as (Neumann)
#' power series:
#' (I-M)^-1  <-  I +M+M^2 +M^3 + ..., (20)
#' where, for brevity,it is referred to I as a 0-degree effect, M as a 1st
#'  degree effect, M2 as a 2nd degree effect, etc. Importantly,
#'  the above decomposition
#' only works if the spectral radius of M, (M), is below one, i.e., (M) < 1,
#'  which we will have to verify below.
#' To conduct this decomposition exercise in our framework, we rewrite
#' the above matrices in the following block matrix form:
#'c(dlogPI,dlogR) <- V(I-M)^-1 and if pho(M) which is the spectral
#'radius is below one,
#'
#'
#'However, in larger set of data sets it seems that does not converge.
#'
#'Thus, we need rewriting the M term. We seperate matrix G^-1 and multiply
#' it to avoid the step of the computation of spectral radius for convergence
#'  insurance objective.


#'we can rewrite it as summation I + M_new + M^2_new +....... .
#'we are going to find dlogC based on summation and V(I-M)^-1 and compare
#'the results to the output of task1 function which is calc_cf_log_m_2 function.
#' \code{calc_cf_log_t2} is a function that calculates the
#' model for one sector assumption for a given set of data and a
#' shock. This function gets sample data in a matrix form
#' ( either arbitrary matrix or real world data)
#' , shocks, and parameters and returns changes of welfare C based on V(I-M)^-1
#' and summation of I+M+M2+..... It also shows all degree effects separately.

#' 1. C_hat_inv returns results based on V(I-M)^-1.
#' 1. C_hat_inv_new returns results based on VG^(-1)*(I-M_new)^-1.
#' 2. C_hat_sum returns result based on VG^(-1)(I+M_new+M^2_new+M^3+_new.....)
#' At the end, C_hat_inv and C_hat_sum should convert to a same number with
#'  an infinitesimal difference.
#'
#'
#' The question is: Where to stop the summation I+M_new+M^2_new...?
#' We set a while loop to make sure we have calculated enough terms
#' and the final answers converge.
#' While loop criteria is to minimize
#' sum(square(V(I-M)^-1)-(I+M2+M3+.....)(i,j))
#' It can be computed based on other criteria as well if we want.


#' @details This script calculates welfare changes as a result of
#'  shock to one sector model:                                   #
#' Here we are using matrix representation of the mode.

#' using options() function we set digits <- 22 to make sure
#' about the accuracy of the model
#' of  convergence.

#' The objective of using option function is to make sure
#' about the convergence of summation
#'
#' @param data a list of model variables given as examples
#'   \describe{
#'   \
#'   \item{R}{matrix of country revenues with
#'            \code{J  <-  length(country) columns}}
#'   \item{pi}{vector of import shares across locations for each
#'            destination combination, where the import share of
#'             origin \code{i} products used in destination
#'             \code{j}
#'   \item{gamma}{gamma with items gamma_j vector of intermediate cost shares
#'    of each country countries for each
#'                    . after receiving this vector as input
#'                    it is going to be turned into diagonal matrix J*J
#' @param shock a list containing a shock with values
#'   \describe{
#'     \item{T_hat}{technological shock}}
#'     \item{tau_hat}{tariff shocks that may affect both importer and exporter}
#' @return after computing all necessary variables,The function returns,V(I-M),
#'  I+M+M2.....,the difference,
#'  dlogc based on V(I-M^-1 and dlogC based on summation.
#'  Otherwise it says:  "The spectral radius is above one
#'  the decomposition cannot be done."
#' and individual degree effects of M M2 ......
#' It also shows the number of repetitions that convergence happens.


#' calc_cf_log_t2 -> function(simple_data,shocks,parameters)

#' three blocks of input:

#' 1. simple_data <-
#' R revenue or GDP vector of size J
#' J number of countries
#' pi <- as.vector , is a vector of size 2*J*J
#' needed to be divided into intermediate and final pi, pi_I:J*J for
#'  intermediate
#' pi_F:J*J for final
#' pi_I <- matrix(pi[1:J*J],nrow <- J,ncol <- J)
#' pi_F <- matrix(pi[(J*J)+1:J*J+J*J],nrow <- J,ncol <- J)
#' gamma <- as.vector of size J



#'2. parameters <-
#'   epsilon



#' 3. shock <-
#'    list(
#'      T_hat  <-  matrix(vector(c(), nrow  <-  J,ncol <- 1)
#'      tau_hat <- matrix(c(),nrow <- 1,ncol <- 2*J*J)) should be divided
#'       into two types of shocks: intermediates and final
#'      tau_hat_I <- matrix(tau_hat[1:(J*J)],nrow <- J,ncol <- J)
#'      tau_hat_F <- matrix(tau_hat[(J*J)+1:(J*J)+(J*J)],nrow <- J,ncol <- J)


#'
#' Outputs are  return(list(c_hat_0,c_new_0,C_hat_inv_msg,
#'  C_hat_inv_new,C_hat_sum_new, nth_msg))
#'
#'  zero effect of version 1 and two, c_hat and c_hat_new based on inverse term
#'  c_hat_new_ sum c_hat based on sum terms
#'  terms of summation polynomial

calc_cf_logl_decomp <- function(J, R, pi_I, pi_F, gamma, T_hat, tau_hat_I,
                                tau_hat_F, epsilon, ktol = 1e-10,
                                kmax_iter = 50) {

  R <- matrix(R, nrow = 1)
  diaggamma <- diag(gamma)
  Igammadiag <- diag(J) - diaggamma
  pi_Igamma <- crossprod(pi_I, Igammadiag)
  pi_Fgamma <- crossprod(pi_F, Igammadiag)
  tpi_I <- t(pi_I)
  tpi_F <- t(pi_F)
  Rrecip <- 1/R

  ################## B & Bp caculations
  B <- epsilon * crossprod(Rrecip, R) *
    (pi_I %*% Igammadiag %*% diag(c(colSums(pi_I))))

  Bp <- epsilon * crossprod(Rrecip,R) *
    (pi_F %*% diaggamma %*% diag(c(colSums(pi_F))))

  #  Matrix A&Ap calculation ######

  A <- crossprod(Rrecip, R) * (pi_I %*% Igammadiag)
  Ap <- crossprod(Rrecip, R) * (pi_F %*% diaggamma)

  Q <-  matrix(rep(R / sum(R), each = J), ncol = J)

  ######################################################################
  ######step3: Creating Block Matrices: computing matrix M #############
  ######################################################################

  eA <- epsilon * A
  eAp <- epsilon * Ap

  M1 <- crossprod(pi_I, Igammadiag)
  M2 <- crossprod(pi_I, diaggamma)
  M3 <- eA %*% M1 + eAp %*% pi_Fgamma - epsilon * Igammadiag
  M4 <- A + Ap + eA %*% tpi_I %*% diaggamma + eAp %*% tpi_F %*% diaggamma -
    epsilon * diaggamma - Q

  M <- matrix(rbind(cbind(M1, M2), cbind(M3,M4)), ncol = 2*J)

  ############################################################################
  #####################                            ###########################
  ##################### Step 5: computing V matrix ###########################
  #####################                            ###########################
  ############################################################################

  vv1 <- -1 / epsilon * tpi_I
  vv2 <- diag(J) - A %*% tpi_I - Ap %*% tpi_F
  V1 <- rbind(vv1, vv2)

  vvv1 <- matrix(diag(tpi_I %*% tau_hat_I), nrow = J)
  vvv2 <- matrix(eA %*% diag(tpi_I %*% tau_hat_I) +
                   eAp %*% diag(tpi_F %*% tau_hat_F) -
                   diag(B %*% t(tau_hat_I)) - diag(Bp %*% t(tau_hat_F)),
                 nrow = J)
  V2 <- rbind(vvv1, vvv2)

  V <- V1 %*% T_hat + V2

  #######################################################################
  #######################################################################
  ###################                                               #####
  ################### step 6: Computing C_hat based on V(I-M)^(-1) #####
  ###################                                               #####
  #######################################################################
  #######################################################################

  z <- solve(diag(nrow(M)) - M)

  pifg <- -pi_Fgamma
  pifI <- diag(J) - tpi_F %*% diaggamma
  pigamma <- cbind(pifg, pifI)

  S1 <- -diag(crossprod(pi_F, tau_hat_F))
  S2 <- 1 / epsilon * tpi_F %*% T_hat

  C_hat_inv <- S1 + S2 + pigamma %*% z %*% V
  C_hat_0 <- S1 + S2 + pigamma %*% diag(2*J) %*% V

  ##########################################################################
  #########################  Step 7:                                   #####
  #########################  computing C_new_inv and C_sum_new         #####
  #########################      using G_inv                           #####
  ##########################################################################

  zer_o <- matrix(0, ncol = J, nrow = J)
  G <- rbind(cbind(diag(J), zer_o),
             cbind(epsilon * Igammadiag, diag(J) + epsilon*diaggamma))
  G_inv <- solve(G)

  M1_n <- crossprod(pi_I, Igammadiag)
  M2_n <- crossprod(pi_I, diaggamma)
  M3_n <- (A %*% tpi_I + Ap %*% tpi_F) %*% (epsilon * Igammadiag)
  M4_n <- A + Ap + eA%*% tpi_I %*% diaggamma + eAp %*% tpi_F %*% diaggamma - Q
  M_n <- matrix(rbind(cbind(M1_n, M2_n), cbind(M3_n, M4_n)), ncol = 2*J)
  M_new <- G_inv %*% M_n

  if(max(abs(eigen(M_new)$value)) > 1) {
    stop("Spectral radius is larger 1 with rho = ",
         max(abs(eigen(M_new)$value)))
  }

  #######################################################################
  ###################*   and stack all in a list               ##########
  ###################*   and add all of them together          ##########
  ###################*   to get (I+M_new+M_new^2+....)         ##########
  ###################*   while the difference of               ##########
  ###################*   sumsquare(V G^-1(I-M)^(-1)            ##########
  ###################*    -(I+M+M2+....)) <- converg3          ##########
  #######################################################################

  z_new <- solve(diag(nrow(M_new)) - M_new)

  # choosing a large number for while loop

  tol <-  1
  i <- 0
  M_sum <- matrix(0, nrow = nrow(M_new), ncol = ncol(M_new))
  decomposition <- as_tibble(
    matrix(0, nrow = J, ncol = kmax_iter),
    .name_repair = ~paste0("effect_", 0:(kmax_iter - 1))
  )

  while(tol > ktol & i < kmax_iter) {
    if (i  ==  0){
      M_P  <-  diag(1, 2 * J)
      decomposition[,i + 1] <- S1 + S2 + pigamma %*% M_P %*% G_inv %*% V
    }
    else {
      M_P <- M_P %*% M_new
      decomposition[,i + 1] <- pigamma %*% M_P %*% G_inv %*% V
    }
    M_sum <- M_sum + M_P
    tol <- max(abs(z_new - M_sum))
    i <- i + 1
    print(paste0("Iteration, ", i, " tolerance reached: ", tol))
  }

  C_hat_inv_new <- S1 + S2 + pigamma %*% z_new %*% G_inv %*% V
  C_hat_sum_new <- S1 + S2 + pigamma %*% M_sum %*% G_inv %*% V

  bind_cols(
    tibble("Welfare0" = C_hat_0[,1],
           "Welfare" = C_hat_inv_new[,1],
           "Welfare_sum_approximation" = C_hat_sum_new[,1]),
    decomposition
  )
}
