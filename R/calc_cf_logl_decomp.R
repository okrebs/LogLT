#' This file is a function,calc_cf_log_t2  that represents task 2.
#' The basic idea behind this task is to rewrite our matrix as (Neumann) power series:
#' (I???M)???1 = I +M+M^2 +M^3 + ..., (20)
#' where, for brevity,it is referred to I as a 0-degree effect, M as a 1st degree effect, M2 as a 2nd degree effect, etc. Importantly, the above decomposition
#' only works if the spectral radius of M, (M), is below one, i.e., (M) < 1, which we will have to verify below.
#' To conduct this decomposition exercise in our framework, we rewrite the above matrices in the following block matrix form:
#'c(dlogPI,dlogR)=V(I-M)^-1 and if pho(M) which is the spectral radius is below one,
#'we can rewrite it as summation I+M+M^2....... .
#'we are going to find dlogC based on summation and V(I-M)^-1 and compare
#'the results to the output of task1 function which is calc_cf_log_m_2 function.

#' \code{calc_cf_log_t2} is a function that calculates the
#' model for one sector assumption for a given set of data and a
#' shock. This function gets sample data in a matrix form
#' ( either arbitrary matrix or real world data)
#' , shocks, and parameters and returns changes of welfare C based on V(I-M)^-1
#' and summation of I+M+M2+..... It also shows all degree effects separately.

#' 1. C_hat_inv returns results based on V(I-M)^-1.
#' 2. C_hat_sum returns result based on V(I+M+M^2+M^3+.....) .
#' At the end, C_hat_inv and C_hat_sum should convert to a same number with
#'  an infinitesimal difference.
#'
#'
#' The question is: Where to stop the summation I+M+M^2...?
#' We set a while loop to make sure we have calculated enough terms and the final answers converge.
#' While loop criteria is to minimize sum(square(V(I-M)^-1)-(I+M2+M3+.....)(i,j))
#' It can be computed based on other criteria as well if we want.


#' @details This script calculates welfare changes as a result of shock to one sector model:                                   #
#' Here we are using matrix representation of the mode.

#' using options() function we set digits=22 to make sure about the accuracy of the model
#' of  convergence.

#' The objective of using option function is to make sure
#' about the convergence of summation
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
#' @return after computing all necessary variables, if the spectral radius of M
#' is below one, The function returns,V(I-M), I+M+M2.....,the difference,
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
#' pi=as.vector , is a vector of size 2*J*J
#' needed to be divided into intermediate and final pi, pi_I:J*J for intermediate
#' pi_F:J*J for final
#' pi_I=matrix(pi[1:J*J],nrow=J,ncol=J)
#' pi_F=matrix(pi[(J*J)+1:J*J+J*J],nrow=J,ncol=J)
#' gamma=as.vector of size J



#'2. parameters <-
#'   epsilon



#' 3. shock <-
#'    list(
#'      T_hat = matrix(vector(c(), nrow = J,ncol=1)
#'      tau_hat=matrix(c(),nrow=1,ncol=2*J*J)) should be divided into two types of shocks: intermediates and final
#'      tau_hat_I=matrix(tau_hat[1:(J*J)],nrow=J,ncol=J)
#'      tau_hat_F=matrix(tau_hat[(J*J)+1:(J*J)+(J*J)],nrow=J,ncol=J)


#' convergency level 10^-22,10^-20,..... can be also given as an input
#'  (for compatibility to other functions, I exclude this)
#' so the user can choose at what number of digits he/she wants
#'  the convergence happens. Here, by default I chose 10^(-20).
#'It is an arbitrary level and it might be ok for someone else to
#'accept convergency at 10^-43,....etc.
#'
#'@criteria The criteria of convergence is sum square of difference matrix which is sigma(i,j)(V(I-M)^(-1)-(I+M+M2+M3+...)(i,j)
#'We can change the criteria based on need in the while loop. The criteria
#'can also be given as input.(I have exclude for compatibility reason) .
#'In another example I used abs(V(I-M)^-1 - (I+M+M2+....)) and the accuracy is
#' worse and happens at 7th digit. But still converge
#'
#' Comparison of results of functions calc_cf_log_m_2 and calc_cf_log_t2
#' ****
#' The results of this function has been compared to the results of
#' function in task one. With the given example, they give the same result.
#' calc_cf_log_m_2(simple_data,shocks,epsilon)
#' C_hat_1               C_hat_2
#' 1 0.0552016985138004401 0.0552016985138004401
#' 2 0.0076433121019108185 0.0076433121019108159



calc_cf_logl_decomp <- function(J, R, pi_I, pi_F, gamma, T_hat, tau_hat_I,
                                tau_hat_F, epsilon, ktol = 1e-10,
                                kmax_iter = 500) {

  # step 1: Computation of necessary matrices ##################################
  diaggamma <- diag(gamma)
  IJ <- diag(J)
  Igammadiag <- IJ - diaggamma

  pi_Igamma <- crossprod(pi_I, Igammadiag)
  pi_Fgamma <- crossprod(pi_F, Igammadiag)
  tpi_I <- t(pi_I)
  tpi_F <- t(pi_F)

  Rrecip <- 1/R
  RR <- crossprod(Rrecip, R)

  A <- RR * (pi_I %*% Igammadiag)
  Ap <-  RR * (pi_F %*% diaggamma)
  B <- epsilon * A
  Bp <- epsilon * Ap

  Q <- matrix(rep(R * gamma, each = J), ncol = J)

  # step3: Creating Block Matrices: computing matrix M #########################
  M1 <- crossprod(pi_I, Igammadiag)
  M2 <- crossprod(pi_I, diaggamma)
  M3 <- B %*% M1 + Bp %*% pi_Fgamma - epsilon * Igammadiag
  M4 <- A + Ap + B %*% tpi_I %*% diaggamma + Bp %*% tpi_F %*% diaggamma -
    epsilon * diaggamma - Q

  M <- matrix(rbind(cbind(M1, M2), cbind(M3, M4)), ncol = 4)

  # step 4: test spectral radius ###############################################
  spectral_R_M  <- max(abs(eigen(M)$value))

  if(spectral_R_M >= 1) {
    stop("The decomposition does not work since the spectral radius of M is ",
         "greater than 1")
  } else {
    # step 5: computing V matrix ###############################################
    vv1 <- -1 / epsilon * tpi_I
    vv2 <- IJ- A %*% tpi_I - Ap%*%tpi_F
    V1 <- rbind(vv1,vv2)

    vvv1 <- matrix(diag(tpi_I %*% tau_hat_I), nrow = 2)
    vvv2 <- matrix(B %*% diag(tpi_I %*% tau_hat_I) +
                     Bp %*% diag(tpi_F %*% tau_hat_F) -
                     diag(B %*% t(tau_hat_I)) -
                     diag(Bp %*% t(tau_hat_F)),
                   nrow=2)
    V2 <- rbind(vvv1, vvv2)
    V <- V1 %*% T_hat + V2

    # step 6: Computing C_hat based on V(I-M)^(-1) #############################
    z <- solve(diag(J * J) - M)
    pifg <- -pi_Fgamma
    pifI <- IJ - tpi_F %*% diaggamma
    pigamma <- cbind(pifg, pifI)

    S1 <- -diag(crossprod(pi_F, tau_hat_F))
    S2 <- 1 / epsilon * tpi_F %*% T_hat
    C_hat_inv <- S1 + S2 + pigamma %*% z %*% V

    # step 7: Computing C_hat = based on I + M2 + M3 +... ######################
    # loop to make sure convergence happens. Converg2 is a big number=1
    # in loop: while converg2>10^(-20) it continues to compute M^i and stack all
    # in a list and add all of them together to get (I+M+M2+....) while the
    # difference of sumsquare(V(I-M)^(-1)-(I+M+M2+....))=converg2
    converg2 <- 1
    i <- 0
    m5 <- list()
    conv2msg <- list()
    nthmsg2 <- list()
    nth2effect <- list()

    while(converg2 > ktol & i < kmax_iter) {

      m5[i + 1] <- list(M %^% i)

      a5 <- apply(simplify2array(m5), c(1,2), sum)
      converg2 <- colSums(matrix(colSums(apply(z - a5, c(1,2), '^', 2))))

      tmp <- pigamma %*%
        matrix(unlist(m5[i + 1]), ncol = J * J, nrow = J*J) %*% V
      if (i == 0) {
        tmp <- S1 + S2 + tmp
      }
      nth2effect[i + 1] <- list(tmp)

      nthmsg2[i + 1] <- list(
        sprintf("The effect of degree %s is %s", i, nth2effect[i + 1])
      )

      i <- i + 1
    }

    ################ c_hat based on summation is :
    C_hat_sum2 <- S1 + S2 + pigamma %*% a5 %*% V

    #  step 8: final messages to be returned and printed #######################

    finalmsg2 <-
      sprintf(paste0("The difference of summation and V(I-M)^-1 using the ",
                     "criteria of convergence is %s and the convergence ",
                     "happens at %s repetition"), converg2, length(m5))
    radius <-
      sprintf(paste0("The decomposition works since the spectral radius is ",
                     ",%s, which is smaller than 1"), spectral_R_M)

    C_hat_inv_msg <- sprintf("dlogC based on V(I_M)^-1 is %s", list(C_hat_inv))

    zmsg <- sprintf("V(I_M)^-1 is")
    criteriamsg <-
      print(paste0("The criteria of convergence is sigma(i,j)(V(I-M)^(-1) - ",
                   "(I+M+M2+M3+...)(i,j) which is the sum of the square of ",
                   "all elements of the difference matrix"))

    a5msg <- sprintf("The sum of all degrees is")
    C_hat_sum_msg2 <-
      sprintf(paste0("dlogC based on the sum of all degree effects using the ",
                     "criteria of convergence is %s"), list(C_hat_sum2))
  }
  list(radius,
       zmsg,
       matrix(z, ncol = J*J, nrow = J*J),
       a5msg,
       matrix(a5, ncol = J*J, nrow = J*J),
       C_hat_inv_msg,
       C_hat_sum_msg2,
       criteriamsg,
       finalmsg2,
       nthmsg2)
}
