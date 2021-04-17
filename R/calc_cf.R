#' Calculate/Simulate model counterfactual
#'
#' TODO: CLEAN THIS UP
#' \code{calc_cf} is a wrapper function for the cpp function that calculates the
#' model counterfactual described in the paper for a given set of data and a
#' shock.
#'
#' @details The main purpose of this function is twofold. Firstly to do some
#' basic checks and input manipulations (e.g. from character to numeric) before
#' calling the cpp code. Secondly, the main set of market clearing equations of
#' the model can be solved as a linear equation system at each step of the
#' iterative counterfactual solution algorithm or one can directly use right
#' hand side values of R. The former converges much faster and is more stable
#' in converging but problematic when the simulation involves shocks with
#' infinitely high trade barriers (preventing a matrix inversion to sole the
#' linear system). Therefore this function also decides which cpp algorithm to
#' pick based upon the shock.
#'
#' @param data a list of model variables possibly derived from \code{convert_io}
#'             which contains
#'   \describe{
#'   \item{location_id}{vector of location ids}
#'   \item{sector_id}{vector of sector ids}
#'   \item{R}{matrix of location-sector revenues with
#'            \code{N = length(location_id)} rows and
#'            \code{J = length(sector_id) columns}}
#'   \item{D}{vector of locations' deficit transfers, i.e. positive for a trade
#'            deficit and negative for a trade surplus}
#'   \item{pi}{vector of import shares across locations for each
#'             sector-destination-use combination, where the import share of
#'             origin \code{o} in sector \code{s} products used in destination
#'             \code{d} and use category \code{u} is at position
#'             \code{o + (d-1)*N + (s-1)*N*N + (u-1)*N*N*J}}
#'   \item{alpha}{matrix of sectoral consumption shares in each location with
#'                \code{N} rows and \code{J} columns}
#'   \item{gamma_jrs}{vector of intermediate cost shares across sectrs for each
#'                    destination-use combination, where the cost share of
#'                    sector \code{s} intermediates used in destination
#'                    \code{d} and use category \code{u} is at position
#'                    \code{d + (s-1)*N + (u-1)*N*J}}
#'   \item{gamma_js}{matrix of location-sector labor cost shares with
#'                \code{N} rows and \code{J} columns}}
#' @param shock a list containing a shock with values
#'   \describe{
#'     \item{T_hat}{technological shock}}
#' @param tol tolerance for algorithm to stop
#' @param zeta dampening factor
#' @param maxiter maximum number of iterations
#' @param method either matrix or RHS
#' @param nthreads number of threads to use in parallel. This relies on OpenMP
#' and thus has no effect on MacOS.
#' @return Returns a list of R_hat, P_hat, multres, multres_tmp, R_div,
#' P_hat_div, C_hat, Q_hat, D_prime, L_prime, L_tilde_hat, c_hat, pi_hat
#' @importFrom magrittr %>%
#' @export calc_cf

calc_cf <- function(data,
                    shock,
                    parameters,
                    tol = 1e-8,
                    zeta = 0.3,
                    maxiter = 10000,
                    method = "automatic",
                    nthreads = 1L) {

  if(method %in% c("linearized_decomp", "linearized", "pseudo_linearized")) {
    if(dim(data[["R"]])[2] > 1) {
      stop("Log linearized model works with 1 sector only")
    }
    # extract from list
    R <- data["R"]
    pi <- data[["pi"]]
    gamma <- as.vector(simple_data[["gamma_js"]])

    # and reformat
    J <- dim(R)[1]
    pi_I <- matrix(pi[1:J^2], nrow = J, ncol = J)
    pi_F <- matrix(pi[(J^2 + 1):(J^2 + J^2)], nrow = J, ncol = J)
    T_hat <- shocks[["T_hat"]]
    tau_hat <- shocks[["tau_hat"]]
    tau_hat_I <- matrix(tau_hat[1:(J^2)], nrow = J, ncol = J)
    tau_hat_F <- matrix(tau_hat[(J^2 + 1):(J^2 + J^2)], nrow = J, ncol = J)

    if(method %in% c("linearized", "pseudo_linearized")) {
      if(method == "pseudo_linearized") {
        pseudo = TRUE
      } else {
        pseudo = FALSE
      }
      calc_cf_logl(J, R, pi_I, pi_F, gamma, T_hat, tau_hat_I, tau_hat_F,
                   parameters$epsilon, use_pseudo = pseudo)
    } else {
      calc_cf_logl_decomp(J, R, pi_I, pi_F, gamma, T_hat, tau_hat_I,
                   tau_hat_F, epsilon, ktol = tol, kmax_iter = maxiter)
    }
  }
  if (!all(parameters$mobility %in% c("perfect", "imperfect", "immobile"))) {
    print(paste0("Invalid mobility argument (allowed are only \"perfect\",)",
                 "\"imperfect\", \"immobile\"). Assuming perfect mobility."))
  }
  parameters$mobility <- case_when(parameters$mobility == "immobile" ~ 0,
                                   parameters$mobility == "imperfect" ~ 1,
                                   parameters$mobility == "perfect" ~ 2,
                                   TRUE ~ 2)

  if(length(parameters$mobility) == 1) {
    parameters$mobility <- rep(parameters$mobility, length(data$location_id))
  } else if (length(parameters$mobility) != length(data$location_id)) {
    stop(paste0("parameters$mobility must be either a single value or one",
                "value for each location"))
  }

  method_auto = "matrix_inv"
  if (method == "automatic") {
    print("Checking for autarky situations in counterfactual")

    J <- length(data$location_id)
    S <- length(data$sector_id)
    tau_hat <- array(shock$tau_hat, dim = c(J, J, S, S + 1))
    j <- 1
    while(j < J & method_auto == "matrix_inv") {
      if(all(is.infinite(tau_hat[-j,j,,])) |
         all(array(data$pi, dim = c(J, J, S, S+1))[-j,j,,] == 0)) {
        method_auto <- "direct"
      }
      j = j + 1
    }
  }

  if(method_auto == "direct" | method == "direct") {
    calc_cf_C(data, shock, parameters, tol, zeta, as.integer(maxiter), nthreads)
  } else {
    calc_cf_mat_C(data, shock, parameters, tol, zeta, as.integer(maxiter), nthreads)
  }
}
