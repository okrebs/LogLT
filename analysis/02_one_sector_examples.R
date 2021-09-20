library(dplyr)
library(tidyr)
library(haven) # for stata output
library(LogLT)

#simple mock data example ######################################################

n_location <- 2
n_sector <- 1

# manually create some data
simple_data <-
  list(
    location_id = seq(n_location),
    sector_id = seq(n_sector),
    R = matrix(c(1, 1.5), nrow = n_location, ncol = n_sector),
    D = rep(0, n_location),
    pi =
      as.vector(array(c(0.6, 0.4, 0.2, 0.8, 0.5, 0.5, 0.4, 0.6),
                      dim = c(n_location, n_location, n_sector, n_sector + 1))),
    alpha = matrix(1, nrow = n_location, ncol = n_sector),
    gamma_jrs = as.vector(array(c(0.5, 0.5),
                                dim = c(n_location, n_sector, n_sector))),
    gamma_js = matrix(c(0.5, 0.5), nrow = n_location, ncol = n_sector)
  )

parameters <-
  list(
    epsilon = 3,
    varphi = 1.5, # does not matter with immobility
    mobility = "immobile"
  )

shock <-
  list(
    T_hat = matrix(c(1.1, 1), nrow = n_location, ncol = n_sector),
    tau_hat = rep(1, n_location * n_location * n_sector * (n_sector + 1)),
    # labor productivity changes:
    delta_hat = matrix(1, nrow = n_location, ncol = n_sector),
    # deficit changes
    varkappa_hat = rep(1, n_location)
  )

# non-linear model
test_ge <- calc_cf(simple_data, shock, parameters)$C_hat - 1
# log linear model, old code
test_ll <- calc_cf(simple_data, shock, parameters, method = "linearized")
# log linear model, decomposition, new code
test_decomp <- calc_cf(simple_data, shock, parameters,
                       method = "linearized_decomp")

# wiot one sector example ######################################################

baseline_data <-
  readRDS("analysis/data/one_sector_baseline_data.rds")

locations <- baseline_data$location_id
J <- n_distinct(baseline_data$location_id)
S <- n_distinct(baseline_data$sector_id)

baseline_shock <-
  list(
    "T_hat" = matrix(1, nrow = J, ncol = S),
    "tau_hat" = rep(1, J * J * (S + 1)),
    "delta_hat" = matrix(1, nrow = J, ncol = S),
    "varkappa_hat" = rep(1, J)
  )

parameters <- list("epsilon" = 5,
                   "varphi" = 1.5,
                   "mobility" = "perfect")

# Run simulations ==============================================================

# Technology shocks ------------------------------------------------------------

welfare_ge <- list()
welfare_ll <- list()
welfare_technology_shock <- list()
shock_size <- c(1, 5, 10, 50, 100)

for (s in 1:length(shock_size)) {
  for (i in 1:J) {
    shock <- baseline_shock
    shock$T_hat[i, ] <- 1 + shock_size[s]/100

    welfare_ge[[i]] <- calc_cf(baseline_data,
                               shock, parameters, tol = 1e-7)$C_hat - 1

    welfare_ll[[i]] <- calc_cf(baseline_data,
                               shock,
                               parameters,
                               method = "linearized_decomp",
                               maxiter = 50)
  }

  welfare_technology_shock[[s]] <-
    expand_grid(shocked_country = locations,
                country = locations) %>%
    mutate(
      shock_size_pct = shock_size[s],
      ge_nonlinear_C_hat = unlist(welfare_ge)
    ) %>%
    bind_cols(bind_rows(welfare_ll))
}

welfare_technology_shock <- bind_rows(welfare_technology_shock)

# output to stata
write_dta(welfare_technology_shock,
          "analysis/data/welfare_technology_shock.dta")

rm(welfare_technology_shock)

# Do the same with trade barrier shocks instead --------------------------------

welfare_ge <- list()
welfare_ll <- list()
welfare_tradebarrier_shock <- list()
shock_size <- c(1, 5, 10, 50, 100)

for (s in 1:length(shock_size)) {
  for (i in 1:J) {
    shock <- baseline_shock
    # barrier shock for all imports of i, except for goods produced in i
    tmp <- array(shock$tau_hat, dim = c(J, J, S + 1))
    tmp[-i, i,] <- 1 + shock_size[s]/100
    shock$tau_hat <- as.vector(tmp)


    welfare_ge[[i]] <- calc_cf(baseline_data,
                               shock, parameters, tol = 1e-7)$C_hat - 1

    welfare_ll[[i]] <- calc_cf(baseline_data,
                               shock,
                               parameters,
                               method = "linearized_decomp",
                               maxiter = 50)
  }

  welfare_tradebarrier_shock[[s]] <-
    expand_grid(shocked_country = locations,
                country = locations) %>%
    mutate(
      shock_size_pct = shock_size[s],
      ge_nonlinear_C_hat = unlist(welfare_ge)
    ) %>%
    bind_cols(bind_rows(welfare_ll))
}

welfare_tradebarrier_shock <- bind_rows(welfare_tradebarrier_shock)

# output to stata
write_dta(welfare_tradebarrier_shock,
          "analysis/data/welfare_tradebarrier_shock.dta")

rm(welfare_tradebarrier_shock)

# Do the same with separate final/intermediate trade barrier shocks instead ----

# Intermediate barriers
welfare_ge <- list()
welfare_ll <- list()
welfare_intermediate_shock <- list()
shock_size <- c(1, 5, 10, 50, 100)

for (s in 1:length(shock_size)) {
  for (i in 1:J) {
    shock <- baseline_shock
    # barrier shock for all imports of i, except for goods produced in i
    tmp <- array(shock$tau_hat, dim = c(J, J, S + 1))
    tmp[-i, i, 1] <- 1 + shock_size[s]/100
    shock$tau_hat <- as.vector(tmp)


    welfare_ge[[i]] <- calc_cf(baseline_data,
                               shock, parameters, tol = 1e-7)$C_hat - 1

    welfare_ll[[i]] <- calc_cf(baseline_data,
                               shock,
                               parameters,
                               method = "linearized_decomp",
                               maxiter = 50)
  }

  welfare_intermediate_shock[[s]] <-
    expand_grid(shocked_country = locations,
                country = locations) %>%
    mutate(
      shock_size_pct = shock_size[s],
      ge_nonlinear_C_hat = unlist(welfare_ge)
    ) %>%
    bind_cols(bind_rows(welfare_ll))
}

welfare_intermediate_shock <- bind_rows(welfare_intermediate_shock)

# output to stata
write_dta(welfare_intermediate_shock,
          "analysis/data/welfare_intermediate_shock.dta")

rm(welfare_intermediate_shock)

# Final barriers
welfare_ge <- list()
welfare_ll <- list()
welfare_final_shock <- list()
shock_size <- c(1, 5, 10, 50, 100)

for (s in 1:length(shock_size)) {
  for (i in 1:J) {
    shock <- baseline_shock
    # barrier shock for all imports of i, except for goods produced in i
    tmp <- array(shock$tau_hat, dim = c(J, J, S + 1))
    tmp[-i, i, 2] <- 1 + shock_size[s]/100
    shock$tau_hat <- as.vector(tmp)


    welfare_ge[[i]] <- calc_cf(baseline_data,
                               shock, parameters, tol = 1e-7)$C_hat - 1

    welfare_ll[[i]] <- calc_cf(baseline_data,
                               shock,
                               parameters,
                               method = "linearized_decomp",
                               maxiter = 50)
  }

  welfare_final_shock[[s]] <-
    expand_grid(shocked_country = locations,
                country = locations) %>%
    mutate(
      shock_size_pct = shock_size[s],
      ge_nonlinear_C_hat = unlist(welfare_ge)
    ) %>%
    bind_cols(bind_rows(welfare_ll))
}

welfare_final_shock <- bind_rows(welfare_final_shock)

# output to stata
write_dta(welfare_final_shock,
          "analysis/data/welfare_final_shock.dta")

rm(welfare_final_shock)

# Individual trade barrier pair shocks -----------------------------------------

welfare_ge <- vector(mode = "list", length = J) # will be list of lists
welfare_ll <- vector(mode = "list", length = J)
welfare_pair_interm_shocks <- list()
shock_size <- c(1, 5, 10, 50, 100)

for (s in 1:length(shock_size)) {
  for (j in 1:J) {
    for (i in 1:J) {
      shock <- baseline_shock
      # barrier shock for imports of j from i
      tmp <- array(shock$tau_hat, dim = c(J, J, S + 1))
      tmp[i, j, 1] <- 1 + shock_size[s]/100
      shock$tau_hat <- as.vector(tmp)


      welfare_ge[[i]][[j]] <- calc_cf(baseline_data,
                                      shock, parameters, tol = 1e-7)$C_hat - 1

      welfare_ll[[i]][[j]] <- calc_cf(baseline_data,
                                      shock,
                                      parameters,
                                      method = "linearized_decomp",
                                      maxiter = 50)
    }
  }

  welfare_pair_interm_shocks[[s]] <-
    expand_grid(shocked_exporter = locations,
                shocked_importer = locations,
                country = locations) %>%
    mutate(
      shock_size_pct = shock_size[s],
      ge_nonlinear_C_hat = unlist(welfare_ge)
    ) %>%
    bind_cols(bind_rows(welfare_ll))
}

welfare_pair_interm_shocks <- bind_rows(welfare_pair_interm_shocks)

# output to stata
write_dta(welfare_pair_interm_shocks,
          "analysis/data/welfare_pair_interm_shocks.dta")

rm(welfare_pair_interm_shocks)

# same for final goods

welfare_ge <- vector(mode = "list", length = J) # will be list of lists
welfare_ll <- vector(mode = "list", length = J)
welfare_pair_final_shocks <- list()
shock_size <- c(1, 5, 10, 50, 100)

for (s in 1:length(shock_size)) {
  for (j in 1:J) {
    for (i in 1:J) {
      shock <- baseline_shock
      # barrier shock for imports of j from i
      tmp <- array(shock$tau_hat, dim = c(J, J, S + 1))
      tmp[i, j, 2] <- 1 + shock_size[s]/100
      shock$tau_hat <- as.vector(tmp)


      welfare_ge[[i]][[j]] <- calc_cf(baseline_data,
                                      shock, parameters, tol = 1e-7)$C_hat - 1

      welfare_ll[[i]][[j]] <- calc_cf(baseline_data,
                                      shock,
                                      parameters,
                                      method = "linearized_decomp",
                                      maxiter = 50)
    }
  }

  welfare_pair_final_shocks[[s]] <-
    expand_grid(shocked_exporter = locations,
                shocked_importer = locations,
                country = locations) %>%
    mutate(
      shock_size_pct = shock_size[s],
      ge_nonlinear_C_hat = unlist(welfare_ge)
    ) %>%
    bind_cols(bind_rows(welfare_ll))
}

welfare_pair_final_shocks <- bind_rows(welfare_pair_final_shocks)

# output to stata
write_dta(welfare_pair_final_shocks,
          "analysis/data/welfare_pair_final_shocks.dta")

rm(welfare_pair_final_shocks)
