library(dplyr)
library(tidyr)
library(haven) # for stata output
# if not yet done install the iotr package
# library(devtools)
# install_github(repo = "okrebs/iotr")
library(iotr) # for WIOD data
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
    T_hat = matrix(c(1.1,1), nrow = n_location, ncol = n_sector),
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
# test_decomp[[6]]
# test_decomp[[7]]

# involved example =============================================================

# get wiot ---------------------------------------------------------------------
wiot <- get_wiot(years = c(2014)#,
                 # set a cache path so you do not have to download this again
                 #cache = "/home/uxb/RawData/WIOD/"
                 ) %>%
  wiot2long() %>%
  filter(Country != "TOT") %>%
  select(origin = Country, sector = RNr, destination, use, flow)

wiot <- wiot %>%
  rm_dynamics(dynamic_categories = c(60, 61), category_to_scale = 57)

# after adapting, aggregate ----------------------------------------------------
wiot <- wiot %>%
  mutate(sector = 1,
         use = ifelse(use >= 57, 2L, 1L)) %>%
  group_by(origin, sector, destination, use) %>%
  summarise(flow = sum(flow), .groups = "drop")

wiot <- wiot %>%
  gen_own_trade(max_replace = 1e-6)

model_data <- convert_io(wiot)

# remove trade imbalances ------------------------------------------------------
locations <- unique(wiot$origin)
J <- length(locations)
S <- length(unique(wiot$sector))
shock_no_deficit <-
  list("T_hat" = matrix(1, nrow = J, ncol = S),
       "tau_hat" = rep(1, J * J * (S + 1)),
       "delta_hat" = matrix(1, nrow = J, ncol = S),
       "varkappa_hat" = rep(0, J) # no deficit
  )
parameters <- list("epsilon" = 5,
                   "varphi" = 1.5,
                   "mobility" = "perfect")

no_deficit <- calc_cf(model_data, shock_no_deficit, parameters, tol = 1e-7)

baseline_data <- model_data
baseline_data$R <- baseline_data$R * no_deficit$R_hat
baseline_data$pi <- baseline_data$pi * no_deficit$pi_hat
baseline_data$D <- no_deficit$D_prime

# Run simulations ==============================================================

test_ge <- list()
test_ll <- list()

for(i in 1:length(locations)) {
  shock <-
    list("T_hat" = matrix(1, nrow = J, ncol = S),
         "tau_hat" = rep(1, J * J * (S + 1)),
         "delta_hat" = matrix(1, nrow = J, ncol = S),
         "varkappa_hat" = rep(1, J)
    )
  shock$T_hat[i,] <- 1.1 # 10% shock in country i

  test_ge[[i]] <- calc_cf(baseline_data, shock, parameters, tol = 1e-7)$C_hat - 1
  # for now the decomposition code does not work
  #test_decomp <- calc_cf(baseline_data, shock, parameters,
                          #method = "linearized_decomp")
  # so I use the old code from task 1
  test_ll[[i]] <- calc_cf(baseline_data, shock, parameters, method = "linearized")
}

test_technology_shock <- expand_grid(shocked_country = locations,
                                     country = locations) %>%
  mutate(ge_nonlinear_C_hat = unlist(test_ge),
         ge_loglinear_dlogC = unlist(test_ll))

#correlation is extremly high
#cor(test_technology_shock$ge_nonlinear_C_hat,
#    test_technology_shock$ge_loglinear_dlogC)

# output to stata
write_dta(test_technology_shock, "analysis/data/test_technology_shock.dta")

# Do the same with trade barrier shocks instead --------------------------------
test_ge <- list()
test_ll <- list()

for(i in 1:length(locations)) {
  shock <-
    list("T_hat" = matrix(1, nrow = J, ncol = S),
         "tau_hat" = rep(1, J * J * (S + 1)),
         "delta_hat" = matrix(1, nrow = J, ncol = S),
         "varkappa_hat" = rep(1, J)
    )
  # 10% barrier shock for all imports of i, except for goods produced in i
  tmp <- array(shock$tau_hat, dim = c(J, J, S + 1))
  tmp[-i, i,] <- 1.1
  shock$tau_hat <- as.vector(tmp)

  test_ge[[i]] <- calc_cf(baseline_data, shock, parameters, tol = 1e-7)$C_hat - 1
  # for now the decomposition code does not work
  #test_decomp <- calc_cf(baseline_data, shock, parameters,
  #method = "linearized_decomp")
  # so I use the old code from task 1
  test_ll[[i]] <- calc_cf(baseline_data, shock, parameters, method = "linearized")
}

test_tradebarrier_shock <- expand_grid(shocked_country = locations,
                                       country = locations) %>%
  mutate(ge_nonlinear_C_hat = unlist(test_ge),
         ge_loglinear_dlogC = unlist(test_ll))

#correlation is extremly high
# cor(test_tradebarrier_shock$ge_nonlinear_C_hat,
#    test_tradebarrier_shock$ge_loglinear_dlogC)

# output to stata
write_dta(test_tradebarrier_shock, "analysis/data/test_tradebarrier_shock.dta")
