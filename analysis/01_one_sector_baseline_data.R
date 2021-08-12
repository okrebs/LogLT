library(dplyr)
library(tidyr)
# if not yet done install the iotr package
# library(devtools)
# install_github(repo = "okrebs/iotr")
library(iotr) # for WIOD data
library(LogLT)

# get wiot ---------------------------------------------------------------------
wiot <- get_wiot(years = c(2014)#,
                 # set a cache path so you do not have to download this again
                 #cache = "/something/something"
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

saveRDS(baseline_data, "analysis/data/one_sector_baseline_data.rds")
