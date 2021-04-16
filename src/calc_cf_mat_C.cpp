#include <iostream>
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// arguments only passed as pointers due to R's internals
// [[Rcpp::export]]
Rcpp::List calc_cf_mat_C(const Rcpp::List data,
                         const Rcpp::List shock,
                         const Rcpp::List parameters,
                         const double tolerance,
                         double zeta,
                         const int maxiter,
                         const int nthreads) {

  // get pointers to actual vectors in the list
  const Rcpp::NumericMatrix R         = data["R"];
  const Rcpp::NumericVector D         = data["D"];
  const Rcpp::NumericVector pi        = data["pi"];
  const Rcpp::NumericMatrix alpha     = data["alpha"];
  const Rcpp::NumericVector gamma_jrs = data["gamma_jrs"];
  const Rcpp::NumericMatrix gamma_js  = data["gamma_js"];

  const Rcpp::NumericVector epsilon   = parameters["epsilon"];
  const double              varphi    = parameters["varphi"];
  const Rcpp::IntegerVector mobility  = parameters["mobility"];

  // get pointers to the shock we are looking at.
  const Rcpp::NumericMatrix T_hat        = shock["T_hat"];
  const Rcpp::NumericVector tau_hat      = shock["tau_hat"];
  const Rcpp::NumericMatrix delta_hat    = shock["delta_hat"];
  const Rcpp::NumericVector varkappa_hat = shock["varkappa_hat"];

  // get number of countries and sectors
  int J = R.nrow();
  int S = R.ncol();
  int U = S + 1;

  // construct initial guesses
  Rcpp::NumericVector P_hat(J * S * U, 1.0);
  Rcpp::NumericMatrix R_hat(J, S);
  R_hat.fill(1.0);

  // initialize helpers
  double norm_sf = 0.0;

  Rcpp::NumericMatrix R_prime(Rcpp::clone(R));
  Rcpp::NumericVector R_prime_new(J * S - 1);
  Rcpp::NumericMatrix R_div(J, S);

  Rcpp::NumericVector D_prime(Rcpp::clone(D));
  D_prime = D * varkappa_hat;

  Rcpp::NumericVector Y_prime(J);
  Rcpp::NumericVector w_hat(J);
  w_hat.fill(1.0);
  Rcpp::NumericVector Y(J);
  for (int j = 0; j < J; j++) {
    for (int s = 0; s < S; s++) {
      Y[j] += gamma_js(j, s) * R(j, s);
    }
  }

  double Y_world_prime = 0.0;
  double Y_world       = 0.0;
  Y_world = Rcpp::sum(Y);


  Rcpp::NumericMatrix L_hat(J, S);
  L_hat.fill(1.0);
  Rcpp::NumericMatrix L_tilde_hat(J, S);
  L_tilde_hat.fill(1.0);
  for (int j = 0; j < J; j++) {
    if (mobility[j] == 0) {
      for (int s = 0; s < S; s++) {
        L_tilde_hat(j, s) = delta_hat(j, s);
      }
    }
  }

  Rcpp::NumericMatrix w_js_hat(J, S);
  w_js_hat.fill(1.0);

  Rcpp::NumericMatrix c_hat(J, S);
  c_hat.fill(1.0);
  Rcpp::NumericVector P_hat_new(J * S * U, 1.0);
  Rcpp::NumericVector P_hat_div(J * S * U);

  Rcpp::NumericVector multres(J * J * S * U);
  Rcpp::NumericVector pi_hat(J * J * S * U, 1.0);
  Rcpp::NumericVector multres_tmp(J * S * U);

  // Armadillo Matrix and Vector for solving the market clearing linear system
  // addtional row and value for normalization
  arma::mat E_part(J * S + 1, J * S, arma::fill::ones);
  arma::vec b(J * S + 1, arma::fill::zeros);

  // normalize to world GDP
  for (int s = 0; s < S; s++) {
    for (int j = 0; j < J; j++) {
      E_part(J * S, j + s * J) = gamma_js(j, s);
    }
  }


  Rcpp::NumericVector tmp(J * S * U);

  Rcpp::Rcout << "Data initialized. Starting iterative procedure." << std::endl;

  // iterative solving of equation system
  double div  = 1 + tolerance;
  double div_new = 0;
  double div_old = 0;
  int    iter = 0;
  int    iter_converged = 0;

  while(tolerance <= div & iter < maxiter)  {

    iter++;

    // w_hat_js
    std::fill(Y_prime.begin(), Y_prime.end(), 0.0);
    for (int j = 0; j < J; j++) {
      for (int s = 0; s < S; s++) {
        Y_prime[j] += gamma_js(j, s) * R_prime(j, s);
      }
    }


    for (int j = 0; j < J; j++) {
      w_hat[j] = Y_prime[j] / Y[j];
    }

    for (int j = 0; j < J; j++) {
      if (mobility[j] != 0) {
        for (int s = 0; s < S; s++) {
          // per capita/ex-ante expected wage equalizes with mobility
          L_hat(j, s) = R_hat(j, s) / w_hat[j];
        }
      }
      if (mobility[j] == 1) {
        for (int s = 0; s < S; s++) {
          L_tilde_hat(j, s) = delta_hat(j, s) * pow(L_hat(j, s), 1 - 1/varphi);
        }
      } else if (mobility[j] == 2) {
        for (int s = 0; s < S; s++) {
          L_tilde_hat(j, s) = delta_hat(j, s) * L_hat(j, s);
        }
      }
    }

    for (int j = 0; j < J; j++) {
      for (int s = 0; s < S; s++) {
        w_js_hat(j, s) = R_hat(j, s) / L_tilde_hat(j, s);
      }
    }

    for (int j = 0; j < J; j++) {
      for (int s = 0; s < S; s++) {
        c_hat(j, s) = pow(w_js_hat(j, s), gamma_js(j, s));
      }
    }
    for (int j = 0; j < J; j++) {
      for (int s = 0; s < S; s++) {
        for (int r = 0; r < S; r++) {
          c_hat(j, s) *= pow(P_hat[j + r * J + s * J * S],
                gamma_jrs[j + r * J + s * J * S]);
        }
      }
    }

    // pi_hat and P_hat_new
    std::fill(multres_tmp.begin(), multres_tmp.end(), 0.0);
#pragma omp parallel for num_threads(nthreads)
    for (int u = 0; u < S + 1; u++) {
      for (int r = 0; r < S; r++) {
        for (int j = 0; j < J; j++) {
          for (int i = 0; i < J; i++) {
            int idx_ijru = i + j * J + r * J * J + u * S * J * J;
            //  if(Rcpp::traits::is_infinite<REALSXP>(tau_hat[idx_ijru])) {
            //    multres[idx_ijru] = 0;
            //    Rcpp::Rcout << "in: " << i << ", " << j << ", " << r << ", " << u << ";  ";
            //  } else {
            multres[idx_ijru] = T_hat(i, r) *
              pow(c_hat(i, r) * tau_hat[idx_ijru], -epsilon[r]);
            //  }
          }
        }
      }
    }

#pragma omp parallel for num_threads(nthreads)
    for (int u = 0; u < S + 1; u++) {
      for (int r = 0; r < S; r++) {
        for (int j = 0; j < J; j++) {
          int idx_jru = j + r * J + u * J * S;
          for (int i = 0; i < J; i++) {
            int idx_ijru = i + j * J + r * J * J + u * S * J * J;
            multres_tmp[idx_jru] += pi[idx_ijru] * multres[idx_ijru];
          }
        }
      }
    }

    for (int idx = 0; idx < J * S * (S + 1); idx++) {
      if(multres_tmp[idx] != 0) {
        for (int i = 0; i < J; i++) {
          pi_hat[i + idx * J] =  multres[i + idx * J] / multres_tmp[idx];
        }
      } else {
        for (int i = 0; i < J; i++) {
          pi_hat[i + idx * J] = 1;
        }
      }
    }

#pragma omp parallel for num_threads(nthreads)
    for (int u = 0; u < S + 1; u++) {
      for (int r = 0; r < S; r++) {
        for (int j = 0; j < J; j++) {
          int idx_jru = j + r * J + u * J * S;
          if(multres_tmp[idx_jru] != 0) {
            P_hat_new[idx_jru] = pow(multres_tmp[idx_jru], - 1 / epsilon[r]);
          } else {
            P_hat_new[idx_jru] = 1;
          }
        }
      }
    }

    // ensure that not only w but also P has converged
    P_hat_div = P_hat_new - P_hat;
    P_hat = P_hat_new;


    // Linear equation system of market clearing to solve for R_prime
    for (int j = 0; j < J; j++) {
      for (int r = 0; r < S; r++) {
        for (int u = 0; u < S; u++) {
          int idx_jru = j + r * J + u * J * S;
          for (int i = 0; i < J; i++) {
            int idx_ijru = i + j * J + r * J * J + u * J * J * S;
            int idx_ijrU = i + j * J + r * J * J + S * J * J * S;
            E_part(i + J * r, j + J * u) =
              pi_hat[idx_ijru] * pi[idx_ijru] * gamma_jrs[idx_jru] +
              pi_hat[idx_ijrU] * pi[idx_ijrU] * alpha(j, r) * gamma_js(j, u);
            if(i == j & r == u) {
              E_part(i + J * r, j + J * u) -= 1;
            }
          }
        }
      }
    }

    b.zeros();
    b(J * S) = Y_world;
    for (int i = 0; i < J; i++) {
      for (int r = 0; r < S; r++) {
        for (int j = 0; j < J; j++) {
          int idx_ijrU = i + j * J + r * J * J + S * J * J * S;
          b(i + J * r) -= pi_hat[idx_ijrU] * pi[idx_ijrU] * alpha(j, r) * D_prime[j];
        }
      }
    }

    R_prime_new = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(arma::solve(E_part, b)));
    R_prime_new.attr("dim") = Rcpp::Dimension(J, S);

    // In select cases we might be up to a rounding error below 0, set this to
    // R's standard precision ~ 1e-16
    for (int i = 0; i < J; i++) {
      for (int r = 0; r < S; r++) {
        if(R(i, r) == 0) {
          R_prime_new(i, r) = 0;
        }
        if(R_prime_new(i, r) < 0) {
          R_prime_new(i, r) = 1.0e-16;
        }
      }
    }

    // scale to world GDP
    Y_world_prime = 0;
    for (int i = 0; i < J; i++) {
      for (int r = 0; r < S; r++) {
        Y_world_prime += gamma_js(i, r) * R_prime_new(i, r);
      }
    }

    for (int j = 0; j < J; j++) {
      for (int s = 0; s < S; s++) {
        R_div(j, s) = R_prime_new(j, s) - R_prime(j, s);
        R_prime(j, s) = R_prime(j, s) +
          zeta * (R_prime_new(j, s) - R_prime(j, s));
        if(R(j, s) != 0) {
          R_hat(j, s) = R_prime(j, s) / R(j, s);
        } else {
          R_hat(j,s) = 1;
        }
      }
    }



    // find largest change in R_hat or P_hat in this iteration
    double divmax = std::max(Rcpp::max(R_div), Rcpp::max(P_hat_div));
    double divmin = std::min(Rcpp::min(R_div), Rcpp::min(P_hat_div));

    // absolute
    if (divmax < -divmin) {
      div = -divmin;
    } else {
      div = divmax;
    }


    Rcpp::checkUserInterrupt();
    Rcpp::Rcout << "Iteration: " << iter <<"  The current div is:  " << div <<
      "\r" << std::flush;

  }

  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "Iterating finished. Caculating additional statistics." << std::endl;

  Rcpp::NumericVector E_j_prime(J);
  Rcpp::NumericVector E_j(J);
  Rcpp::NumericVector C_hat(J, 1.0);


  E_j = Y + D;
  E_j_prime = Y_prime + D_prime;
  for (int j = 0; j < J; j++) {
    for (int r = 0; r < S; r++) {
      C_hat[j] *=  1 / pow(P_hat(j + r * J + J * S * S), alpha(j, r));
    }
    C_hat[j] *= E_j_prime[j] / E_j[j];
  }

  Rcpp::NumericVector M_prime(J * S * S);
  M_prime.fill(0.0);
  for (int j = 0; j < J; j++) {
    for (int s = 0; s < S; s++) {
      if (R_prime(j, s) != 0) {
        for (int r = 0; r < S; r++) {
          M_prime[j + r * J + s * J * S] =
            gamma_jrs[j + r * J + s * J * S] * R_prime(j, s);
        }
      }
    }
  }

  Rcpp::NumericMatrix E_prime(J, S);
  for(int j = 0; j < J; j++) {
    for (int s = 0; s < S; s++) {
      E_prime(j, s) = alpha(j, s) * E_j_prime[j];
    }
  }

  Rcpp::NumericMatrix Q_hat(J, S);
  Rcpp::NumericMatrix Q_tmp(J, S);
  for (int i = 0; i < J; i++) {
    for (int r = 0; r < S; r++) {
      for (int j = 0; j < J; j++) {
        for (int u = 0; u < S; u++) {
          int idx_ijru = i + j * J + r * J * J + u * J * J * S;
          Q_hat(i, r) += pow(pi_hat[idx_ijru] * pi[idx_ijru],
                1 - 1 / epsilon[r])  * M_prime[j + r * J + u * J * S];
          Q_tmp(i, r) += pow(pi[idx_ijru], 1 - 1 / epsilon[r])  *
            gamma_jrs[j + r * J + u * J * S] * R_prime(j, u);
        }
        int idx_ijrS = i + j * J + r * J * J + S * J * J * S;
        Q_hat(i, r) += pow(pi_hat[idx_ijrS] * pi[idx_ijrS],
              1 - 1 / epsilon[r]) * E_prime(j, r);
        Q_tmp(i, r) += pow(pi[idx_ijrS], 1 - 1 / epsilon[r]) * alpha(j, r) *
          E_j[j];
      }
      if (Q_tmp(i, r) != 0) {
        Q_hat(i, r) /= Q_tmp(i, r);
        Q_hat(i, r) *= pow(T_hat(i, r), 1 /   epsilon[r]) / c_hat(i, r);
      } else {
        Q_hat(i, r) = 1.0;
      }
    }
  }


  return Rcpp::List::create(Rcpp::Named("R_hat") = R_hat,
                            Rcpp::Named("P_hat") = P_hat,
                            Rcpp::Named("multres") = multres,
                            Rcpp::Named("multres_tmp") = multres_tmp,
                            Rcpp::Named("R_div") = R_div,
                            Rcpp::Named("P_hat_div") = P_hat_div,
                            Rcpp::Named("C_hat") = C_hat,
                            Rcpp::Named("Q_hat") = Q_hat,
                            Rcpp::Named("D_prime") = D_prime,
                            Rcpp::Named("Y_prime") = Y_prime,
                            Rcpp::Named("L_hat") = L_hat,
                            Rcpp::Named("L_tilde_hat") = L_tilde_hat,
                            Rcpp::Named("c_hat") = c_hat,
                            Rcpp::Named("pi_hat") = pi_hat);
}
