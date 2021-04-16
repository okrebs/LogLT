#include <iostream>
#include <Rcpp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// arguments onyl passed as pointers due to R's internals
// [[Rcpp::export]]
Rcpp::List calc_cf_C(const Rcpp::List data,
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
  double norm_sf       = 0.0;

  Rcpp::NumericMatrix R_prime(Rcpp::clone(R));
  Rcpp::NumericMatrix R_hat_new(J, S);
  Rcpp::NumericMatrix R_hat_old(J, S);
  Rcpp::NumericMatrix R_div(J, S);

  Rcpp::NumericVector D_prime(Rcpp::clone(D));
  D_prime = D * varkappa_hat;

  Rcpp::NumericVector Y_prime(J);
  Rcpp::NumericVector Y(J);
  for (int j = 0; j < J; j++) {
    for (int s = 0; s < S; s++) {
      Y[j] += gamma_js(j, s) * R(j, s);
    }
  }

  double Y_world_prime = 0.0;
  double Y_world       = 0.0;
  Y_world = Rcpp::sum(Y);

  Rcpp::NumericVector M_prime(J * S * S);
  Rcpp::NumericMatrix E_prime(J, S);

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
  Rcpp::NumericVector P_hat_old(J * S * U, 1.0);
  Rcpp::NumericVector P_hat_div(J * S * U);

  Rcpp::NumericVector multres(J * J * S * U);
  Rcpp::NumericVector pi_hat(J * J * S * U, 1.0);
  Rcpp::NumericVector multres_tmp(J * S * U);

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

    // Y_prime, M_prime, E_prime, c_hat, L_hat, L_tilde_hat, w_js_hat
    std::fill(Y_prime.begin(), Y_prime.end(), 0.0);
    for (int j = 0; j < J; j++) {
      for (int s = 0; s < S; s++) {

        Y_prime[j] +=  gamma_js(j, s) * R_prime(j, s);

        if (R_prime(j, s) != 0) {
          for (int r = 0; r < S; r++) {
            M_prime[j + r * J + s * J * S] =
              gamma_jrs[j + r * J + s * J * S] * R_prime(j, s);
          }
        }
      }

      for (int s = 0; s < S; s++) {

        E_prime(j, s) = alpha(j, s) * (Y_prime[j] + D_prime[j]);

        if (R(j, s) != 0) {
          if (mobility[j] == 1) {
            L_hat(j, s) = R_hat(j, s) / (Y_prime[j] / Y[j]);
            L_tilde_hat(j, s) = delta_hat(j, s) * pow(L_hat(j, s), 1 - 1 / varphi);
          } else if (mobility[j] == 2) {
            L_hat(j, s) = R_hat(j, s) / (Y_prime[j] / Y[j]);
            L_tilde_hat(j, s) = delta_hat(j, s) * L_hat(j, s);
          }
          w_js_hat(j, s) = R_hat(j, s) / L_tilde_hat(j, s);
          c_hat(j, s) = pow(w_js_hat(j, s), gamma_js(j, s));
          for (int r = 0; r < S; r++) {
            c_hat(j, s) *= pow(P_hat[j + r * J + s * J * S],
                  gamma_jrs[j + r * J + s * J * S]);
          }
        }
      }
    }

    // pi_hat and P_hat_new
    std::fill(multres_tmp.begin(), multres_tmp.end(), 0.0);
#pragma omp parallel for num_threads(nthreads)
    for (int u = 0; u < S + 1; u++) {
      for (int r = 0; r < S; r++) {
        for (int j = 0; j < J; j++) {
          int idx_jru = j + r * J + u * J * S;
          for (int i = 0; i < J; i++) {
            int idx_ijru = i + j * J + r * J * J + u * S * J * J;
            multres[idx_ijru] = T_hat(i, r) *
              pow(c_hat(i, r) * tau_hat[idx_ijru], -epsilon[r]);
            multres_tmp[idx_jru] += pi[idx_ijru] * multres[idx_ijru];
          }
          if(multres_tmp[idx_jru] != 0) {
            P_hat_new[idx_jru] = pow(multres_tmp[idx_jru], - 1 / epsilon[r]);
          } else {
            P_hat_new[idx_jru] = 1;
          }
          P_hat_div[idx_jru] = P_hat_new[idx_jru] - P_hat[idx_jru];

          for (int i = 0; i < J; i++) {
            int idx_ijru = i + j * J + r * J * J + u * S * J * J;
            if(multres_tmp[idx_jru] != 0) {
              pi_hat[idx_ijru] =  multres[idx_ijru] / multres_tmp[idx_jru];
            } else {
              pi_hat[idx_ijru] = 1;
            }
          }
        }
      }
    }

    // update revenues
    R_hat_new.fill(0.0);
#pragma omp parallel for num_threads(nthreads)
    for (int i = 0; i < J; i++) {
      for (int r = 0; r < S; r++) {
        for (int j = 0; j < J; j++) {
          for (int u = 0; u < S; u++) {
            int idx_ijru = i + j * J + r * J * J + u * J * J * S;
            R_hat_new(i, r) += pi_hat[idx_ijru] * pi[idx_ijru]  *
              M_prime[j + r * J + u * J * S];
          }
          int idx_ijrS = i + j * J + r * J * J + S * J * J * S;
          R_hat_new(i, r) += pi_hat[idx_ijrS] * pi[idx_ijrS] * E_prime(j, r);
        }
        if (R(i, r) != 0) {
          R_hat_new(i, r) *= 1 / R(i, r);
        } else {
          R_hat_new(i, r) = 1;
        }
        // get changes for convergence test
        R_div(i, r) = (R_hat_new(i, r) - R_hat(i, r)) * R(i, r);
        R_hat_new(i, r) = R_hat(i, r) + zeta * (R_hat_new(i, r) - R_hat(i, r));
      }
    }

    Y_world_prime = 0;
    for (int i = 0; i < J; i++) {
      for (int r = 0; r < S; r++) {
        Y_world_prime += gamma_js(i, r) * R(i, r) * R_hat_new(i, r);
      }
    }

    // scale to world GDP
    norm_sf =  Y_world / Y_world_prime;
    for (int j = 0; j < J; j++) {
      for (int s = 0; s < S; s++) {
        R_prime(j, s) = R(j, s) * R_hat_new(j, s) * norm_sf;
        if(R(j, s) != 0) {
          R_hat_new(j, s) = R_prime(j, s) / R(j, s);
        } else {
          R_hat_new(j,s) = 1;
        }
      }
    }

    // find largest change in R_hat or P_hat in this iteration
    double divmax = std::max(Rcpp::max(R_div), Rcpp::max(P_hat_div));
    double divmin = std::min(Rcpp::min(R_div), Rcpp::min(P_hat_div));

    // absolute
    if (divmax < -divmin) {
      div_new = -divmin;
    } else {
      div_new = divmax;
    }

    if(div_new > div & iter > 10) {
      iter_converged = 0;
      if (zeta > 0.02) {
        zeta = round(zeta * 8000) / 10000;
        //  Rcpp::Rcout << std::endl << "This iteration did not converge. " <<
        //    "Decreasing dampening by 20pct to: " << zeta <<
        //      " and rerun previous iteration" << std::endl;
        std::swap(div, div_old);
        std::swap(P_hat, P_hat_old);
        std::swap(R_hat, R_hat_old);
      } else {
        Rcpp::Rcout << "This iteration did not converge." <<
          "Dampening already very low. Let's try our luck and continue." <<
            std::endl;
        std::swap(div, div_old);
        std::swap(div, div_new);
        std::swap(P_hat, P_hat_old);
        std::swap(P_hat, P_hat_new);
        std::swap(R_hat, R_hat_old);
        std::swap(R_hat, R_hat_new);
      }
    } else {
      std::swap(div, div_old);
      std::swap(div, div_new);
      std::swap(P_hat, P_hat_old);
      std::swap(P_hat, P_hat_new);
      std::swap(R_hat, R_hat_old);
      std::swap(R_hat, R_hat_new);
      iter_converged += 1;
      if (iter_converged == 10) {
        iter_converged = 0;
        if (zeta < 0.7) {
          zeta = round(zeta * 13000) / 10000;
          //     Rcpp:Rcout << std::endl << "Converged for 10 iterations." <<
          //    "Trying a dampening increase by 30% to: " << zeta << std::endl;
        }
      }
    }

    if(iter % 100 == 0) {
      Rcpp::checkUserInterrupt();
    }
    Rcpp::Rcout << "Iteration: " << iter <<"  The current div is:  " << div <<
      " and dampening: " << zeta << "\r" << std::flush;

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


  return List::create(Named("R_hat") = R_hat,
                      Named("P_hat") = P_hat,
                      Named("multres") = multres,
                      Named("multres_tmp") = multres_tmp,
                      Named("R_div") = R_div,
                      Named("P_hat_div") = P_hat_div,
                      Named("C_hat") = C_hat,
                      Named("Q_hat") = Q_hat,
                      Named("D_prime") = D_prime,
                      Named("Y_prime") = Y_prime,
                      Named("L_hat") = L_hat,
                      Named("L_tilde_hat") = L_tilde_hat,
                      Named("c_hat") = c_hat,
                      Named("pi_hat") = pi_hat);
}

