#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "densities.h"

// [[Rcpp::export]]
double lp_2pl_logit(arma::vec & x, arma::mat & data) {
  double eps = 1e-16;
  int n = data.n_rows;
  int j = data.n_cols;
  double theta_prob = 0.0;
  double lp = 0.0;
  for(int jj = 0; jj < j; jj++) {
    for(int i = 0; i < n; i++) {
      theta_prob = 1.0 / (1.0 + std::exp(-((x(jj) * x(i + (2*j))) - x(jj + j))));
      lp += std::log(((theta_prob * data(i, jj)) + ((1.0 - theta_prob) * (1.0 - data(i, jj)))) + eps);
      if(jj == 0) {
        lp += log_normd_dx(x(i + (2*j)), 0.0, 1.0);
      }
    }
    lp += logtruncnorm_dens_dx(x(jj), 0.0, 2.0, 0.0, 10.0);
    lp += log_normd_dx(x(jj + j), 0.0, 5.0);
  }
  return lp;
}

// [[Rcpp::export]]
double add_2(double x) {
  if(std::isnan(x)) {
    return 0.0;
  } else {
    return 10.1;
  }
}

// // [[Rcpp::export]]
// double lp_2pl_logit_mmis(arma::vec & x, arma::mat & data, arma::vec & x_missing, arma::mat & missing) {
//   double eps = 1e-16;
//   int n = data.n_rows;
//   int j = data.n_cols;
//   double theta_prob = 0.0;
//   double lp = 0.0;
//   for(int jj = 0; jj < j; jj++) {
//     for(int i = 0; i < n; i++) {
//       if(std::isnan(data(i, jj))) {
//         missing(i, jj) = 
//       }
//       theta_prob = 1.0 / (1.0 + std::exp(-((x(jj) * x(i + (2*j))) - x(jj + j))));
//       lp += std::log(((theta_prob * data(i, jj)) + ((1.0 - theta_prob) * (1.0 - data(i, jj)))) + eps);
//       if(jj == 0) {
//         lp += log_normd_dx(x(i + (2*j)), 0.0, 1.0);
//       }
//     }
//     lp += logtruncnorm_dens_dx(x(jj), 0.0, 2.0, 0.0, 10.0);
//     lp += log_normd_dx(x(jj + j), 0.0, 5.0);
//   }
//   return lp;
// }

// // [[Rcpp::export]]
// double lp_2pl_logit_reg(arma::vec & x, arma::mat & data, int p) {
//   int n = data.n_rows;
//   int dj = data.n_cols;
//   int j = dj - p;
//   double lp = 0.0;
//   for(int i = j; i < dj; i++) {
//     for(int nn = 0; nn < n; nn++) {
//       x(2*j + nn) += x(2*j + n + i - j) * data(nn,i);
//     }
//     lp += log_normd_dx(x(2*j + n + i - j), 0.0, 10.0);
//   }
//   arma::vec x_ = x.subvec(0, 2*j + n - 1);
//   arma::mat data_ = data.cols(0, j - 1);
//   lp += lp_2pl_logit(x_, data_);
//   return lp;
// }

// [[Rcpp::export]]
double lp_2pl_logit_reg(arma::vec & x, arma::mat & data, int p) {
  double eps = 1e-16;
  int n = data.n_rows;
  int j = data.n_cols;
  int jnp = j - p;
  double theta_prob = 0.0;
  double lp = 0.0;
  for(int jj = 0; jj < jnp; jj++) {
    for(int i = 0; i < n; i++) {
      theta_prob = 1.0 / (1.0 + std::exp(-((x(jj) * x(i + (2*jnp))) - x(jj + jnp))));
      lp += std::log(((theta_prob * data(i, jj)) + ((1.0 - theta_prob) * (1.0 - data(i, jj)))) + eps);
      if(jj == 0) {
        double sum_bx = 0.0;
        for(int pp = 0; pp < p; pp++) {
          sum_bx += x(2*jnp + n + pp) * data(i,jnp + pp);
        }
        lp += log_normd_dx(x(i + (2*jnp)) - sum_bx, 0.0, 1.0);
      }
    }
    lp += logtruncnorm_dens_dx(x(jj), 0.0, 2.0, 0.0, 10.0);
    lp += log_normd_dx(x(jj + jnp), 0.0, 5.0);
  }
  for(int i = jnp; i < j; i++) {
    // for(int nn = 0; nn < n; nn++) {
    //   x(2*jnp + nn) += x(2*jnp + n + i - jnp) * data(nn,i);
    // }
    lp += log_normd_dx(x(2*jnp + n + i - jnp), 0.0, 10.0);
  }
  return lp;
}

// [[Rcpp::export]]
double lp_lm(arma::vec & x, arma::mat & data, int n_mis) {
  int n = data.n_rows;
  int j = data.n_cols;
  //arma::vec ru = arma::randu(n_mis, arma::distr_param(0, 1));
  double lp = 0.0;
  double imp_data = 0.0;
  for(int nn = 0; nn < n; nn++) {
    int miss_index = 0;
    if(std::isnan(data(nn,0))) {
      //imp_data = ru(miss_index) > (1.0 - x(2 + miss_index));
      imp_data = x(2 + miss_index);
      lp += log_normd_dx(imp_data, x(0), x(1));
      miss_index += 1;
    } else {
      lp += log_normd_dx(data(nn,0), x(0), x(1));
    }
  }
  for(int jj = 0; jj < j; jj++) {
    lp += log_normd_dx(x(jj), 0.0, 5.0);
  }
  lp += log_unifd_dx(x(1), 0.0, 10.0);
  for(int i = 0; i < n_mis; i++) {
    lp += log_normd_dx(x(2 + i), 0.0, 5.0);
  }
  return lp;
}

// [[Rcpp::export]]
double lp_2pl_ho2l(arma::vec & x, arma::mat & data, int p) {
  double dim_1order = data(0,0);
  arma::vec dim_lengths = arma::vec(dim_1order);
  for(int i = 0; i < dim_1order; i++) {
    dim_lengths(i) = data(0, 1 + i);
  } // dont actually need this
  double j = data.n_cols - dim_1order - 1;
  double eps = 1e-16;
  int n = data.n_rows;
  double theta_prob = 0.0;
  double lp = 0.0;
  int d_index = 0;
  for(int d = 0; d < dim_1order; d++) {
    for(int jj = d_index; jj < (dim_lengths(d) + d_index); jj++) {
      for(int i = 0; i < n; i++) {
        theta_prob = 1.0 / (1.0 + std::exp(-((x(jj) * x(i + 2*j + d*n)) - x(jj + j))));
        lp += std::log(((theta_prob * data(i, jj + dim_1order + 1)) + ((1.0 - theta_prob) * (1.0 - data(i, jj + dim_1order + 1)))) + eps);
        if(jj == d_index) {
          lp += log_normd_dx(x(i + 2*j + d*n) - x(2*j + (dim_1order + 1)*n + d) * x(i + dim_1order*n + 2*j), 0.0, 1.0);
          if(jj == 0) {
            lp += log_normd_dx(x(i + dim_1order*n + 2*j), 0.0, 1.0);
          }
        }
      }
      lp += logtruncnorm_dens_dx(x(jj), 0.0, 2.0, 0.0, 10.0);
      lp += log_normd_dx(x(jj + j), 0.0, 5.0);
    }
    d_index += dim_lengths(d);
    lp += log_unifd_dx(x((2*j) + (dim_1order + 1)*n + d), -1.0, 1.0);
  }
  return lp;
}
