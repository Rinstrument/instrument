// #include <RcppArmadillo.h>
// #include "RcppArmadilloExtensions/sample.h"
// # include <RcppArmadilloExtensions/fixprob.h>
// // [[Rcpp::depends(RcppArmadillo)]]
// #include "densities.h"
// #include "random_distributions.h"
// #include "logadd.h"

// // [[Rcpp::export]]
// double lp_2pl_logit(arma::vec & x, arma::mat & data) {
//   double eps = 1e-16;
//   int n = data.n_rows;
//   int j = data.n_cols;
//   double theta_prob = 0.0;
//   double lp = 0.0;
//   for(int jj = 0; jj < j; jj++) {
//     for(int i = 0; i < n; i++) {
//       theta_prob = 1.0 / (1.0 + std::exp(-((x(jj) * x(i + (2*j))) - x(jj + j))));
//       // lp += std::log(((theta_prob * data(i, jj)) + ((1.0 - theta_prob) * (1.0 - data(i, jj)))) + eps);
//       lp += data(i, jj)*std::log(theta_prob + eps) + (1.0 - data(i, jj))*std::log(1.0 - theta_prob + eps);
//       if(jj == 0) {
//         lp += log_normd_dx(x(i + (2*j)), 0.0, 2.0);
//       }
//     }
//     lp += logtruncnorm_dens_dx(x(jj), 0.0, 1.0, 0.0, 5.0);
//     lp += log_normd_dx(x(jj + j), 0.0, 2.0);
//   }
//   return lp;
// }

// // [[Rcpp::export]]
// double add_2(double x) {
//   if(std::isnan(x)) {
//     return 0.0;
//   } else {
//     return 10.1;
//   }
// }

// // // [[Rcpp::export]]
// // double lp_2pl_logit_mmis(arma::vec & x, arma::mat & data, arma::vec & x_missing, arma::mat & missing) {
// //   double eps = 1e-16;
// //   int n = data.n_rows;
// //   int j = data.n_cols;
// //   double theta_prob = 0.0;
// //   double lp = 0.0;
// //   for(int jj = 0; jj < j; jj++) {
// //     for(int i = 0; i < n; i++) {
// //       if(std::isnan(data(i, jj))) {
// //         missing(i, jj) = 
// //       }
// //       theta_prob = 1.0 / (1.0 + std::exp(-((x(jj) * x(i + (2*j))) - x(jj + j))));
// //       lp += std::log(((theta_prob * data(i, jj)) + ((1.0 - theta_prob) * (1.0 - data(i, jj)))) + eps);
// //       if(jj == 0) {
// //         lp += log_normd_dx(x(i + (2*j)), 0.0, 1.0);
// //       }
// //     }
// //     lp += logtruncnorm_dens_dx(x(jj), 0.0, 2.0, 0.0, 10.0);
// //     lp += log_normd_dx(x(jj + j), 0.0, 5.0);
// //   }
// //   return lp;
// // }

// // // [[Rcpp::export]]
// // double lp_2pl_logit_reg(arma::vec & x, arma::mat & data, int p) {
// //   int n = data.n_rows;
// //   int dj = data.n_cols;
// //   int j = dj - p;
// //   double lp = 0.0;
// //   for(int i = j; i < dj; i++) {
// //     for(int nn = 0; nn < n; nn++) {
// //       x(2*j + nn) += x(2*j + n + i - j) * data(nn,i);
// //     }
// //     lp += log_normd_dx(x(2*j + n + i - j), 0.0, 10.0);
// //   }
// //   arma::vec x_ = x.subvec(0, 2*j + n - 1);
// //   arma::mat data_ = data.cols(0, j - 1);
// //   lp += lp_2pl_logit(x_, data_);
// //   return lp;
// // }

// // [[Rcpp::export]]
// double lp_2pl_logit_reg(arma::vec & x, arma::mat & data, int p) {
//   double eps = 1e-16;
//   int n = data.n_rows;
//   int j = data.n_cols;
//   int jnp = j - p;
//   double theta_prob = 0.0;
//   double lp = 0.0;
//   for(int jj = 0; jj < jnp; jj++) {
//     for(int i = 0; i < n; i++) {
//       theta_prob = 1.0 / (1.0 + std::exp(-((x(jj) * x(i + (2*jnp))) - x(jj + jnp))));
//       lp += std::log(((theta_prob * data(i, jj)) + ((1.0 - theta_prob) * (1.0 - data(i, jj)))) + eps);
//       if(jj == 0) {
//         double sum_bx = 0.0;
//         for(int pp = 0; pp < p; pp++) {
//           sum_bx += x(2*jnp + n + pp) * data(i,jnp + pp);
//         }
//         lp += log_normd_dx(x(i + (2*jnp)) - sum_bx, 0.0, 1.0);
//       }
//     }
//     lp += logtruncnorm_dens_dx(x(jj), 0.0, 2.0, 0.0, 10.0);
//     lp += log_normd_dx(x(jj + jnp), 0.0, 5.0);
//   }
//   for(int i = jnp; i < j; i++) {
//     // for(int nn = 0; nn < n; nn++) {
//     //   x(2*jnp + nn) += x(2*jnp + n + i - jnp) * data(nn,i);
//     // }
//     lp += log_normd_dx(x(2*jnp + n + i - jnp), 0.0, 10.0);
//   }
//   return lp;
// }

// // [[Rcpp::export]]
// double lp_lm(arma::vec & x, arma::mat & data, int n_mis) {
//   int n = data.n_rows;
//   int j = data.n_cols;
//   //arma::vec ru = arma::randu(n_mis, arma::distr_param(0, 1));
//   double lp = 0.0;
//   double imp_data = 0.0;
//   for(int nn = 0; nn < n; nn++) {
//     int miss_index = 0;
//     if(std::isnan(data(nn,0))) {
//       //imp_data = ru(miss_index) > (1.0 - x(2 + miss_index));
//       imp_data = x(2 + miss_index);
//       lp += log_normd_dx(imp_data, x(0), x(1));
//       miss_index += 1;
//     } else {
//       lp += log_normd_dx(data(nn,0), x(0), x(1));
//     }
//   }
//   for(int jj = 0; jj < j; jj++) {
//     lp += log_normd_dx(x(jj), 0.0, 5.0);
//   }
//   lp += log_unifd_dx(x(1), 0.0, 10.0);
//   for(int i = 0; i < n_mis; i++) {
//     lp += log_normd_dx(x(2 + i), 0.0, 5.0);
//   }
//   return lp;
// }

// // [[Rcpp::export]]
// double lp_2pl_ho2l(arma::vec & x, arma::mat & data, int p) {
//   double dim_1order = data(0,0);
//   arma::vec dim_lengths = arma::vec(dim_1order);
//   for(int i = 0; i < dim_1order; i++) {
//     dim_lengths(i) = data(0, 1 + i);
//   } // dont actually need this
//   double j = data.n_cols - dim_1order - 1;
//   //Rcpp::Rcout << "j: " << j << std::endl;
//   double eps = 1e-16;
//   int n = data.n_rows;
//   double theta_prob = 0.0;
//   double lp = 0.0;
//   int d_index = 0;
//   for(int d = 0; d < dim_1order; d++) {
//     //Rcpp::Rcout << "d_index: " << d_index << std::endl;
//     //Rcpp::Rcout << "(dim_lengths(d) + d_index): " << (dim_lengths(d) + d_index) << std::endl;
//     for(int jj = d_index; jj < (dim_lengths(d) + d_index); jj++) {
//       for(int i = 0; i < n; i++) {
//         theta_prob = 1.0 / (1.0 + std::exp(-((x(jj) * x(i + 2*j + d*n)) - x(jj + j))));
//         lp += std::log(((theta_prob * data(i, jj + dim_1order + 1)) + ((1.0 - theta_prob) * (1.0 - data(i, jj + dim_1order + 1)))) + eps);
//         if(jj == d_index) {
//           lp += log_normd_dx(x(i + 2*j + d*n), x(2*j + (dim_1order + 1)*n + d) * x(i + dim_1order*n + 2*j), 1.0);
//           if(jj == 0) {
//             lp += log_normd_dx(x(i + dim_1order*n + 2*j), 0.0, 1.0);
//           }
//         }
//       }
//       lp += logtruncnorm_dens_dx(x(jj), 1.0, 2.0, 0.0, 8.0);
//       lp += log_normd_dx(x(jj + j), 0.0, 5.0);
//     }
//     d_index += dim_lengths(d);
//     lp += log_unifd_dx(x(2*j + (dim_1order + 1)*n + d), -1.0, 1.0);
//     // lp += logtruncnorm_dens_dx(x(2*j + (dim_1order + 1)*n + d), 0.0, 5.0, -10.0, 10.0);
//        //     Rcpp::Rcout << "2*j + (dim_1order + 1)*n + d: " << 2*j + (dim_1order + 1)*n + d << std::endl;

//       //  Rcpp::Rcout << "x(2*j + (dim_1order + 1)*n + d): " << x(2*j + (dim_1order + 1)*n + d) << std::endl;

//   }
//   //Rcpp::Rcout << "lp: " << lp << std::endl;

//   return lp;
// }



// // [[Rcpp::export]]
// arma::vec sir_irt(arma::mat & data, int m_iterations, int n_iterations) {
//   //int iterations = 1e6;
//   double eps = 1e-16;
//   int n = data.n_rows;
//   int j = data.n_cols;
//   int p = 2*j + n;
//   arma::ivec samp_indices = arma::ivec(p);
//   //arma::ivec sindex = arma::randi(m_iterations, arma::distr_param(+(0), +(p-1)));
//   arma::ivec sindex = arma::regspace<arma::ivec>(0,  m_iterations - 1);
//   //Rcpp::Rcout << sindex << std::endl;
//   arma::mat x = arma::mat(p, m_iterations, arma::fill::zeros);
//   arma::ivec seeds = arma::randi<arma::ivec>(p*m_iterations, arma::distr_param(+(1), +(123456789)));
//   // compute standardized importance weights
//   arma::vec w = arma::vec(m_iterations, arma::fill::zeros);

  
//   for(int iter = 0; iter < m_iterations; iter++) {
//     double lp = 0.0;
//     for(int jj = 0; jj < j; jj++) {
//       for(int i = 0; i < n; i++) {
//         if(jj == 0) {
//           x(i + (2*j), iter) = normal_01_sample(seeds(p*iter + i + (2*j)))*3.0;
//           lp += log_normd_dx(x(i + (2*j), iter), 0.0, 1.0);
//         }
//       }
//       x(jj, iter) = truncated_normal_ab_sample(0.5, 5.0, 0.0, 20.0, seeds(p*iter + jj));
//       x(jj + j, iter) = normal_01_sample(seeds(p*iter + j + jj))*10.0;
//       lp += logtruncnorm_dens_dx(x(jj, iter), 0.5, 5.0, 0.0, 10.0);
//       lp += log_normd_dx(x(jj + j, iter), 0.0, 5.0);
//     }

//     double theta_prob = 0.0;
//     double lk = 0.0;
//     for(int jj = 0; jj < j; jj++) {
//       for(int i = 0; i < n; i++) {
//         // theta_prob = 1.0 / (1.0 + std::exp(-((x(jj, iter) * x(i + (2*j), iter)) - x(jj + j, iter))));
//         // theta_prob = 1.0 / (1.0 + std::exp(-((x(jj, iter)*x(i + (2*j), iter) - x(jj + j, iter)))));
//         theta_prob = 1.0 / (1.0 + std::exp(-(x(jj, iter)*(x(i + (2*j), iter) - x(jj + j, iter)))));

//         // lp += std::log(((theta_prob * data(i, jj)) + ((1.0 - theta_prob) * (1.0 - data(i, jj)))) + eps);
//         lk += data(i, jj)*std::log(theta_prob + eps) + (1.0 - data(i, jj))*std::log(1.0 - theta_prob + eps);
//       }
//     }

//     // w(iter) = lk - lp;
//     w(iter) = lk;

//     // Rcpp::Rcout << "lk: " <<  lk << std::endl;
//     // Rcpp::Rcout << "lp: " << lp << std::endl;
//     // Rcpp::Rcout << "w: " << w(iter) << std::endl;
    
//   }

//   // standardize importance weights

//   // resample with replacement
  
//   // w = w + std::abs(arma::min(w)) + 1.0;
//   w = arma::exp(w / arma::min(w));
//   w = w / arma::accu(w);
//   Rcpp::Rcout << "weights: " << w << std::endl;

//   arma::vec p_means = arma::vec(p, arma::fill::zeros);
//   for(int iter = 0; iter < n_iterations; iter++) {
//     samp_indices = Rcpp::RcppArmadillo::sample(sindex, m_iterations, true, w);
//     for(int col_index = 0; col_index < m_iterations; col_index++) {
//       p_means += x.col(samp_indices(col_index));
//     }
//   }
//   p_means = p_means / (n_iterations * m_iterations);
//   //Rcpp::Rcout << samp_indices << std::endl;

//   return p_means;
// }



// // [[Rcpp::export]]
// arma::mat sir_logit(arma::mat & data, int m_iterations, int n_iterations) {
//   //int iterations = 1e6;
//   double eps = 1e-16;
//   int n = data.n_rows;
//   int j = data.n_cols;
//   int p = 1;
//   arma::ivec samp_indices = arma::ivec(p);
//   //arma::ivec sindex = arma::randi(m_iterations, arma::distr_param(+(0), +(p-1)));
//   arma::ivec sindex = arma::regspace<arma::ivec>(0,  m_iterations - 1);
//   //Rcpp::Rcout << sindex << std::endl;
//   arma::mat x = arma::mat(p, m_iterations, arma::fill::zeros);
//   arma::ivec seeds = arma::randi<arma::ivec>(p*m_iterations, arma::distr_param(+(1), +(999999999)));
//   // compute standardized importance weights
//   arma::vec w = arma::vec(m_iterations, arma::fill::zeros);

  
//   for(int iter = 0; iter < m_iterations; iter++) {
//     double lp = 0.0;
//     // // for(int jj = 0; jj < j; jj++) {
//     //   for(int i = 0; i < n; i++) {
//     //     // if(jj == 0) {
//     //       x(i + (2*j), iter) = normal_01_sample(seeds(p*iter + i + (2*j)))*3.0;
//     //       lp += log_normd_dx(x(i + (2*j), iter), 0.0, 1.0);
//     //     // }
//     //   // }
//     //   x(jj, iter) = truncated_normal_ab_sample(0.5, 5.0, 0.0, 20.0, seeds(p*iter + jj));
//     //   x(jj + j, iter) = normal_01_sample(seeds(p*iter + j + jj))*10.0;
//     //   lp += logtruncnorm_dens_dx(x(jj, iter), 0.5, 5.0, 0.0, 10.0);
//     //   lp += log_normd_dx(x(jj + j, iter), 0.0, 5.0);
//     // }

//     x(0, iter) = normal_01_sample(seeds(p*iter))*10.0;
//     lp += log_normd_dx(x(0, iter), 0.0, 10.0);

//     double theta_prob = 0.0;
//     double lk = 0.0;
//     // for(int jj = 0; jj < j; jj++) {
//       for(int i = 0; i < n; i++) {
//         theta_prob = 1.0 / (1.0 + std::exp(-(x(0, iter))));
//         lk += data(i, 0)*std::log(theta_prob + eps) + (1.0 - data(i, 0))*std::log(1.0 - theta_prob + eps);
//       }
    
//     w(iter) = lk; //- lp;
//   }

//   double md = arma::max(w);
//   arma::vec lw = w - md;
//   for(int l_i = 0; l_i < m_iterations; l_i++) {
//     if(std::isnan(lw(l_i))) {
//       lw(l_i) = arma::min(lw);
//     }
//   }
//   arma::vec probs = arma::exp(lw - logadd(lw));

//   arma::mat p_means = arma::mat(p, n_iterations, arma::fill::zeros);
//   for(int iter = 0; iter < n_iterations; iter++) {
//     // Rcpp::RcppArmadillo::FixProb(probs, m_iterations, true);
//     samp_indices = Rcpp::RcppArmadillo::sample(sindex, m_iterations, true, probs);
//     for(int col_index = 0; col_index < m_iterations; col_index++) {
//       p_means.col(iter) += x.col(samp_indices(col_index));
//     }
//     p_means.col(iter) = p_means.col(iter) / m_iterations;
//   }
//   return probs;
// }



// // [[Rcpp::export]]
// arma::vec sir_irt2(arma::mat & data, int m_iterations, int n_iterations) {
//   //int iterations = 1e6;
//   double eps = 1e-16;
//   int n = data.n_rows;
//   int j = data.n_cols;
//   int p = 2*j + n;
//   arma::ivec samp_indices = arma::ivec(p);
//   //arma::ivec sindex = arma::randi(m_iterations, arma::distr_param(+(0), +(p-1)));
//   arma::ivec sindex = arma::regspace<arma::ivec>(0,  m_iterations - 1);
//   //Rcpp::Rcout << sindex << std::endl;
//   arma::mat x = arma::mat(p, m_iterations, arma::fill::zeros);
//   arma::ivec seeds = arma::randi<arma::ivec>(p*m_iterations, arma::distr_param(+(1), +(999999999)));
//   // compute standardized importance weights
//   arma::vec w = arma::vec(m_iterations, arma::fill::zeros);

//   for(int iter = 0; iter < m_iterations; iter++) {
//     double lp = 0.0;
//     for(int jj = 0; jj < j; jj++) {
//       for(int i = 0; i < n; i++) {
//         if(jj == 0) {
//           x(i + (2*j), iter) = normal_01_sample(seeds(p*iter + i + (2*j)))*1.0;
//         }
//       }
//       x(jj, iter) = truncated_normal_ab_sample(0.5, 2.0, 0.0, 10.0, seeds(p*iter + jj));
//       x(jj + j, iter) = normal_01_sample(seeds(p*iter + j + jj))*5.0;
//     }

//     double theta_prob = 0.0;
//     double lk = 0.0;
//     for(int jj = 0; jj < j; jj++) {
//       for(int i = 0; i < n; i++) {
//         theta_prob = 1.0 / (1.0 + std::exp(-((x(jj, iter)*x(i + (2*j), iter)) - x(jj + j, iter))));
//         // theta_prob = 1.0 / (1.0 + std::exp(-(x(jj, iter)*(x(i + (2*j), iter) - x(jj + j, iter)))));
//         lk += data(i, jj)*std::log(theta_prob + eps) + (1.0 - data(i, jj))*std::log(1.0 - theta_prob + eps);
//       }
//     }
    
//     w(iter) = lk; //- lp;
//   }

//   double md = arma::max(w);
//   Rcpp::Rcout << "md: " << md << std::endl;
//   arma::vec lw = w - md;
//   for(int l_i = 0; l_i < m_iterations; l_i++) {
//     if(std::isnan(lw(l_i))) {
//       lw(l_i) = arma::min(lw);
//     }
//   }
//   arma::vec probs = arma::exp(lw - logadd(lw));

//   arma::mat p_means = arma::mat(p, n_iterations, arma::fill::zeros);
//   for(int iter = 0; iter < n_iterations; iter++) {
//     Rcpp::RcppArmadillo::FixProb(probs, m_iterations, true);
//     samp_indices = Rcpp::RcppArmadillo::sample(sindex, m_iterations, true, probs);
//     for(int col_index = 0; col_index < m_iterations; col_index++) {
//       p_means.col(iter) += x.col(samp_indices(col_index));
//     }
//     p_means.col(iter) = p_means.col(iter) / m_iterations;
//   }
//   return lw;
// }



// // [[Rcpp::export]]
// arma::vec sir_irt3(arma::mat & data, int m_iterations, int n_iterations) {
//   //int iterations = 1e6;
//   double eps = 1e-16;
//   int n = data.n_rows;
//   int j = data.n_cols;
//   int p = 2*j + n;
//   arma::ivec samp_indices = arma::ivec(p);
//   //arma::ivec sindex = arma::randi(m_iterations, arma::distr_param(+(0), +(p-1)));
//   arma::ivec sindex = arma::regspace<arma::ivec>(0,  m_iterations - 1);
//   //Rcpp::Rcout << sindex << std::endl;
//   arma::mat x = arma::mat(p, m_iterations, arma::fill::zeros);
//   arma::ivec seeds = arma::randi<arma::ivec>(p*m_iterations, arma::distr_param(+(1), +(999999999)));
//   // compute standardized importance weights
//   arma::vec w = arma::vec(m_iterations, arma::fill::zeros);

//   for(int iter = 0; iter < m_iterations; iter++) {
//     double lp = 0.0;
//     for(int jj = 0; jj < j; jj++) {
//       for(int i = 0; i < n; i++) {
//         if(jj == 0) {
//           x(i + (2*j), iter) = normal_01_sample(seeds(p*iter + i + (2*j)))*10.0;
//         }
//       }
//       x(jj, iter) = truncated_normal_ab_sample(0.5, 10.0, 0.0, 30.0, seeds(p*iter + jj));
//       x(jj + j, iter) = normal_01_sample(seeds(p*iter + j + jj))*10.0;
//     }

//     double theta_prob = 0.0;
//     double lk = 0.0;
//     for(int jj = 0; jj < j; jj++) {
//       for(int i = 0; i < n; i++) {
//         theta_prob = 1.0 / (1.0 + std::exp(-((x(jj, iter)*x(i + (2*j), iter)) - x(jj + j, iter))));
//         // theta_prob = 1.0 / (1.0 + std::exp(-(x(jj, iter)*(x(i + (2*j), iter) - x(jj + j, iter)))));
//         lk += data(i, jj)*std::log(theta_prob + eps) + (1.0 - data(i, jj))*std::log(1.0 - theta_prob + eps);
//       }
//     }
    
//     w(iter) = lk; //- lp;
//   }

//   double md = arma::max(w);
//   Rcpp::Rcout << "md: " << md << std::endl;
//   arma::vec lw = w - md;
//   for(int l_i = 0; l_i < m_iterations; l_i++) {
//     if(std::isnan(lw(l_i))) {
//       lw(l_i) = arma::min(lw);
//     }
//   }
//   arma::vec probs = arma::exp(lw - logadd(lw));
//   // for(int l_i = 0; l_i < m_iterations; l_i++) {
//   //   if(probs(l_i) == 1.0) {
//   //     probs(l_i) = 2e-06;
//   //   }
//   // }

//   arma::mat p_means = arma::mat(p, n_iterations, arma::fill::zeros);
//   for(int iter = 0; iter < n_iterations; iter++) {
//     Rcpp::RcppArmadillo::FixProb(probs, m_iterations, true);
//     samp_indices = Rcpp::RcppArmadillo::sample(sindex, m_iterations, true, probs);
//     for(int col_index = 0; col_index < m_iterations; col_index++) {
//       p_means.col(iter) += x.col(samp_indices(col_index));
//     }
//     p_means.col(iter) = p_means.col(iter) / m_iterations;
//   }
//   return probs;
// }