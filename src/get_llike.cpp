// #define RCPPDIST_DONT_USE_ARMA
// #include <RcppDist.h>
// using namespace Rcpp;
// // [[Rcpp::plugins(cpp11)]]
// // [[Rcpp::depends(Rcpp, RcppDist)]]
// 
// // [[Rcpp::export]]
// double get_llike(Rcpp::NumericVector pars, Rcpp::NumericVector pars_theta,
//                  Rcpp::NumericVector pars_index, Rcpp::NumericVector pars_theta_index, 
//                  int model) {
//   
//   
//   // model = 0: GPCM (2p) regression
//   
//   if(model == 0) {
//     
//     double llk = 0;
//     d.attr("dim") = Dimension(nDmax, J);
//     NumericMatrix dne = as<NumericMatrix>(d);
//     
//     theta.attr("dim") = Dimension(n, nT);
//     NumericMatrix theta_mat = as<NumericMatrix>(theta);
//     NumericMatrix d2(nDmax + 1, J);
//     for(int i = 1; i <= nDmax; i++) {
//       for(int j = 0; j < J; j++) {
//         if(i <= lJ(j)) {
//           d2(i, j) = dne(i - 1, j);
//           // double dt = d2(i, j);
//           // llk = llk + std::log(R::dnorm(dt, 0.0, 2.0, false) + eps);
//         }
//       }
//     }
//     
//     for(int i = 0; i < n; i++) {
//       for(int j = 0; j < J; j++) {
//         NumericVector r = cumsum((theta_mat(i, tJ(j)) * a(j)) - d2(_, j));
//         NumericVector er = exp(r);
//         NumericVector pr = er / Rcpp::sum(er);
//         int xval = x(i, j);
//         double lpr = std::log(pr(xval) + eps);
//         llk += lpr;
//       }
//     }
//     
//     
//   }
//   
//   return llk;
// }