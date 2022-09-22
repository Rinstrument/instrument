// [[Rcpp::export]]
double llk_1logit(arma::mat & x, arma::mat & data) {
  
  int n = data.n_rows;
  int j = data.n_cols;
  
  arma::vec p_vec = arma::vec(1);
  p_vec = 1 / (1 + arma::exp(-x));
  double llk = arma::accu(data * arma::log(p_vec + 0.0000000000001) + 
                            (1.0 - data) * arma::log(1 - p_vec + 0.0000000000001)) + 
    arma::accu(arma::log(arma::normpdf(x, 0.0, 10.0) + 0.0000000000001));
  
  return llk;
}