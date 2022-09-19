#define RCPPDIST_DONT_USE_ARMA
#include <RcppDist.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(Rcpp, RcppDist)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include "abs2.h"

// [[Rcpp::export]]
double lgp2PR(IntegerMatrix & x,
              NumericMatrix & z,
              NumericVector lambda,
              int nB,
              int nT,
              int n,
              int J,
              NumericVector tJ,
              int nDmax,
              NumericVector lJ,
              NumericVector theta,
              NumericVector d,
              NumericVector a,
              NumericVector beta,
              double eps) {

  double llk = 0;
  d.attr("dim") = Dimension(nDmax, J);
  NumericMatrix dne = as<NumericMatrix>(d);

  theta.attr("dim") = Dimension(n, nT);
  NumericMatrix theta_mat = as<NumericMatrix>(theta);
  NumericMatrix d2(nDmax + 1, J);
  for(int i = 1; i <= nDmax; i++) {
    for(int j = 0; j < J; j++) {
      if(i <= lJ(j)) {
        d2(i, j) = dne(i - 1, j);
        double dt = d2(i, j);
        llk = llk + std::log(R::dnorm(dt, 0.0, 2.0, false) + eps);
      }
    }
  }
  
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < J; j++) {
      NumericVector r = cumsum((theta_mat(i, tJ(j)) * a(j)) - d2(_, j));
      NumericVector er = exp(r);
      NumericVector pr = er / Rcpp::sum(er);
      int xval = x(i, j);
      double lpr = std::log(pr(xval) + eps);
      llk += lpr;
    }
  }

  for(int i = 0; i < n; i++) {
    double thetagt = theta_mat(i, nT - 1);
    for(int tt = 0; tt < (nT - 1); tt++) {
      double thetait = theta_mat(i, tt);
      double lambdat = lambda(tt);
      llk = llk + std::log(R::dnorm(thetait - (thetagt * lambdat), 0.0, 1.0, false) + eps);
    }
  }

  for(int i = 0; i < n; i++) {
    double thetagt = theta_mat(i, nT - 1);
    NumericVector zi = z(i, _ );
    double zzsum = 0;
    for(int zz = 0; zz < nB; zz++) {
      zzsum += zi(zz) * beta(zz);
    }
    llk = llk + std::log(R::dnorm(thetagt - zzsum, 0.0, 1.0, false) + eps);
  }

  for(int l = 0; l < nT - 1; l++) {
    double lambdal = lambda(l);
    //llk = llk + std::log(R::dnorm(lambdal, 0.0, 1.0, false) + eps);
    llk = llk + std::log(d_truncnorm(lambdal, 0.0, 10.0, -10, 10) + eps);
  }

  for(int j = 0; j < J; j++) {
    double at = a(j);
    llk = llk + std::log(d_truncnorm(at, 0.0, 2.0, 0.0, 10.0) + eps);
  }
  
  for(int zz = 0; zz < nB; zz++) {
    double betaz = beta(zz);
    llk = llk + std::log(R::dnorm(betaz, 0.0, 10.0, false) + eps);
  }

  return llk;
}

// [[Rcpp::export]]
double lt2PR(IntegerMatrix & x,
             NumericMatrix & z,
             int iter,
             int burn,
             double delta,
             NumericMatrix & post,
             NumericVector & mean_theta,
             NumericVector & mean_theta_sq,
             NumericMatrix & draw,
             NumericMatrix & draw_theta,
             NumericVector ix,
             NumericVector ixe,
             int npar,
             int ntheta,
             int n,
             int nB,
             int J,
             int nDmax,
             NumericVector lJ,
             int nT,
             NumericVector tJ,
             NumericVector & accept,
             double eps,
             bool display_progress = true) {
  
  Progress p(iter, display_progress);
  
  NumericVector oldpars = draw(0, _ );
  NumericVector oldpars_theta = draw_theta(0, _ );

  for(int it = 1; it < iter; it++) {
    NumericVector prop = Rcpp::rnorm(npar, 0.0, delta);
    NumericVector newpars = oldpars + prop;
    
    NumericVector prop_theta = Rcpp::rnorm(ntheta, 0.0, delta);
    NumericVector newpars_theta = oldpars_theta + prop_theta;
    
    for(int q = ix(2) - 1; q < ixe(2); q ++){
      oldpars(q) = abs2(oldpars(q));
      newpars(q) = abs2(newpars(q));
    }
    
    double numer = lgp2PR(x,
                          z,
                          newpars[Range(ix(0) - 1, ixe(0) - 1)],
                          nB,
                          nT,
                          n,
                          J,
                          tJ,
                          nDmax,
                          lJ,
                          newpars_theta,
                          newpars[Range(ix(1) - 1, ixe(1) - 1)],
                          newpars[Range(ix(2) - 1, ixe(2) - 1)],
                          newpars[Range(ix(3) - 1, ixe(3) - 1)],
                          eps);

    double denom = lgp2PR(x,
                          z,
                          oldpars[Range(ix(0) - 1, ixe(0) - 1)],
                          nB,
                          nT,
                          n,
                          J,
                          tJ,
                          nDmax,
                          lJ,
                          oldpars_theta,
                          oldpars[Range(ix(1) - 1, ixe(1) - 1)],
                          oldpars[Range(ix(2) - 1, ixe(2) - 1)],
                          oldpars[Range(ix(3) - 1, ixe(3) - 1)],
                          eps);
    
    double acceptp = std::exp(numer - denom);
    double acceptit = (acceptp > R::runif(0.0, 1.0));
    
    if(acceptit == true) {
      oldpars = newpars;
      oldpars_theta = newpars_theta;
      if(it >= burn) {
        post(it - burn, _ ) = newpars;
        //llikelihood(it - burn) = get_llike(newpars, newpars_theta);
        mean_theta = mean_theta + newpars_theta;
        mean_theta_sq = mean_theta_sq + (newpars_theta * newpars_theta);
      }
      accept[it] = 1;
    } else {
      if(it >= burn) {
        post(it - burn, _ ) = oldpars;
        mean_theta = mean_theta + oldpars_theta;
        mean_theta_sq = mean_theta_sq + (oldpars_theta * oldpars_theta);
      }
      accept[it] = 0;
    }
    
    p.increment();
  }
  
  return 1.0;
}
