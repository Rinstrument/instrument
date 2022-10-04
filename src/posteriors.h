#ifndef POSTERIOR
#define POSTERIOR

double lp_2pl_logit(arma::vec & x, arma::mat & data);
double lp_2pl_logit_reg(arma::vec & x, arma::mat & data, int p);

#endif