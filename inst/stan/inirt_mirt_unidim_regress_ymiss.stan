data {
  int<lower=1> N;              // number of participants
  int<lower=1> J;              // number of questions
  int<lower=1> K;              // number of regression parameters
  int<lower=2> Ncateg;         // Number of categories in ordered logistic responses
  int<lower=1> N_long_obs;              // number of long observations
  int<lower=1,upper=N> nn[N_long_obs];  // participant for observation n
  int<lower=1,upper=J> jj[N_long_obs];  // question for observation n
  int<lower=0,upper=Ncateg> y[N_long_obs];   // correctness for observation n
  matrix[N_long_obs, K] x;   // correctness for observation n
  int<lower=1> D;        // number of first-order dimensions
  int<lower=1> L;        // number of non-zero loadings
  int<lower=1> Lbeta;    // number of regression parameters
  int<lower=1> beta_dstart[D]; // beta start index for each dimension
  int<lower=1> beta_dend[D];   // beta end index for each dimension
  real weights[N_long_obs]; // weights for each observation
}
parameters {
  matrix[N, D] theta;              // ability
  ordered[Ncateg-1] delta[J];      // difficulty for k old: vector[J] delta;  
  // vector<lower=0>[J] alpha;     // discrimination of k
  vector<lower=0>[L] alpha_l;      // distrimination over multiple dimensions
  vector[Lbeta] beta_l;            // regression parameters for each dimension
}
transformed parameters {
  matrix[D, J] alpha; // connstrain the upper traingular elements to zero 
  matrix[Lbeta, D] beta; // organize regression parameters into a matrix

    for(j in 1:J) {
      for(d in (j+1):D) {
        alpha[d, j] = 0;
      }
    } 
    {
      int index = 0;
      for(d in 1:D) {
        for(j in d:J) {
          index = index + 1;
          alpha[d, j] = alpha_l[index];
        }
      }
    }
    
    {
      int bindex = 0;
      int b_lower = 0;
      int b_upper = 0;
      for(d in 1:D) {
        b_lower = beta_dstart[d];
        b_upper = beta_dend[d];
        for(i in b_lower:b_upper) {
          bindex = bindex + 1;
          beta[i, d] = beta_l[bindex];
        }
      }
    }
}
model {
  to_vector(theta) ~ normal(0, 1);
  for(j in 1:J) {
    for(i in 1:(Ncateg-1)) {
      delta[j][i] ~ normal(0, 1);
    }
  }
  alpha_l ~ lognormal(0, 0.3);
  beta_l ~ normal(0, 5);
  {
    vector[N_long_obs] nu;
    for (i in 1:N_long_obs) {
      nu[i] = ((theta[nn[i], ] + (x[nn[i], ] * beta))*col(alpha, jj[i]));
      target += ordered_logistic_lpmf(y[i] | nu[i], delta[jj[i]]) * weights[i];
    }
  }
}
