data {
  int<lower=1> N;              // number of participants
  int<lower=1> J;              // number of questions
  int<lower=1> K;              // number of regression parameters
  int<lower=2> Ncateg_max;         // Max Number of categories in ordered logistic responses
  int<lower=2,upper=Ncateg_max> Ncategi[J]; // Number of categories for each item
  int<lower=1> N_long_obs;              // number of long observations
  int<lower=1,upper=N> nn[N_long_obs];  // participant for observation n
  int<lower=1,upper=J> jj[N_long_obs];  // question for observation n
  int<lower=0,upper=Ncateg_max> y[N_long_obs];   // correctness for observation n
  int<lower=1> D;        // number of first-order dimensions
  int<lower=1> nDelta;        // total number of delta parameters
  int<lower=1> L;        // number of non-zero loadings
  real weights[N_long_obs]; // weights for each observation
}
parameters {
  matrix[N, D] theta;              // ability
  vector[nDelta] delta_l;          // difficulty
  vector<lower=0>[L] alpha_l;      // distrimination over multiple dimensions
}
transformed parameters {
  matrix[D, J] alpha; // connstrain the upper traingular elements to zero 
  vector[Ncateg_max-1] delta_trans[J]; // Make excess categories infinite

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
    int idx = 0;
    int d_index = 0;
    for(j in 1:J) {
      vector[Ncategi[j]-1] ds_ind = sort_asc(delta_l[(d_index+1):(d_index+Ncategi[j]-1)]);
      for(i in 1:(Ncategi[j]-1)) {
        d_index = d_index + 1;
        delta_trans[j][i] = ds_ind[i];
      }
      for(i in (Ncategi[j]):(Ncateg_max-1)) {
        delta_trans[j][i] = 1e7 + idx;
        idx = idx + 1;
      }
    }
  }
}
model {
  to_vector(theta) ~ normal(0, 1);
  alpha_l ~ lognormal(0, 0.3);
  delta_l ~ normal(0, 1);
  {
    vector[N_long_obs] nu;
    for (i in 1:N_long_obs) {
      nu[i] = ((theta[nn[i], ])*col(alpha, jj[i]));
      target += ordered_logistic_lpmf(y[i] | nu[i], delta_trans[jj[i]]) * weights[i];
    }
  }
}
