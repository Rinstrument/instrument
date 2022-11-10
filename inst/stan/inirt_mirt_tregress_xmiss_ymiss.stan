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
  matrix[N_long_obs, K] x;   // design matrix for predictors in latent regression model
  int<lower=1> D;        // number of first-order dimensions
  int<lower=1> nDelta;        // total number of delta parameters
  int<lower=1> L;        // number of non-zero loadings (alpha parameters)
  int<lower=1> Lbeta;    // number of regression parameters
  int<lower=1> beta_dstart[D]; // beta start index for each dimension
  int<lower=1> beta_dend[D];   // beta end index for each dimension
  int<lower=1> nobeta_dstart[D]; // beta start index for each dimension
  int<lower=1> nobeta_dend[D];   // beta end index for each dimension
  real weights[N_long_obs]; // weights for each observation
  // int<lower=1> Lxmiss;         // number of missing x values
  matrix[N_long_obs, K] x_miss;    // missing x index matrix (1 if missing, 0 else)
  int reg_miss[N, K];       // id value of missing x within a matrix, 0 else
  int<lower=0,upper=1> x_in_row_is_missing[N_long_obs]; // any missing x's in given row? for efficiency
}
parameters {
  matrix[N, D] theta;              // ability
  vector[nDelta] delta_l;          // difficulty
  vector<lower=0>[L] alpha_l;      // distrimination over multiple dimensions
  vector[Lbeta] beta_l;            // regression parameters for each dimension
}
transformed parameters {
  matrix[D, J] alpha; // connstrain the upper traingular elements to zero 
  matrix[Lbeta, D] beta; // organize regression parameters into a matrix
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
    int b_lower = 0;
    int b_upper = 0;
    for(d in 1:D) {
      b_lower = nobeta_dstart[d];
      b_upper = nobeta_dend[d];
      for(i in b_lower:b_upper) {
        beta[i, d] = 0;
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
  beta_l ~ normal(0, 5);
  {
    vector[N_long_obs] nu;
    for (i in 1:N_long_obs) {
      row_vector[D] xb;
      if(x_in_row_is_missing[i]) {
        for(k in 1:K) {
          for(d in 1:D) {
            if(x_miss[i, k]) {
              xb[d] += 0.0;
            } else {
              xb[d] += x[nn[i], k] * beta[k,d];
            }
          }
        }
      } else {
        xb = x[nn[i], ] * beta;
      }
      nu[i] = ((theta[nn[i], ] + xb)*col(alpha, jj[i]));
      target += ordered_logistic_lpmf(y[i] | nu[i], delta_trans[jj[i]]) * weights[i];
    }
  }
}
