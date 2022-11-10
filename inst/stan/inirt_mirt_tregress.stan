data {
  int<lower=1> N;              // number of participants
  int<lower=1> J;              // number of questions
  int<lower=1> K;              // number of regression parameters
  int<lower=2> Ncateg_max;         // Max Number of categories in ordered logistic responses
  int<lower=2,upper=Ncateg_max> Ncategi[J]; // Number of categories for each item
  int<lower=1> N_long;              // number of long observations
  int<lower=1,upper=N> nn[N_long];  // participant for observation n
  int<lower=1,upper=J> jj[N_long];  // question for observation n
  int<lower=0,upper=Ncateg_max> y[N_long];   // correctness for observation n
  matrix[N_long, K] x;   // correctness for observation n
  int<lower=1> D;        // number of first-order dimensions
  int<lower=1> nDelta;        // total number of delta parameters
  int<lower=1> L;        // number of non-zero loadings
  int<lower=1> Lbeta;    // number of regression parameters
  int<lower=1> beta_dstart[D]; // beta start index for each dimension
  int<lower=1> beta_dend[D];   // beta end index for each dimension
  int<lower=1> nobeta_dstart[D]; // beta start index for each dimension
  int<lower=1> nobeta_dend[D];   // beta end index for each dimension
  real weights[N_long]; // weights for each observation
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
    vector[N_long] nu;
    for (i in 1:N_long) {
      nu[i] = ((theta[nn[i], ] + (x[nn[i], ] * beta))*col(alpha, jj[i]));
      target += ordered_logistic_lpmf(y[i] | nu[i], delta_trans[jj[i]]) * weights[i];
    }
  }
}
