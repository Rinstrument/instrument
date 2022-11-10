data {
  int<lower=1> N;              // number of participants
  int<lower=1> J;              // number of questions
  int<lower=2> Ncateg_max;         // Max Number of categories in ordered logistic responses
  int<lower=2,upper=Ncateg_max> Ncategi[J]; // Number of categories for each item
  int<lower=1> N_long;              // number of long observations
  int<lower=1,upper=N> nn[N_long];  // participant for observation n
  int<lower=1,upper=J> jj[N_long];  // question for observation n
  int<lower=0,upper=Ncateg_max> y[N_long];   // correctness for observation n
  int<lower=1> D;        // number of first-order dimensions
  int<lower=1> nDelta;        // total number of delta parameters
  int<lower=1> L;        // number of non-zero loadings (alpha parameters)
  int<lower=1> alpha_dstart[D]; // alpha start index for each dimension
  int<lower=1> alpha_dend[D];   // alpha end index for each dimension
  int<lower=1,upper=D> lambda_ind[N_long]; // which 1st order dim does each obs. belong to? 
  real weights[N_long]; // weights for each observation
}
parameters {
  matrix[N, D] theta_resid;        // residual 1st order ability
  vector[N] theta_g;               // general ability
  vector[nDelta] delta_l;          // difficulty 
  vector<lower=0>[L] alpha_l;      // distrimination over multiple dimensions
  vector[D] lambda;
}
transformed parameters {
  matrix[D, J] alpha; // connstrain the upper traingular elements to zero 
  vector[Ncateg_max-1] delta_trans[J]; // Make excess categories infinite

  {
    for(d in 1:D) {
      for(j in 1:J) {
        alpha[d, j] = 0;
      }
    }
  }
  {
    int aindex = 0;
    int a_lower = 0;
    int a_upper = 0;
    for(d in 1:D) {
      a_lower = alpha_dstart[d];
      a_upper = alpha_dend[d];
      for(i in a_lower:a_upper) {
        aindex += 1;
        alpha[d, i] = alpha_l[aindex];
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
  to_vector(theta_resid) ~ normal(0, 1);
  theta_g ~ normal(0, 1);
  alpha_l ~ lognormal(0, 0.3);
  delta_l ~ normal(0, 1);
  lambda ~ normal(0, 5);
  {
    vector[N_long] nu;
    for (i in 1:N_long) {
      nu[i] = ((lambda[lambda_ind[i]]*theta_g[nn[i]] + theta_resid[nn[i], lambda_ind[i]])*alpha[lambda_ind[i], jj[i]]);
      target += ordered_logistic_lpmf(y[i] | nu[i], delta_trans[jj[i]]) * weights[i];
    }
  }
}
