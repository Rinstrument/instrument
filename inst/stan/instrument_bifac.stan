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
  int<lower=1> L;        // number of non-zero loadings (alpha parameters)
  int<lower=0,upper=1> has_treg;  // do theta regression?
  int<lower=1> alpha_dstart[D]; // alpha start index for each dimension
  int<lower=1> alpha_dend[D];   // alpha end index for each dimension
  int<lower=1,upper=D> lambda_ind[N_long]; // which 1st order dim does each obs. belong to? 
  int<lower=1> beta_dstart; // beta start index for each dimension
  int<lower=1> beta_dend;   // beta end index for each dimension
  real weights[N_long]; // weights for each observation
  matrix[N_long, K] x_miss;    // missing x index matrix (1 if missing, 0 else)
  int reg_miss[N, K];       // id value of missing x within a matrix, 0 else
  int<lower=0,upper=1> x_in_row_is_missing[N_long]; // any missing x's in given row? for efficiency
  int<lower=1> nDelta_r;            // number of delta structural regression parameters
  int<lower=1> nAlpha_r;            // number of alpha structural regression parameters
  matrix[N, nDelta_r] d_design;     // delta structural design matrix
  matrix[N, nAlpha_r] a_design;     // alpha structural design matrix
}
parameters {
  matrix[N, D] theta_resid;        // residual 1st order ability
  vector[N] theta_g;               // general ability
  vector[nDelta] delta_l;          // difficulty 
  vector[nDelta_r] delta_r_l;      // structural regression, delta
  vector<lower=0>[L] alpha_l;      // distrimination over multiple dimensions
  vector[nAlpha_r] alpha_r_l;      // structural regression, alpha
  vector[K] beta_l;            // regression parameters for each dimension
}
transformed parameters {
  matrix[D, J] alpha; // connstrain the upper traingular elements to zero 
  matrix[K, 1] beta; // organize regression parameters into a matrix
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
    for(i in 1:K) {
      beta[i, 1] = 0;
    }
  }
  {
    int bindex = 0;
    int b_lower = 0;
    int b_upper = 0;
    b_lower = beta_dstart;
    b_upper = beta_dend;
    for(i in b_lower:b_upper) {
      bindex = bindex + 1;
      beta[i, 1] = beta_l[bindex];
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

  if(has_treg) {
    beta_l ~ normal(0, 5);
  }

  {
    vector[N_long] nu;
    for (i in 1:N_long) {
      real xb = 0.0;
      if(has_treg) {
        if(x_in_row_is_missing[i]) {
          for(k in 1:K) {
            if(x_miss[i, k]) {
              xb += 0.0;
            } else {
              xb += x[nn[i], k] * beta[k,1];
            }
          }
        } else {
          xb = x[nn[i], ] * beta[,1];
        }
        nu[i] = (theta_g[nn[i]] + xb + theta_resid[nn[i], lambda_ind[i]])*(alpha[lambda_ind[i], jj[i]] + a_design[nn[i], ]*alpha_r_l);
      } else {
        nu[i] = (theta_g[nn[i]] + theta_resid[nn[i], lambda_ind[i]])*(alpha[lambda_ind[i], jj[i]] + a_design[nn[i], ]*alpha_r_l);
      }
      target += ordered_logistic_lpmf(y[i] | nu[i], delta_trans[jj[i]] + d_design[nn[i], ]*delta_r_l) * weights[i];
    }
  }
}
