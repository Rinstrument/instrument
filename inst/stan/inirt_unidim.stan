data {
  int<lower=1> N;              // number of participants
  int<lower=1> J;              // number of questions
  int<lower=1> K;              // number of regression parameters
  // int<lower=0,upper=1> any_rand; // any random effects?
  // int<lower=1> Krand;              // number of random effect pars (assuming there are any) Krand is part of K, i.e., K = Knon_rand + Krand
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
  int<lower=0,upper=1> has_treg;  // do theta regression?
  // int<lower=1> Lzeta;    // number of rand eff reg. parameters
  int<lower=1> beta_dstart[D]; // beta start index for each dimension
  int<lower=1> beta_dend[D];   // beta end index for each dimension
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
  matrix[N, D] theta;              // ability
  vector[nDelta] delta_l;          // difficulty
  vector[nDelta_r] delta_r_l;      // structural regression, delta
  vector<lower=0>[L] alpha_l;      // distrimination over multiple dimensions
  vector[nAlpha_r] alpha_r_l;      // structural regression, alpha
  vector[K] beta_l;            // regression parameters for each dimension
  // vector[Lzeta] zeta_l;          // random regression pars
}
transformed parameters {
  matrix[D, J] alpha;                  // connstrain the upper traingular elements to zero 
  matrix[K, D] beta;               // organize regression parameters into a matrix
  vector[Ncateg_max-1] delta_trans[J]; // Make excess categories infinite

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
  alpha_r_l ~ normal(0, 1);
  delta_l ~ normal(0, 1);
  delta_r_l ~ normal(0, 1);

  if(has_treg) {
    beta_l ~ normal(0, 5);
  }
  
  {
    vector[N_long] nu;
    for (i in 1:N_long) {
      if(has_treg) {
        real xb = 0.0;
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
        nu[i] = (theta[nn[i], ] + xb)*(col(alpha, jj[i]) + a_design[nn[i], ]*alpha_r_l);
      } else {
        nu[i] = (theta[nn[i], ])*(col(alpha, jj[i]) + a_design[nn[i], ]*alpha_r_l);
      }
      target += ordered_logistic_lpmf(y[i] | nu[i], delta_trans[jj[i]] + d_design[nn[i], ]*delta_r_l) * weights[i];
    }
  }
}
