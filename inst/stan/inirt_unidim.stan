// Unidimensional IRT with latent regression (inirt package)
// Author: Michael Kleinsasser
// Description:
// Stan program meant to be used by the inirt::inirt() R function
// Example (test if it compiles to c++):
// mod = rstan::stan_model(file = "./inst/stan/inirt_unidim.stan", verbose = TRUE)
data {
  int<lower=1> N;              // number of participants
  int<lower=1> J;              // number of questions
  int<lower=0> K;              // number of regression parameters
  int<lower=0,upper=1> any_rand; // any random effects?
  int<lower=0,upper=1> any_rand_ind; // any independent random effects?
  int<lower=0,upper=1> any_rand_cor; // any correlated random effects?
  int<lower=0,upper=1> any_rand_ind_a; // any independent random effects?
  int<lower=0,upper=1> any_rand_cor_a; // any correlated random effects?
  int<lower=0,upper=1> any_rand_ind_d; // any independent random effects?
  int<lower=0,upper=1> any_rand_cor_d; // any correlated random effects?
  int<lower=2> Ncateg_max;         // Max Number of categories in ordered logistic responses
  int<lower=2,upper=Ncateg_max> Ncategi[J]; // Number of categories for each item
  int<lower=1> N_long;              // number of long observations
  int<lower=1,upper=N> nn[N_long];  // participant for observation n
  int<lower=1,upper=J> jj[N_long];  // question for observation n
  int<lower=0,upper=Ncateg_max> y[N_long];   // correctness for observation n
  matrix[N, K] x;   // fixed effect design matrix for observation n
  int<lower=1> D;        // number of first-order dimensions
  int<lower=1> nDelta;        // total number of delta parameters
  int<lower=1> L;        // number of non-zero loadings
  int<lower=0,upper=1> has_treg;  // do theta regression?
  int<lower=1> beta_dstart[has_treg ? D : 0]; // beta start index for each dimension
  int<lower=1> beta_dend[has_treg ? D : 0];   // beta end index for each dimension
  int<lower=1> zeta_dstart[any_rand_ind ? D : 0]; // zeta start index for each dimension
  int<lower=1> zeta_dend[any_rand_ind ? D : 0];   // zeta end index for each dimension
  real weights[N_long]; // weights for each observation
  matrix[N_long, K] x_miss;    // missing x index matrix (1 if missing, 0 else)
  // int reg_miss[N, K];       // id value of missing x within a matrix, 0 else
  // int<lower=0,upper=1> x_in_row_is_missing[N_long]; // any missing x's in given row? for efficiency
  int<lower=1> nDelta_r;            // number of delta structural regression parameters
  int<lower=1> nAlpha_r;            // number of alpha structural regression parameters
  matrix[N, nDelta_r] d_design;     // delta structural design matrix
  matrix[N, nAlpha_r] a_design;     // alpha structural design matrix
  
  int<lower=0> Lzeta;        // Number of uncorrelated random eff. parms
  int<lower=0> Laeta;
  int<lower=0> Ldeta;
  
  int<lower=0> u_Lzeta_cor;      // number of corr ranef par vectors
  int<lower=0> l_Lzeta_cor;      // length of corr ranef par vectors
  int<lower=0> u_Laeta_cor;
  int<lower=0> l_Laeta_cor; 
  int<lower=0> u_Ldeta_cor;
  int<lower=0> l_Ldeta_cor; 
  
  int<lower=0> Lzeta_cor;    // total number of correlated random eff. parms
  int<lower=0> Laeta_cor;
  int<lower=0> Ldeta_cor;
  matrix[N, Lzeta] z;   // design matrix for the uncorrelated random effects
  matrix[N, Laeta] ar;
  matrix[N, Ldeta] dr;

  int<lower=0> Lzeta_sd;     // number of sd pars
  int<lower=0> zeta_sd_ind[Lzeta]; // sd index for each column of z
  int<lower=0> cor_z_item_ind[Lzeta_cor]; // item index for each column of z for correlated random effs.
  int<lower=0> cor_z_item_elem_ind[Lzeta_cor]; // element within item index for each column of z for correlated random effs.

  int<lower=0> Laeta_sd;
  int<lower=0> alindex[Laeta];
  int<lower=0> aeta_sd_ind[Laeta];
  int<lower=0> cor_a_item_ind[Laeta_cor];
  int<lower=0> cor_a_item_elem_ind[Laeta_cor];

  int<lower=0> Ldeta_sd;
  int<lower=0> dlindex[Ldeta];
  int<lower=0> deta_sd_ind[Ldeta];
  int<lower=0> cor_d_item_ind[Ldeta_cor];
  int<lower=0> cor_d_item_elem_ind[Ldeta_cor];

  matrix[N, Lzeta_cor] z_c;   // design matrix for the correlated random effects
  matrix[N, Laeta_cor] a_c;
  matrix[N, Ldeta_cor] d_c;
}
transformed data {
  vector[l_Lzeta_cor] zeros_Lzeta_cor;
  zeros_Lzeta_cor = rep_vector(0, l_Lzeta_cor);

  vector[l_Laeta_cor] zeros_Laeta_cor;
  zeros_Laeta_cor = rep_vector(0, l_Laeta_cor);

  vector[l_Ldeta_cor] zeros_Ldeta_cor;
  zeros_Ldeta_cor = rep_vector(0, l_Ldeta_cor);
}
parameters {
  matrix[N, D] theta;              // ability

  vector[nDelta] delta_l;          // difficulty
  vector[nDelta_r] delta_r_l;      // structural regression, delta
  vector[L] alpha_l;      // distrimination over multiple dimensions
  // real<lower=0> sigma_alpha;
  vector[nAlpha_r] alpha_r_l;      // structural regression, alpha
  vector[K] beta_l;            // regression parameters for each dimension

  vector[Lzeta] zeta_l;          // random regression pars
  vector<lower=0>[Lzeta_sd] zeta_l_sd;          // random regression pars
  vector[Laeta] aeta_l;
  vector<lower=0>[Laeta_sd] aeta_l_sd; 
  vector[Ldeta] deta_l;
  vector<lower=0>[Ldeta_sd] deta_l_sd; 

  corr_matrix[l_Lzeta_cor] Omega;        // prior correlation
  vector<lower=0>[l_Lzeta_cor] tau;              // prior scale
  corr_matrix[l_Laeta_cor] Omega_a;     
  vector<lower=0>[l_Laeta_cor] tau_a;
  corr_matrix[l_Ldeta_cor] Omega_d;     
  vector<lower=0>[l_Ldeta_cor] tau_d;

  vector[l_Lzeta_cor] zeta_c[u_Lzeta_cor];          // random regression pars
  vector[l_Laeta_cor] aeta_c[u_Laeta_cor];
  vector[l_Ldeta_cor] deta_c[u_Ldeta_cor];

  // real g_phi; // gamma parameter for the alpha ~ item parameter regression model
}
transformed parameters {
  matrix[D, J] alpha;                  // connstrain the upper traingular elements to zero 
  matrix[K, D] beta;               // organize regression parameters into a matrix            beta and zeta could potentially be eliminated??
  matrix[Lzeta, D] zeta;               // organize ranef regression parameters into a matrix
  vector[Ncateg_max-1] delta_trans[J]; // Make excess categories infinite
  vector[N_long] db;
  vector[N_long] ab;
  vector[N_long] xb;
  // vector[L] g_mu;
  // vector[L] g_alpha;
  // vector[L] g_beta;

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
    if(has_treg) {
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

  {
    if(any_rand_ind) {
      int zindex = 0;
      int z_lower = 0;
      int z_upper = 0;
      for(d in 1:D) {
        z_lower = zeta_dstart[d];
        z_upper = zeta_dend[d];
        for(i in z_lower:z_upper) {
          zindex = zindex + 1;
          zeta[i, d] = zeta_l[zindex]*zeta_l_sd[zeta_sd_ind[i]];
        }
      }
    }
  }

  // ensures delta parameters satisfy requirements of the ordered logistic 
  // distribution according to STAN's definition
  // Basically, sort and make NA values inf
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

  {
    for(i in 1:N_long) {
      db[i] = d_design[nn[i], ]*delta_r_l;
      if(any_rand_ind_d) {
        for(k in 1:Ldeta) {
          db[i] += dr[nn[i], k] * deta_l[dlindex[k]]*deta_l_sd[deta_sd_ind[k]];
        }
      }
      if(any_rand_cor_d) {
        for(k in 1:Ldeta_cor) {
          db[i] += d_c[nn[i], k] * deta_c[cor_d_item_ind[k]][cor_d_item_elem_ind[k]];
        }
      }
    }
  }

  {
    for(i in 1:N_long) {
      ab[i] = a_design[nn[i], ]*alpha_r_l;
      if(any_rand_ind_a) {
        for(k in 1:Laeta) {
          ab[i] += ar[nn[i], k] * aeta_l[alindex[k]]*aeta_l_sd[aeta_sd_ind[k]];
        }
      }
      if(any_rand_cor_a) {
        for(k in 1:Laeta_cor) {
          ab[i] += a_c[nn[i], k] * aeta_c[cor_a_item_ind[k]][cor_a_item_elem_ind[k]];
        }
      }
    }
  }

  {
    for(i in 1:N_long) {
      if(has_treg) {
        for(k in 1:K) {
          if(x_miss[i, k] == 0) {
            xb[i] += x[nn[i], k] * beta[k,1];
          }
        }
      } else {
        xb[i] = 0.0;
      }
      if(any_rand_cor) {
        for(k in 1:Lzeta_cor) {
          xb[i] += z_c[nn[i], k] * zeta_c[cor_z_item_ind[k]][cor_z_item_elem_ind[k]];
        }
      }
      if(any_rand_ind) {
        for(k in 1:Lzeta) {
          xb[i] += z[nn[i], k] * zeta[k,1];
        }
      }
    }
  }

  // Transformed parameters for the regression on alpha parameters
  // alpha ~ x1 + x2 + ... Part of the item covariate portion
  // {
  //   g_mu = exp(ab);
  //   g_alpha = g_mu .* g_mu / g_phi;
  //   g_beta = g_mu / g_phi;
  // }
  
}
model {
  to_vector(theta) ~ normal(0, 1);
  // alpha_l ~ lognormal(0, 0.3);
  // alpha_l ~ lognormal(0, sigma_alpha); //cauchy(0, 5);
  // sigma_alpha ~ cauchy(0, 5);
  // alpha_l ~ gamma(g_alpha, g_beta);
  alpha_l ~ normal(0, 0.5);
  alpha_r_l ~ normal(0, 0.5); //cauchy(0, 5);
  // sigma_alpha ~ cauchy(0, 5);
  
  delta_l ~ normal(0, 1);
  delta_r_l ~ normal(0, 1);

  if(has_treg) {
    beta_l ~ normal(0, 5);
  }

  if(any_rand_ind) {
    zeta_l  ~ normal(0, 1);
    zeta_l_sd ~ cauchy(0, 5);
  }

  if(any_rand_cor) {
    tau ~ cauchy(0, 2.5);
    Omega ~ lkj_corr(1); //l_Lzeta_cor
    for(i in 1:u_Lzeta_cor) {
      zeta_c[i] ~ multi_normal(zeros_Lzeta_cor, quad_form_diag(Omega, tau));
    }
  }

  if(any_rand_ind_a) {
    aeta_l  ~ normal(0, 1);
    aeta_l_sd ~ cauchy(0, 5); // may need to hammer down on this?
  }

  if(any_rand_cor_a) {
    tau_a ~ cauchy(0, 2.5); // may need to hammer down on this?
    Omega_a ~ lkj_corr(1);
    for(i in 1:u_Laeta_cor) {
      aeta_c[i] ~ multi_normal(zeros_Laeta_cor, quad_form_diag(Omega_a, tau_a));
    }
  }

  if(any_rand_ind_d) {
    deta_l  ~ normal(0, 1);
    deta_l_sd ~ cauchy(0, 5);
  }

  if(any_rand_cor_d) {
    tau_d ~ cauchy(0, 2.5);
    Omega_d ~ lkj_corr(1);
    for(i in 1:u_Ldeta_cor) {
      deta_c[i] ~ multi_normal(zeros_Ldeta_cor, quad_form_diag(Omega_d, tau_d));
    }
  }
  
  {
    vector[N_long] nu;
    for (i in 1:N_long) {
      // likelihood for the model
      nu[i] = (theta[nn[i], ] + xb[i])*(exp(col(alpha, jj[i]) + ab[i]));
      target += ordered_logistic_lpmf(y[i] | nu[i], delta_trans[jj[i]] + db[i]) * weights[i];
    }
  }
}












// real xb = 0.0;
      // real ab = a_design[nn[i], ]*alpha_r_l;
      // real db = d_design[nn[i], ]*delta_r_l;

      // for(k in 1:K) {
      //   if(x_miss[i, k] == 0) {
      //     xb += x[nn[i], k] * beta[k,1];
      //   }
      // }
      // // if(any_rand) { // handle random effects
      // if(any_rand_cor) {
      //   for(k in 1:Lzeta_cor) {
      //     xb += z_c[nn[i], k] * zeta_c[cor_z_item_ind[k]][cor_z_item_elem_ind[k]];
      //   }
      // }
      // if(any_rand_ind) {
      //   for(k in 1:Lzeta) {
      //     xb += z[nn[i], k] * zeta[k,1];
      //   }
      // }
      // // }
      // if(any_rand_ind_a) {
      //   for(k in 1:Laeta) {
      //     ab += ar[nn[i], k] * aeta_l[alindex[k]]*aeta_l_sd[aeta_sd_ind[k]];
      //   }
      // }
      // if(any_rand_cor_a) {
      //   for(k in 1:Laeta_cor) {
      //     ab += a_c[nn[i], k] * aeta_c[cor_a_item_ind[k]][cor_a_item_elem_ind[k]];
      //   }
      // }
      // if(any_rand_ind_d) {
      //   for(k in 1:Ldeta) {
      //     db += dr[nn[i], k] * deta_l[dlindex[k]]*deta_l_sd[deta_sd_ind[k]];
      //   }
      // }
      // if(any_rand_cor_d) {
      //   for(k in 1:Ldeta_cor) {
      //     db += d_c[nn[i], k] * deta_c[cor_d_item_ind[k]][cor_d_item_elem_ind[k]];
      //   }
      // }
