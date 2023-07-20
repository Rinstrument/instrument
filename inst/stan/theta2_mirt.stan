// Unidimensional IRT with latent regression (theta2 package)
// Author: Michael Kleinsasser
// Description:
// Stan program meant to be used by the theta2::theta2() R function
// Example (test if it compiles to c++):
// mod = rstan::stan_model(file = "./inst/stan/theta2_mirt.stan", verbose = TRUE)
// rstan::stanc(file = "./inst/stan/theta2_mirt.stan", verbose = TRUE)
// https://github.com/henrixapp/muq2/blob/35d366b07cf1929c03e1ac8b5e6f1f355e12a760/external/include/stan/prob/distributions/univariate/discrete/ordered_logistic.hpp
functions {
  real ordered_logistic_log_irt_vec(array[] int y, vector nu, matrix cut, 
    vector eta, array[] int K, int nlong, array[] int itype) {
    real val = 0.0;
    int K_i = 0;
    int y_i = 0;
    real nu_i = 0.0;
    real eta_i = 0.0;
    for(i in 1:nlong) {
      K_i = K[i];
      y_i = y[i];
      nu_i = nu[i];
      if(itype[i] < 3) {
        if (y_i == 1) {
          val += -log1p_exp(nu_i - cut[i, 1]);
          //val += log(inv_logit(cut[i, 1] - nu_i));
        } else if(y_i == K_i) {
          val += -log1p_exp(cut[i, K_i-1] - nu_i);
          //val += log(inv_logit(nu_i - cut[i, K_i-1]));
        } else {
          val += log_inv_logit_diff(cut[i, y_i] - nu_i, cut[i, y_i-1] - nu_i);
        }
      } else {
        eta_i = eta[i];
        if (y_i == 1) {
        //  val += log(eta_i + (1.0 - eta_i)*inv_logit(cut[i, 1] - nu_i));
         val += log(eta_i + (1.0 - eta_i)*(1.0 - inv_logit(nu_i - cut[i, 1]))); //nu_i - cut[i, 1]
         //val += log(eta_i + (1.0 - eta_i)*(inv_logit(cut[i, 1] - nu_i)));
        } else if(y_i == K_i) {
         val += log(eta_i + (1.0 - eta_i)*inv_logit(nu_i - cut[i, K_i-1])); //nu_i - cut[i, K_i-1]
         //val += log(eta_i + (1.0 - eta_i)*(inv_logit(nu_i - cut[i, K_i-1])));
        } else {
          //val += log(eta_i + (1.0 - eta_i)*(inv_logit(cut[i, y_i] - nu_i))) - log(eta_i + (1-eta_i)*(inv_logit(cut[i, y_i-1] - nu_i)));
          //val += log(    (1.0 - eta_i)*(inv_logit(cut[i, y_i] - nu_i)) -   (1.0-eta_i)*(inv_logit(cut[i, y_i-1] - nu_i))       );
          val += log(  (1.0 - eta_i)*(inv_logit(cut[i, y_i-1] - nu_i))    -    (1.0 - eta_i)*(inv_logit(cut[i, y_i] - nu_i))     );
        }
      }
    }
    return val;
  }
}

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
  int<lower=1> Ncategi_jj[N_long];
  int<lower=1,upper=N> nn[N_long];  // participant for observation n
  int<lower=1,upper=J> jj[N_long];  // question for observation n
  int<lower=0,upper=Ncateg_max> y[N_long];   // correctness for observation n
  int<lower=1,upper=3> itype[N_long];      // item type for each item (1pl=1, 2pl=2, 3pl=3)
  int<lower=0,upper=1> any_eta3pl;       // any eta parameters?
  int<lower=0> nEta3pl; // number of eta parameters
  int<lower=0> find_eta3pl[N_long];      // find the correct eta parameter for long format of data
  matrix[N, K] x;   // fixed effect design matrix for observation n
  matrix[N_long, K] xLong;   // fixed effect design matrix for observation n (long format)
  int<lower=1> D;        // number of first-order dimensions
  int<lower=0> DAlpha;  // Copy of D except that it is zero if itype == "1pl" - no alpha parameters estimated in this case
  // int<lower=0> exploratory;
  int<lower=0> alpha_dstart[D];
  int<lower=0> alpha_dend[D];
  int<lower=1> nDelta;        // total number of delta parameters
  // int<lower=1,upper=D> lambda_ind[N_long]; // which 1st order dim does each obs. belong to? 
  int<lower=0> L;        // number of non-zero loadings
  int<lower=0,upper=1> has_treg;  // do theta regression?
  int<lower=0> beta_dstart[has_treg ? D : 0]; // beta start index for each dimension
  int<lower=0> beta_dend[has_treg ? D : 0];   // beta end index for each dimension
  int<lower=0> zeta_dstart[any_rand_ind ? D : 0]; // zeta start index for each dimension
  int<lower=0> zeta_dend[any_rand_ind ? D : 0];   // zeta end index for each dimension
  real fweights[N_long]; // weights for each observation
  //matrix[N, K] x_miss;    // missing x index matrix (1 if missing, 0 else)
  // int reg_miss[N, K];       // id value of missing x within a matrix, 0 else
  // int<lower=0,upper=1> x_in_row_is_missing[N_long]; // any missing x's in given row? for efficiency
  int<lower=0> nDelta_r;            // number of delta structural regression parameters
  int<lower=0> nAlpha_r;            // number of alpha structural regression parameters
  int<lower=0> LMean;        // want at least intercept for alpha pars?
  int<lower=0> deltaMean;    // want at least intercept for delta pars?
  matrix[N, nDelta_r] d_design;     // delta structural design matrix
  matrix[N, nAlpha_r] a_design;     // alpha structural design matrix

  int<lower=0> which_dim_fixed_reg[D];
  int<lower=0> which_dim_fixed_reg_sort[D];
  
  int<lower=0> which_dim_ind_reg[D];
  int<lower=0> which_dim_ind_reg_sort[D];
  // whether to use extra slots for independend ranef matrices
  int<lower=0,upper=1> rand_ind_g1 ;
  int<lower=0,upper=1> rand_ind_g2 ;
  int<lower=0,upper=1> rand_ind_g3 ;
  int<lower=0,upper=1> rand_ind_g4 ;
  int<lower=0,upper=1> rand_ind_g5 ;
  int<lower=0,upper=1> rand_ind_g6 ;
  int<lower=0,upper=1> rand_ind_g7 ;
  int<lower=0,upper=1> rand_ind_g8 ;
  int<lower=0,upper=1> rand_ind_g9 ;
  int<lower=0,upper=1> rand_ind_g10;
  int<lower=0,upper=1> rand_ind_g11;
  int<lower=0,upper=1> rand_ind_g12;
  int<lower=0,upper=1> rand_ind_g13;
  int<lower=0,upper=1> rand_ind_g14;
  int<lower=0,upper=1> rand_ind_g15;
  int<lower=0,upper=1> rand_ind_g16;
  int<lower=0,upper=1> rand_ind_g17;
  int<lower=0,upper=1> rand_ind_g18;
  int<lower=0,upper=1> rand_ind_g19;
  int<lower=0,upper=1> rand_ind_g20;
  int<lower=0,upper=1> rand_ind_g21;
  int<lower=0,upper=1> rand_ind_g22;
  int<lower=0,upper=1> rand_ind_g23;
  int<lower=0,upper=1> rand_ind_g24;
  int<lower=0,upper=1> rand_ind_g25;
  int<lower=0,upper=1> rand_ind_g26;
  int<lower=0,upper=1> rand_ind_g27;
  int<lower=0,upper=1> rand_ind_g28;
  int<lower=0,upper=1> rand_ind_g29;
  int<lower=0,upper=1> rand_ind_g30;
  int<lower=0,upper=1> rand_ind_g31;

  int<lower=0> Lzeta   ;        // Number of uncorrelated random eff. parms
  int<lower=0> Lzeta_2 ;
  int<lower=0> Lzeta_3 ;
  int<lower=0> Lzeta_4 ;
  int<lower=0> Lzeta_5 ;
  int<lower=0> Lzeta_6 ; 
  int<lower=0> Lzeta_7 ;
  int<lower=0> Lzeta_8 ;
  int<lower=0> Lzeta_9 ;
  int<lower=0> Lzeta_10;
  int<lower=0> Lzeta_11;
  int<lower=0> Lzeta_12;
  int<lower=0> Lzeta_13;
  int<lower=0> Lzeta_14;
  int<lower=0> Lzeta_15;
  int<lower=0> Lzeta_16;
  int<lower=0> Lzeta_17;
  int<lower=0> Lzeta_18;
  int<lower=0> Lzeta_19;
  int<lower=0> Lzeta_20;
  int<lower=0> Lzeta_21;
  int<lower=0> Lzeta_22;
  int<lower=0> Lzeta_23;
  int<lower=0> Lzeta_24;
  int<lower=0> Lzeta_25;
  int<lower=0> Lzeta_26;
  int<lower=0> Lzeta_27;
  int<lower=0> Lzeta_28;
  int<lower=0> Lzeta_29;
  int<lower=0> Lzeta_30;
  int<lower=0> Lzeta_31;
  int<lower=0> Lzeta_32;

  int<lower=0> Laeta;
  int<lower=0> Ldeta;

  int<lower=0> which_dim_cor_reg[D];
  int<lower=0,upper=1> rand_cor_g1;
  int<lower=0,upper=1> rand_cor_g2;

  int<lower=0> u_Lzeta_cor;      // number of corr ranef par vectors per dimension with regression
  int<lower=0> u_Lzeta_cor_2;
  int<lower=0> u_Lzeta_cor_3;

  int<lower=0> l_Lzeta_cor;      // length of corr ranef par vectors per dimension with regression
  int<lower=0> l_Lzeta_cor_2;
  int<lower=0> l_Lzeta_cor_3;

  int<lower=0> u_Laeta_cor;
  int<lower=0> l_Laeta_cor; 
  int<lower=0> u_Ldeta_cor;
  int<lower=0> l_Ldeta_cor; 
  
  int<lower=0> Lzeta_cor;

  int<lower=0> Lzeta_cor_2;
  int<lower=0> Lzeta_cor_3;
  // int<lower=0> Lzeta_cor[n_rand_cor_sets];    // total number of correlated random eff. parms

  int<lower=0> Laeta_cor;
  int<lower=0> Ldeta_cor;

  matrix[N, Lzeta] z;   // design matrix for the uncorrelated random effects
  matrix[N, Lzeta_2] z_2;
  matrix[N, Lzeta_3] z_3;

  matrix[N_long, Lzeta   ] zLong   ;
  matrix[N_long, Lzeta_2 ] zLong_2 ;
  matrix[N_long, Lzeta_3 ] zLong_3 ;
  matrix[N_long, Lzeta_4 ] zLong_4 ;
  matrix[N_long, Lzeta_5 ] zLong_5 ;
  matrix[N_long, Lzeta_6 ] zLong_6 ;
  matrix[N_long, Lzeta_7 ] zLong_7 ;
  matrix[N_long, Lzeta_8 ] zLong_8 ;
  matrix[N_long, Lzeta_9 ] zLong_9 ;
  matrix[N_long, Lzeta_10] zLong_10;
  matrix[N_long, Lzeta_11] zLong_11;
  matrix[N_long, Lzeta_12] zLong_12;
  matrix[N_long, Lzeta_13] zLong_13;
  matrix[N_long, Lzeta_14] zLong_14;
  matrix[N_long, Lzeta_15] zLong_15;
  matrix[N_long, Lzeta_16] zLong_16;
  matrix[N_long, Lzeta_17] zLong_17;
  matrix[N_long, Lzeta_18] zLong_18;
  matrix[N_long, Lzeta_19] zLong_19;
  matrix[N_long, Lzeta_20] zLong_20;
  matrix[N_long, Lzeta_21] zLong_21;
  matrix[N_long, Lzeta_22] zLong_22;
  matrix[N_long, Lzeta_23] zLong_23;
  matrix[N_long, Lzeta_24] zLong_24;
  matrix[N_long, Lzeta_25] zLong_25;
  matrix[N_long, Lzeta_26] zLong_26;
  matrix[N_long, Lzeta_27] zLong_27;
  matrix[N_long, Lzeta_28] zLong_28;
  matrix[N_long, Lzeta_29] zLong_29;
  matrix[N_long, Lzeta_30] zLong_30;
  matrix[N_long, Lzeta_31] zLong_31;
  matrix[N_long, Lzeta_32] zLong_32;

  matrix[N, Laeta] ar;
  matrix[N, Ldeta] dr;

  int<lower=0> Lzeta_sd   ;     // number of sd pars
  int<lower=0> Lzeta_sd_2 ;
  int<lower=0> Lzeta_sd_3 ;
  int<lower=0> Lzeta_sd_4 ;
  int<lower=0> Lzeta_sd_5 ;
  int<lower=0> Lzeta_sd_6 ;
  int<lower=0> Lzeta_sd_7 ;
  int<lower=0> Lzeta_sd_8 ;
  int<lower=0> Lzeta_sd_9 ;
  int<lower=0> Lzeta_sd_10;
  int<lower=0> Lzeta_sd_11;
  int<lower=0> Lzeta_sd_12;
  int<lower=0> Lzeta_sd_13;
  int<lower=0> Lzeta_sd_14;
  int<lower=0> Lzeta_sd_15;
  int<lower=0> Lzeta_sd_16;
  int<lower=0> Lzeta_sd_17;
  int<lower=0> Lzeta_sd_18;
  int<lower=0> Lzeta_sd_19;
  int<lower=0> Lzeta_sd_20;
  int<lower=0> Lzeta_sd_21;
  int<lower=0> Lzeta_sd_22;
  int<lower=0> Lzeta_sd_23;
  int<lower=0> Lzeta_sd_24;
  int<lower=0> Lzeta_sd_25;
  int<lower=0> Lzeta_sd_26;
  int<lower=0> Lzeta_sd_27;
  int<lower=0> Lzeta_sd_28;
  int<lower=0> Lzeta_sd_29;
  int<lower=0> Lzeta_sd_30;
  int<lower=0> Lzeta_sd_31;
  int<lower=0> Lzeta_sd_32;

  int<lower=0> zeta_sd_ind   [Lzeta   ]; // sd index for each column of z
  int<lower=0> zeta_sd_ind_2 [Lzeta_2 ];
  int<lower=0> zeta_sd_ind_3 [Lzeta_3 ];
  int<lower=0> zeta_sd_ind_4 [Lzeta_4 ];
  int<lower=0> zeta_sd_ind_5 [Lzeta_5 ];
  int<lower=0> zeta_sd_ind_6 [Lzeta_6 ];
  int<lower=0> zeta_sd_ind_7 [Lzeta_7 ];
  int<lower=0> zeta_sd_ind_8 [Lzeta_8 ];
  int<lower=0> zeta_sd_ind_9 [Lzeta_9 ];
  int<lower=0> zeta_sd_ind_10[Lzeta_10];
  int<lower=0> zeta_sd_ind_11[Lzeta_11];
  int<lower=0> zeta_sd_ind_12[Lzeta_12];
  int<lower=0> zeta_sd_ind_13[Lzeta_13];
  int<lower=0> zeta_sd_ind_14[Lzeta_14];
  int<lower=0> zeta_sd_ind_15[Lzeta_15];
  int<lower=0> zeta_sd_ind_16[Lzeta_16];
  int<lower=0> zeta_sd_ind_17[Lzeta_17];
  int<lower=0> zeta_sd_ind_18[Lzeta_18];
  int<lower=0> zeta_sd_ind_19[Lzeta_19];
  int<lower=0> zeta_sd_ind_20[Lzeta_20];
  int<lower=0> zeta_sd_ind_21[Lzeta_21];
  int<lower=0> zeta_sd_ind_22[Lzeta_22];
  int<lower=0> zeta_sd_ind_23[Lzeta_23];
  int<lower=0> zeta_sd_ind_24[Lzeta_24];
  int<lower=0> zeta_sd_ind_25[Lzeta_25];
  int<lower=0> zeta_sd_ind_26[Lzeta_26];
  int<lower=0> zeta_sd_ind_27[Lzeta_27];
  int<lower=0> zeta_sd_ind_28[Lzeta_28];
  int<lower=0> zeta_sd_ind_29[Lzeta_29];
  int<lower=0> zeta_sd_ind_30[Lzeta_30];
  int<lower=0> zeta_sd_ind_31[Lzeta_31];
  int<lower=0> zeta_sd_ind_32[Lzeta_32];

  matrix<lower=0>[Lzeta,    Lzeta_sd]    zeta_sd_ind_diag;
  matrix<lower=0>[Lzeta_2,  Lzeta_sd_2]  zeta_sd_ind_diag_2;
  matrix<lower=0>[Lzeta_3,  Lzeta_sd_3]  zeta_sd_ind_diag_3;
  matrix<lower=0>[Lzeta_4,  Lzeta_sd_4]  zeta_sd_ind_diag_4;
  matrix<lower=0>[Lzeta_5,  Lzeta_sd_5]  zeta_sd_ind_diag_5;
  matrix<lower=0>[Lzeta_6,  Lzeta_sd_6]  zeta_sd_ind_diag_6;
  matrix<lower=0>[Lzeta_7,  Lzeta_sd_7]  zeta_sd_ind_diag_7;
  matrix<lower=0>[Lzeta_8,  Lzeta_sd_8]  zeta_sd_ind_diag_8;
  matrix<lower=0>[Lzeta_9,  Lzeta_sd_9]  zeta_sd_ind_diag_9;
  matrix<lower=0>[Lzeta_10, Lzeta_sd_10] zeta_sd_ind_diag_10;
  matrix<lower=0>[Lzeta_11, Lzeta_sd_11] zeta_sd_ind_diag_11;
  matrix<lower=0>[Lzeta_12, Lzeta_sd_12] zeta_sd_ind_diag_12;
  matrix<lower=0>[Lzeta_13, Lzeta_sd_13] zeta_sd_ind_diag_13;
  matrix<lower=0>[Lzeta_14, Lzeta_sd_14] zeta_sd_ind_diag_14;
  matrix<lower=0>[Lzeta_15, Lzeta_sd_15] zeta_sd_ind_diag_15;
  matrix<lower=0>[Lzeta_16, Lzeta_sd_16] zeta_sd_ind_diag_16;
  matrix<lower=0>[Lzeta_17, Lzeta_sd_17] zeta_sd_ind_diag_17;
  matrix<lower=0>[Lzeta_18, Lzeta_sd_18] zeta_sd_ind_diag_18;
  matrix<lower=0>[Lzeta_19, Lzeta_sd_19] zeta_sd_ind_diag_19;
  matrix<lower=0>[Lzeta_20, Lzeta_sd_20] zeta_sd_ind_diag_20;
  matrix<lower=0>[Lzeta_21, Lzeta_sd_21] zeta_sd_ind_diag_21;
  matrix<lower=0>[Lzeta_22, Lzeta_sd_22] zeta_sd_ind_diag_22;
  matrix<lower=0>[Lzeta_23, Lzeta_sd_23] zeta_sd_ind_diag_23;
  matrix<lower=0>[Lzeta_24, Lzeta_sd_24] zeta_sd_ind_diag_24;
  matrix<lower=0>[Lzeta_25, Lzeta_sd_25] zeta_sd_ind_diag_25;
  matrix<lower=0>[Lzeta_26, Lzeta_sd_26] zeta_sd_ind_diag_26;
  matrix<lower=0>[Lzeta_27, Lzeta_sd_27] zeta_sd_ind_diag_27;
  matrix<lower=0>[Lzeta_28, Lzeta_sd_28] zeta_sd_ind_diag_28;
  matrix<lower=0>[Lzeta_29, Lzeta_sd_29] zeta_sd_ind_diag_29;
  matrix<lower=0>[Lzeta_30, Lzeta_sd_30] zeta_sd_ind_diag_30;
  matrix<lower=0>[Lzeta_31, Lzeta_sd_31] zeta_sd_ind_diag_31;
  matrix<lower=0>[Lzeta_32, Lzeta_sd_32] zeta_sd_ind_diag_32;
  //int<lower=1> zeta_sd_ind_ones[Lzeta_sd];

  // int<lower=0> Lzeta_sd_2;
  // int<lower=0> zeta_sd_ind_2[Lzeta_2];

  // int<lower=0> Lzeta_sd_3;
  // int<lower=0> zeta_sd_ind_3[Lzeta_3];

  int<lower=0> cor_z_item_ind[Lzeta_cor]; // item index for each column of z for correlated random effs.
  int<lower=0> cor_z_item_elem_ind[Lzeta_cor]; // element within item index for each column of z for correlated random effs.
  int<lower=0> cor_z_item_ind_2[Lzeta_cor_2];
  int<lower=0> cor_z_item_elem_ind_2[Lzeta_cor_2];
  int<lower=0> cor_z_item_ind_3[Lzeta_cor_3];
  int<lower=0> cor_z_item_elem_ind_3[Lzeta_cor_3];

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

  // extra slots of z_c
  matrix[N, Lzeta_cor_2] z_c_2;
  matrix[N, Lzeta_cor_3] z_c_3;

  matrix[N_long, Lzeta_cor] z_cLong; // long form equivalent of z_c
  matrix[N_long, Lzeta_cor_2] z_cLong_2;
  matrix[N_long, Lzeta_cor_3] z_cLong_3;

  matrix[N, Laeta_cor] a_c;
  matrix[N, Ldeta_cor] d_c;

}

transformed data {
  vector[l_Lzeta_cor] zeros_Lzeta_cor;
  zeros_Lzeta_cor = rep_vector(0, l_Lzeta_cor);

  vector[l_Lzeta_cor_2] zeros_Lzeta_cor_2;
  zeros_Lzeta_cor_2 = rep_vector(0, l_Lzeta_cor_2);
  vector[l_Lzeta_cor_3] zeros_Lzeta_cor_3;
  zeros_Lzeta_cor_3 = rep_vector(0, l_Lzeta_cor_3);

  vector[l_Laeta_cor] zeros_Laeta_cor;
  zeros_Laeta_cor = rep_vector(0, l_Laeta_cor);

  vector[l_Ldeta_cor] zeros_Ldeta_cor;
  zeros_Ldeta_cor = rep_vector(0, l_Ldeta_cor);
}

parameters {
  matrix[N, D] theta;              // general ability
  // real<lower=0,upper=1> sig_thetag_reg;
  // vector<lower=0>[D] sig_thetag_reg;          // residual uncertainty in thetag regression
  // matrix[N, D] theta_resid;        // residual 1st order ability

  // vector[D] lambda; // loadings for second order on first, e.g., theta1 = lambda * theta + error
  // <lower=-1,upper=1>

  vector[nDelta] delta_l;          // difficulty
  vector[nDelta_r] delta_r_l;      // structural regression, delta
  vector[L] alpha_l;      // distrimination over multiple dimensions
  vector[nAlpha_r] alpha_r_l;      // structural regression, alpha
  vector<lower=0,upper=1>[nEta3pl] eta3pl_l; //

  vector<lower=-1,upper=1>[K] beta_l;            // regression parameters for each dimension

  vector[Lzeta   ] zeta_l   ;          // random regression pars
  vector[Lzeta_2 ] zeta_l_2 ;
  vector[Lzeta_3 ] zeta_l_3 ;
  vector[Lzeta_4 ] zeta_l_4 ;
  vector[Lzeta_5 ] zeta_l_5 ;
  vector[Lzeta_6 ] zeta_l_6 ;
  vector[Lzeta_7 ] zeta_l_7 ;
  vector[Lzeta_8 ] zeta_l_8 ;
  vector[Lzeta_9 ] zeta_l_9 ;
  vector[Lzeta_10] zeta_l_10;
  vector[Lzeta_11] zeta_l_11;
  vector[Lzeta_12] zeta_l_12;
  vector[Lzeta_13] zeta_l_13;
  vector[Lzeta_14] zeta_l_14;
  vector[Lzeta_15] zeta_l_15;
  vector[Lzeta_16] zeta_l_16;
  vector[Lzeta_17] zeta_l_17;
  vector[Lzeta_18] zeta_l_18;
  vector[Lzeta_19] zeta_l_19;
  vector[Lzeta_20] zeta_l_20;
  vector[Lzeta_21] zeta_l_21;
  vector[Lzeta_22] zeta_l_22;
  vector[Lzeta_23] zeta_l_23;
  vector[Lzeta_24] zeta_l_24;
  vector[Lzeta_25] zeta_l_25;
  vector[Lzeta_26] zeta_l_26;
  vector[Lzeta_27] zeta_l_27;
  vector[Lzeta_28] zeta_l_28;
  vector[Lzeta_29] zeta_l_29;
  vector[Lzeta_30] zeta_l_30;
  vector[Lzeta_31] zeta_l_31;
  vector[Lzeta_32] zeta_l_32;
  
  vector<lower=0>[Lzeta_sd   ] zeta_l_sd;          // random regression pars
  vector<lower=0>[Lzeta_sd_2 ] zeta_l_sd_2;
  vector<lower=0>[Lzeta_sd_3 ] zeta_l_sd_3;
  vector<lower=0>[Lzeta_sd_4 ] zeta_l_sd_4;
  vector<lower=0>[Lzeta_sd_5 ] zeta_l_sd_5;
  vector<lower=0>[Lzeta_sd_6 ] zeta_l_sd_6;
  vector<lower=0>[Lzeta_sd_7 ] zeta_l_sd_7;
  vector<lower=0>[Lzeta_sd_8 ] zeta_l_sd_8;
  vector<lower=0>[Lzeta_sd_9 ] zeta_l_sd_9;
  vector<lower=0>[Lzeta_sd_10] zeta_l_sd_10;
  vector<lower=0>[Lzeta_sd_11] zeta_l_sd_11;
  vector<lower=0>[Lzeta_sd_12] zeta_l_sd_12; 
  vector<lower=0>[Lzeta_sd_13] zeta_l_sd_13;
  vector<lower=0>[Lzeta_sd_14] zeta_l_sd_14;
  vector<lower=0>[Lzeta_sd_15] zeta_l_sd_15;
  vector<lower=0>[Lzeta_sd_16] zeta_l_sd_16;
  vector<lower=0>[Lzeta_sd_17] zeta_l_sd_17;
  vector<lower=0>[Lzeta_sd_18] zeta_l_sd_18;
  vector<lower=0>[Lzeta_sd_19] zeta_l_sd_19;
  vector<lower=0>[Lzeta_sd_20] zeta_l_sd_20;
  vector<lower=0>[Lzeta_sd_21] zeta_l_sd_21;
  vector<lower=0>[Lzeta_sd_22] zeta_l_sd_22;
  vector<lower=0>[Lzeta_sd_23] zeta_l_sd_23;
  vector<lower=0>[Lzeta_sd_24] zeta_l_sd_24;
  vector<lower=0>[Lzeta_sd_25] zeta_l_sd_25;
  vector<lower=0>[Lzeta_sd_26] zeta_l_sd_26;
  vector<lower=0>[Lzeta_sd_27] zeta_l_sd_27;
  vector<lower=0>[Lzeta_sd_28] zeta_l_sd_28;
  vector<lower=0>[Lzeta_sd_29] zeta_l_sd_29;
  vector<lower=0>[Lzeta_sd_30] zeta_l_sd_30;
  vector<lower=0>[Lzeta_sd_31] zeta_l_sd_31;
  vector<lower=0>[Lzeta_sd_32] zeta_l_sd_32;
  
  // vector<lower=0>[Lzeta_sd_2] zeta_l_sd_2;

  // vector[Lzeta_3] zeta_l_3;
  // vector<lower=0>[Lzeta_sd_3] zeta_l_sd_3;

  vector[Laeta] aeta_l;
  vector<lower=0>[Laeta_sd] aeta_l_sd; 
  vector[Ldeta] deta_l;
  vector<lower=0>[Ldeta_sd] deta_l_sd; 

  corr_matrix[l_Lzeta_cor] Omega;        // prior correlation
  vector<lower=0>[l_Lzeta_cor] tau;              // prior scale

  // extra slots in the case of multiple separate regressions
  corr_matrix[l_Lzeta_cor_2] Omega_2;
  vector<lower=0>[l_Lzeta_cor_2] tau_2;
  corr_matrix[l_Lzeta_cor_3] Omega_3;
  vector<lower=0>[l_Lzeta_cor_3] tau_3;

  corr_matrix[l_Laeta_cor] Omega_a;
  vector<lower=0>[l_Laeta_cor] tau_a;
  corr_matrix[l_Ldeta_cor] Omega_d;
  vector<lower=0>[l_Ldeta_cor] tau_d;

  vector[l_Lzeta_cor] zeta_c[u_Lzeta_cor];          // random regression pars

  // extra slots in the case of multiple separate regressions
  vector[l_Lzeta_cor_2] zeta_c_2[u_Lzeta_cor_2];
  vector[l_Lzeta_cor_3] zeta_c_3[u_Lzeta_cor_3];

  vector[l_Laeta_cor] aeta_c[u_Laeta_cor];
  vector[l_Ldeta_cor] deta_c[u_Ldeta_cor];

}

transformed parameters {

  matrix[D*(DAlpha ? 1 : 0), J] alpha;                  // connstrain the upper traingular elements to zero 
  matrix[K, D] betat; // may generalize to [K,D]              // organize regression parameters into a matrix            beta and zeta could potentially be eliminated??

  matrix[Lzeta,   1]  zeta  ;               // organize ranef regression parameters into a matrix
  matrix[Lzeta_2, 1]  zeta_2;
  matrix[Lzeta_3, 1]  zeta_3;
  matrix[Lzeta_4, 1]  zeta_4;
  matrix[Lzeta_5, 1]  zeta_5;
  matrix[Lzeta_6, 1]  zeta_6;
  matrix[Lzeta_7, 1]  zeta_7;
  matrix[Lzeta_8, 1]  zeta_8;
  matrix[Lzeta_9, 1]  zeta_9;
  matrix[Lzeta_10, 1] zeta_10;
  matrix[Lzeta_11, 1] zeta_11;
  matrix[Lzeta_12, 1] zeta_12;
  matrix[Lzeta_13, 1] zeta_13;
  matrix[Lzeta_14, 1] zeta_14;
  matrix[Lzeta_15, 1] zeta_15;
  matrix[Lzeta_16, 1] zeta_16;
  matrix[Lzeta_17, 1] zeta_17;
  matrix[Lzeta_18, 1] zeta_18;
  matrix[Lzeta_19, 1] zeta_19;
  matrix[Lzeta_20, 1] zeta_20;
  matrix[Lzeta_21, 1] zeta_21;
  matrix[Lzeta_22, 1] zeta_22;
  matrix[Lzeta_23, 1] zeta_23;
  matrix[Lzeta_24, 1] zeta_24;
  matrix[Lzeta_25, 1] zeta_25;
  matrix[Lzeta_26, 1] zeta_26;
  matrix[Lzeta_27, 1] zeta_27;
  matrix[Lzeta_28, 1] zeta_28;
  matrix[Lzeta_29, 1] zeta_29;
  matrix[Lzeta_30, 1] zeta_30;
  matrix[Lzeta_31, 1] zeta_31;
  matrix[Lzeta_32, 1] zeta_32;



  matrix[J, Ncateg_max-1] delta; // Make excess categories infinite
  vector[N_long] db;
  vector[N_long*(DAlpha ? 1 : 0)] ab;
  matrix[N_long, D] xb;
  vector[N_long] nu;
  matrix[N_long, Ncateg_max-1] c;
  vector<lower=0,upper=1>[N_long] eta3pl;

  vector[Lzeta   ] zeta_l_sd_elong   ;
  vector[Lzeta_2 ] zeta_l_sd_elong_2 ;
  vector[Lzeta_3 ] zeta_l_sd_elong_3 ;
  vector[Lzeta_4 ] zeta_l_sd_elong_4 ;
  vector[Lzeta_5 ] zeta_l_sd_elong_5 ;
  vector[Lzeta_6 ] zeta_l_sd_elong_6 ;
  vector[Lzeta_7 ] zeta_l_sd_elong_7 ;
  vector[Lzeta_8 ] zeta_l_sd_elong_8 ;
  vector[Lzeta_9 ] zeta_l_sd_elong_9 ;
  vector[Lzeta_10] zeta_l_sd_elong_10;
  vector[Lzeta_11] zeta_l_sd_elong_11;
  vector[Lzeta_12] zeta_l_sd_elong_12;
  vector[Lzeta_13] zeta_l_sd_elong_13;
  vector[Lzeta_14] zeta_l_sd_elong_14;
  vector[Lzeta_15] zeta_l_sd_elong_15;
  vector[Lzeta_16] zeta_l_sd_elong_16;
  vector[Lzeta_17] zeta_l_sd_elong_17;
  vector[Lzeta_18] zeta_l_sd_elong_18;
  vector[Lzeta_19] zeta_l_sd_elong_19;
  vector[Lzeta_20] zeta_l_sd_elong_20;
  vector[Lzeta_21] zeta_l_sd_elong_21;
  vector[Lzeta_22] zeta_l_sd_elong_22;
  vector[Lzeta_23] zeta_l_sd_elong_23;
  vector[Lzeta_24] zeta_l_sd_elong_24;
  vector[Lzeta_25] zeta_l_sd_elong_25;
  vector[Lzeta_26] zeta_l_sd_elong_26;
  vector[Lzeta_27] zeta_l_sd_elong_27;
  vector[Lzeta_28] zeta_l_sd_elong_28;
  vector[Lzeta_29] zeta_l_sd_elong_29;
  vector[Lzeta_30] zeta_l_sd_elong_30;
  vector[Lzeta_31] zeta_l_sd_elong_31;
  vector[Lzeta_32] zeta_l_sd_elong_32;

  // new
  matrix[Lzeta_cor,     1]  zeta_cMat;
  matrix[Lzeta_cor_2,   1]  zeta_cMat_2;
  matrix[Lzeta_cor_3,   1]  zeta_cMat_3;

  {
    // elongate zeta_l_sd for fast dot product
    zeta_l_sd_elong    = (zeta_sd_ind_diag    * zeta_l_sd   ); //.* zeta_sd_ind_ones;
    zeta_l_sd_elong_2  = (zeta_sd_ind_diag_2  * zeta_l_sd_2 );
    zeta_l_sd_elong_3  = (zeta_sd_ind_diag_3  * zeta_l_sd_3 );
    zeta_l_sd_elong_4  = (zeta_sd_ind_diag_4  * zeta_l_sd_4 );
    zeta_l_sd_elong_5  = (zeta_sd_ind_diag_5  * zeta_l_sd_5 );
    zeta_l_sd_elong_6  = (zeta_sd_ind_diag_6  * zeta_l_sd_6 );
    zeta_l_sd_elong_7  = (zeta_sd_ind_diag_7  * zeta_l_sd_7 );
    zeta_l_sd_elong_8  = (zeta_sd_ind_diag_8  * zeta_l_sd_8 );
    zeta_l_sd_elong_9  = (zeta_sd_ind_diag_9  * zeta_l_sd_9 );
    zeta_l_sd_elong_10 = (zeta_sd_ind_diag_10 * zeta_l_sd_10);
    zeta_l_sd_elong_11 = (zeta_sd_ind_diag_11 * zeta_l_sd_11);
    zeta_l_sd_elong_12 = (zeta_sd_ind_diag_12 * zeta_l_sd_12);
    zeta_l_sd_elong_13 = (zeta_sd_ind_diag_13 * zeta_l_sd_13);
    zeta_l_sd_elong_14 = (zeta_sd_ind_diag_14 * zeta_l_sd_14);
    zeta_l_sd_elong_15 = (zeta_sd_ind_diag_15 * zeta_l_sd_15);
    zeta_l_sd_elong_16 = (zeta_sd_ind_diag_16 * zeta_l_sd_16);
    zeta_l_sd_elong_17 = (zeta_sd_ind_diag_17 * zeta_l_sd_17);
    zeta_l_sd_elong_18 = (zeta_sd_ind_diag_18 * zeta_l_sd_18);
    zeta_l_sd_elong_19 = (zeta_sd_ind_diag_19 * zeta_l_sd_19);
    zeta_l_sd_elong_20 = (zeta_sd_ind_diag_20 * zeta_l_sd_20);
    zeta_l_sd_elong_21 = (zeta_sd_ind_diag_21 * zeta_l_sd_21);
    zeta_l_sd_elong_22 = (zeta_sd_ind_diag_22 * zeta_l_sd_22);
    zeta_l_sd_elong_23 = (zeta_sd_ind_diag_23 * zeta_l_sd_23);
    zeta_l_sd_elong_24 = (zeta_sd_ind_diag_24 * zeta_l_sd_24);
    zeta_l_sd_elong_25 = (zeta_sd_ind_diag_25 * zeta_l_sd_25);
    zeta_l_sd_elong_26 = (zeta_sd_ind_diag_26 * zeta_l_sd_26);
    zeta_l_sd_elong_27 = (zeta_sd_ind_diag_27 * zeta_l_sd_27);
    zeta_l_sd_elong_28 = (zeta_sd_ind_diag_28 * zeta_l_sd_28);
    zeta_l_sd_elong_29 = (zeta_sd_ind_diag_29 * zeta_l_sd_29);
    zeta_l_sd_elong_30 = (zeta_sd_ind_diag_30 * zeta_l_sd_30);
    zeta_l_sd_elong_31 = (zeta_sd_ind_diag_31 * zeta_l_sd_31);
    zeta_l_sd_elong_32 = (zeta_sd_ind_diag_32 * zeta_l_sd_32);

    // zeta_l_sd_elong = zeta_sd_ind_diag * zeta_sd_ind_ones;
    // for(ze in 1:Lzeta) {
    //   zeta_l_sd_elong[ze] = zeta_l_sd[zeta_sd_ind[ze]];
    // }
  }

  {
    for(j in 1:J) {
      for(d in (j+1):D) {
        alpha[d, j] = negative_infinity();
      }
    }
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
      for(k in 1:K) {
        for(d in 1:D) {
          betat[k, d] = 0.0;
        }
      }
      int bindex = 0;
      int b_lower = 0;
      int b_upper = 0;
      for(d in 1:D) {
        b_lower = beta_dstart[d];
        b_upper = beta_dend[d];
        for(i in b_lower:b_upper) {
          bindex = bindex + 1;
          betat[i, d] = beta_l[bindex];
        }
      }
    }
  }

  {
    if(any_rand_ind) {
      // int zindex = 0;
      // int z_lower = 0;
      // int z_upper = 0;
      // for(d in 1:D) {
      //   z_lower = zeta_dstart[d];
      //   z_upper = zeta_dend[d];
      //   for(i in z_lower:z_upper) {
      //     zindex = zindex + 1;
      //     zeta[i, which_dim_ind_reg[d]] = zeta_l[zindex]*zeta_l_sd[zeta_sd_ind[i]];
      //   } // try to vectorize this zeta math with dot_product()
      // }

      // z_lower = zeta_dstart[1];
      // z_upper = zeta_dend[1];
      // for(i in z_lower:z_upper) {
        // zindex = zindex + 1;
                       zeta[,    1] = zeta_l    .* zeta_l_sd_elong   ;
      if(rand_ind_g1 ) zeta_2[,  1] = zeta_l_2  .* zeta_l_sd_elong_2 ;
      if(rand_ind_g2 ) zeta_3[,  1] = zeta_l_3  .* zeta_l_sd_elong_3 ;
      if(rand_ind_g3 ) zeta_4[,  1] = zeta_l_4  .* zeta_l_sd_elong_4 ;
      if(rand_ind_g4 ) zeta_5[,  1] = zeta_l_5  .* zeta_l_sd_elong_5 ;
      if(rand_ind_g5 ) zeta_6[,  1] = zeta_l_6  .* zeta_l_sd_elong_6 ;
      if(rand_ind_g6 ) zeta_7[,  1] = zeta_l_7  .* zeta_l_sd_elong_7 ;
      if(rand_ind_g7 ) zeta_8[,  1] = zeta_l_8  .* zeta_l_sd_elong_8 ;
      if(rand_ind_g8 ) zeta_9[,  1] = zeta_l_9  .* zeta_l_sd_elong_9 ;
      if(rand_ind_g9 ) zeta_10[, 1] = zeta_l_10 .* zeta_l_sd_elong_10;
      if(rand_ind_g10) zeta_11[, 1] = zeta_l_11 .* zeta_l_sd_elong_11;
      if(rand_ind_g11) zeta_12[, 1] = zeta_l_12 .* zeta_l_sd_elong_12;
      if(rand_ind_g12) zeta_13[, 1] = zeta_l_13 .* zeta_l_sd_elong_13;
      if(rand_ind_g13) zeta_14[, 1] = zeta_l_14 .* zeta_l_sd_elong_14;
      if(rand_ind_g14) zeta_15[, 1] = zeta_l_15 .* zeta_l_sd_elong_15;
      if(rand_ind_g15) zeta_16[, 1] = zeta_l_16 .* zeta_l_sd_elong_16;
      if(rand_ind_g16) zeta_17[, 1] = zeta_l_17 .* zeta_l_sd_elong_17;
      if(rand_ind_g17) zeta_18[, 1] = zeta_l_18 .* zeta_l_sd_elong_18;
      if(rand_ind_g18) zeta_19[, 1] = zeta_l_19 .* zeta_l_sd_elong_19;
      if(rand_ind_g19) zeta_20[, 1] = zeta_l_20 .* zeta_l_sd_elong_20;
      if(rand_ind_g20) zeta_21[, 1] = zeta_l_21 .* zeta_l_sd_elong_21;
      if(rand_ind_g21) zeta_22[, 1] = zeta_l_22 .* zeta_l_sd_elong_22;
      if(rand_ind_g22) zeta_23[, 1] = zeta_l_23 .* zeta_l_sd_elong_23;
      if(rand_ind_g23) zeta_24[, 1] = zeta_l_24 .* zeta_l_sd_elong_24;
      if(rand_ind_g24) zeta_25[, 1] = zeta_l_25 .* zeta_l_sd_elong_25;
      if(rand_ind_g25) zeta_26[, 1] = zeta_l_26 .* zeta_l_sd_elong_26;
      if(rand_ind_g26) zeta_27[, 1] = zeta_l_27 .* zeta_l_sd_elong_27;
      if(rand_ind_g27) zeta_28[, 1] = zeta_l_28 .* zeta_l_sd_elong_28;
      if(rand_ind_g28) zeta_29[, 1] = zeta_l_29 .* zeta_l_sd_elong_29;
      if(rand_ind_g29) zeta_30[, 1] = zeta_l_30 .* zeta_l_sd_elong_30;
      if(rand_ind_g30) zeta_31[, 1] = zeta_l_31 .* zeta_l_sd_elong_31;
      if(rand_ind_g31) zeta_32[, 1] = zeta_l_32 .* zeta_l_sd_elong_32;
      // }

      // if(rand_ind_g1) {
      //   zindex = 0;
      //   z_lower = zeta_dstart[2];
      //   z_upper = zeta_dend[2];
      //   for(i in z_lower:z_upper) {
      //     zindex = zindex + 1;
      //     zeta_2[i, 1] = zeta_l_2[zindex]*zeta_l_sd_2[zeta_sd_ind_2[i]];
      //   }
      // }

      // if(rand_ind_g2) {
      //   zindex = 0;
      //   z_lower = zeta_dstart[3];
      //   z_upper = zeta_dend[3];
      //   for(i in z_lower:z_upper) {
      //     zindex = zindex + 1;
      //     zeta_3[i, 1] = zeta_l_3[zindex]*zeta_l_sd_3[zeta_sd_ind_3[i]];
      //   }
      // }

    }
  }

  // new 
  if(any_rand_cor) {

    for(k in 1:Lzeta_cor) {
      zeta_cMat[k, 1] = zeta_c[cor_z_item_ind[k]][cor_z_item_elem_ind[k]];
    }

    if(rand_cor_g1) {
      for(k in 1:Lzeta_cor_2) {
        zeta_cMat_2[k, 1] = zeta_c_2[cor_z_item_ind_2[k]][cor_z_item_elem_ind_2[k]];
      }
    }

    if(rand_cor_g2) {
      for(k in 1:Lzeta_cor_3) {
        zeta_cMat_3[k, 1] = zeta_c_3[cor_z_item_ind_3[k]][cor_z_item_elem_ind_3[k]];
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
      for(i in 1:(Ncategi[j]-1)) { // faster to combine this into a single for loop?
        d_index = d_index + 1;
        delta[j, i] = ds_ind[i];
      }
      for(i in (Ncategi[j]):(Ncateg_max-1)) {
        delta[j, i] = 1e7;// + j + i;
        idx = idx + 1;
      }
    }
  }

  {
    for(i in 1:N_long) {
      if(deltaMean) {
        db[i] = d_design[nn[i], ]*delta_r_l;
      } else {
        db[i] = 0.0;
      }
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

      if(L) {
        if(LMean) {
          ab[i] = a_design[nn[i], ]*alpha_r_l;
        } else {
          ab[i] = 0.0;
        }
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

      if(has_treg) {
        for(d in 1:D) {
          xb[i, d] = 0.0;
          //for(k in 1:K) {  // vectorize this for loop with dot_product
          //}
          //xb[i, d] += dot_product(x[nn[i], ], betat[, d]);
          xb[i, d] += dot_product(xLong[i, ], betat[, d]); // can this be vectorized over dimension D?
        }

      } else {
        for(d in 1:D) {
          xb[i, d] = 0.0; // can this be replaces with rep_vec(0)?
        }
      }
      
      if(any_rand_cor) {
        // edit here ------------------------------
        // prev version:
        // for(k in 1:Lzeta_cor) {
        //   xb[i, which_dim_cor_reg[1]] += z_c[nn[i], k] * zeta_c[cor_z_item_ind[k]][cor_z_item_elem_ind[k]];
                      
        // }
        xb[i, which_dim_cor_reg[1]] += dot_product(   z_cLong[i, ]    ,   zeta_cMat[ , 1]    );
        // edit here ------------------------------
        // dot_product(   z_cLong[i, ]    ,   zeta_cMat[ , 1]    );
        //+= dot_product(zLong[i,      ], zeta[,     1]) ;

        if(rand_cor_g1) {
          xb[i, which_dim_cor_reg[2]] += dot_product(   z_cLong_2[i, ]    ,   zeta_cMat_2[ , 1]    );
        }

        if(rand_cor_g2) {
          xb[i, which_dim_cor_reg[3]] += dot_product(   z_cLong_3[i, ]    ,   zeta_cMat_3[ , 1]    );
        }
        
        // for(rg in 1:n_rand_cor_sets) {
        //   for(k in 1:lzeta_cor[rg]) {
        //     xb[i, rg] += z_c[nn[i], k] * zeta_c[cor_z_item_ind[k]][cor_z_item_elem_ind[k]];
        //   }
        // }

        // for(k in 1:Lzeta_cor) {
        //   for(d in 1:D) {
        //     xb[i, d] += z_c[nn[i], k] * zeta_c[cor_z_item_ind[k]][cor_z_item_elem_ind[k]];
        //   }
        // }
      }

      if(any_rand_ind) {
      
        // for(k in 1:Lzeta) {

        // PREVIOUS:
        // xb[i, which_dim_ind_reg_sort[1]] += dot_product(z[nn[i], ], zeta[, 1]);
        
        // }

        // if(rand_ind_g1) {
        //   // for(k in 1:Lzeta_2) {
        //     // xb[i, which_dim_ind_reg[2]] += z[nn[i], k] * zeta[k, 2];
        //   xb[i, which_dim_ind_reg_sort[2]] += dot_product(z_2[nn[i], ], zeta_2[, 1]);
        //   // }
        // }

        // if(rand_ind_g2) {
        //   // for(k in 1:Lzeta_3) {
        //     // xb[i, which_dim_ind_reg[3]] += z[nn[i], k] * zeta[k, 3];
        //   xb[i, which_dim_ind_reg_sort[3]] += dot_product(z_3[nn[i], ], zeta_3[, 1]);
        //   // }
        // }

                         xb[i, which_dim_ind_reg_sort[1]]  += dot_product(zLong[i,      ], zeta[,     1]) ;
        if(rand_ind_g1)  xb[i, which_dim_ind_reg_sort[2]]  += dot_product(zLong_2[i,    ], zeta_2[,   1]) ;
        if(rand_ind_g2)  xb[i, which_dim_ind_reg_sort[3]]  += dot_product(zLong_3[i,    ], zeta_3[,   1]) ;
        if(rand_ind_g3)  xb[i, which_dim_ind_reg_sort[4]]  += dot_product(zLong_4[i,    ], zeta_4[,   1]) ;
        if(rand_ind_g4)  xb[i, which_dim_ind_reg_sort[5]]  += dot_product(zLong_5[i,    ], zeta_5[,   1]) ;
        if(rand_ind_g5)  xb[i, which_dim_ind_reg_sort[6]]  += dot_product(zLong_6[i,    ], zeta_6[,   1]) ;
        if(rand_ind_g6)  xb[i, which_dim_ind_reg_sort[7]]  += dot_product(zLong_7[i,    ], zeta_7[,   1]) ;
        if(rand_ind_g7)  xb[i, which_dim_ind_reg_sort[8]]  += dot_product(zLong_8[i,    ], zeta_8[,   1]) ;
        if(rand_ind_g8)  xb[i, which_dim_ind_reg_sort[9]]  += dot_product(zLong_9[i,    ], zeta_9[,   1]) ;
        if(rand_ind_g9)  xb[i, which_dim_ind_reg_sort[10]] += dot_product(zLong_10[i,   ], zeta_10[,   1]);
        if(rand_ind_g10) xb[i, which_dim_ind_reg_sort[11]] += dot_product(zLong_11[i,   ], zeta_11[,   1]);
        if(rand_ind_g11) xb[i, which_dim_ind_reg_sort[12]] += dot_product(zLong_12[i,   ], zeta_12[,   1]);
        if(rand_ind_g12) xb[i, which_dim_ind_reg_sort[13]] += dot_product(zLong_13[i,   ], zeta_13[,   1]);
        if(rand_ind_g13) xb[i, which_dim_ind_reg_sort[14]] += dot_product(zLong_14[i,   ], zeta_14[,   1]);
        if(rand_ind_g14) xb[i, which_dim_ind_reg_sort[15]] += dot_product(zLong_15[i,   ], zeta_15[,   1]);
        if(rand_ind_g15) xb[i, which_dim_ind_reg_sort[16]] += dot_product(zLong_16[i,   ], zeta_16[,   1]);
        if(rand_ind_g16) xb[i, which_dim_ind_reg_sort[17]] += dot_product(zLong_17[i,   ], zeta_17[,   1]);
        if(rand_ind_g17) xb[i, which_dim_ind_reg_sort[18]] += dot_product(zLong_18[i,   ], zeta_18[,   1]);
        if(rand_ind_g18) xb[i, which_dim_ind_reg_sort[19]] += dot_product(zLong_19[i,   ], zeta_19[,   1]);
        if(rand_ind_g19) xb[i, which_dim_ind_reg_sort[20]] += dot_product(zLong_20[i,   ], zeta_20[,   1]);
        if(rand_ind_g20) xb[i, which_dim_ind_reg_sort[21]] += dot_product(zLong_21[i,   ], zeta_21[,   1]);
        if(rand_ind_g21) xb[i, which_dim_ind_reg_sort[22]] += dot_product(zLong_22[i,   ], zeta_22[,   1]);
        if(rand_ind_g22) xb[i, which_dim_ind_reg_sort[23]] += dot_product(zLong_23[i,   ], zeta_23[,   1]);
        if(rand_ind_g23) xb[i, which_dim_ind_reg_sort[24]] += dot_product(zLong_24[i,   ], zeta_24[,   1]);
        if(rand_ind_g24) xb[i, which_dim_ind_reg_sort[25]] += dot_product(zLong_25[i,   ], zeta_25[,   1]);
        if(rand_ind_g25) xb[i, which_dim_ind_reg_sort[26]] += dot_product(zLong_26[i,   ], zeta_26[,   1]);
        if(rand_ind_g26) xb[i, which_dim_ind_reg_sort[27]] += dot_product(zLong_27[i,   ], zeta_27[,   1]);
        if(rand_ind_g27) xb[i, which_dim_ind_reg_sort[28]] += dot_product(zLong_28[i,   ], zeta_28[,   1]);
        if(rand_ind_g28) xb[i, which_dim_ind_reg_sort[29]] += dot_product(zLong_29[i,   ], zeta_29[,   1]);
        if(rand_ind_g29) xb[i, which_dim_ind_reg_sort[30]] += dot_product(zLong_30[i,   ], zeta_30[,   1]);
        if(rand_ind_g30) xb[i, which_dim_ind_reg_sort[31]] += dot_product(zLong_31[i,   ], zeta_31[,   1]);
        if(rand_ind_g31) xb[i, which_dim_ind_reg_sort[32]] += dot_product(zLong_32[i,   ], zeta_32[,   1]);

      }

      c[i, ] = delta[jj[i], ] + db[i];
      
      nu[i] = dot_product(theta[nn[i], ] + xb[i, ], exp(col(alpha, jj[i])));

      if(any_eta3pl) {
        eta3pl[i] = (itype[i] == 3) ? eta3pl_l[find_eta3pl[i]] : 0.0;
      } else {
        eta3pl[i] = 0.0;
      }
    }

  }

}

model {
  // prior N(0,1) on general theta dimension (second-order theta)
  to_vector(theta) ~ normal(0, 1);
  // standard deviation of residual in the factor regression model:
  // theta_d = lambda_d * theta_g + {theta_resid_d}. Estimate one per dimension
  //sig_thetag_reg ~ normal(0, 1); //uniform(0, 1);
  // to_vector(theta_resid) ~ normal(0, sig_thetag_reg_rep);
  // for(d in 1:D) {
  //   to_vector(theta_resid[,d]) ~ normal(0, sig_thetag_reg[d]);
  // }
  // How does general theta load on the first-order thetas
  // theta_d = {lambda_d} * theta_g + theta_resid_d
  // lambda1 is restricted to range 0.4-1.0 using a transformation above
  // This is for identifiability
  // Otherwise lambda is scaled to correlation range -1,1

  //lambda ~ normal(0, 1);
  // lambda[1] ~ normal(0, 5);
  // for(i in 2:D) {
  //   lambda[i] ~ uniform(-1, 1);
  // }

  // If 2pl model (i.e., discrimination params), then sample from a normal(-0.5,1)
  // and transform to a lognormal distribution above since discriminations
  // must be positive. Transformation preferred so that regression is made 
  // possible in this context.
  // if(L) {
  alpha_l ~ normal(-0.5, 1.0); // should these priors be wider?
  if(LMean) {
    alpha_r_l ~ normal(-0.5, 1.0);
  }
  // }
  // difficulty parameters are always estimated and they have vague normal
  // priors
  delta_l ~ normal(0, 5);
  if(deltaMean) {
    delta_r_l ~ normal(0, 5);
  }
  // 3pl models have a guessing parameter. Beta prior so that support is between
  // 0 and 1
  if(any_eta3pl) {
    eta3pl_l ~ beta(1, 19); // 5,23  2,3  1,19
  }
  // fixed effects for theta regression model
  if(has_treg) {
    beta_l ~ normal(0, 5);
  }
  // independent random effects for theta regression model. Estimate vector
  // of random effects and standard deviation in each random effect distribution.
  if(any_rand_ind) {
                        zeta_l     ~ normal(0, 1); zeta_l_sd    ~ cauchy(0, 5);
    if(rand_ind_g1 )  { zeta_l_2   ~ normal(0, 1); zeta_l_sd_2  ~ cauchy(0, 5); }
    if(rand_ind_g2 )  { zeta_l_3   ~ normal(0, 1); zeta_l_sd_3  ~ cauchy(0, 5); }
    if(rand_ind_g3 )  { zeta_l_4   ~ normal(0, 1); zeta_l_sd_4  ~ cauchy(0, 5); }
    if(rand_ind_g4 )  { zeta_l_5   ~ normal(0, 1); zeta_l_sd_5  ~ cauchy(0, 5); }
    if(rand_ind_g5 )  { zeta_l_6   ~ normal(0, 1); zeta_l_sd_6  ~ cauchy(0, 5); }
    if(rand_ind_g6 )  { zeta_l_7   ~ normal(0, 1); zeta_l_sd_7  ~ cauchy(0, 5); }
    if(rand_ind_g7 )  { zeta_l_8   ~ normal(0, 1); zeta_l_sd_8  ~ cauchy(0, 5); }
    if(rand_ind_g8 )  { zeta_l_9   ~ normal(0, 1); zeta_l_sd_9  ~ cauchy(0, 5); }
    if(rand_ind_g9 )  { zeta_l_10  ~ normal(0, 1); zeta_l_sd_10 ~ cauchy(0, 5); }
    if(rand_ind_g10 ) { zeta_l_11  ~ normal(0, 1); zeta_l_sd_11 ~ cauchy(0, 5); }
    if(rand_ind_g11 ) { zeta_l_12  ~ normal(0, 1); zeta_l_sd_12 ~ cauchy(0, 5); }
    if(rand_ind_g12 ) { zeta_l_13  ~ normal(0, 1); zeta_l_sd_13 ~ cauchy(0, 5); }
    if(rand_ind_g13 ) { zeta_l_14  ~ normal(0, 1); zeta_l_sd_14 ~ cauchy(0, 5); }
    if(rand_ind_g14 ) { zeta_l_15  ~ normal(0, 1); zeta_l_sd_15 ~ cauchy(0, 5); }
    if(rand_ind_g15 ) { zeta_l_16  ~ normal(0, 1); zeta_l_sd_16 ~ cauchy(0, 5); }
    if(rand_ind_g16 ) { zeta_l_17  ~ normal(0, 1); zeta_l_sd_17 ~ cauchy(0, 5); }
    if(rand_ind_g17 ) { zeta_l_18  ~ normal(0, 1); zeta_l_sd_18 ~ cauchy(0, 5); }
    if(rand_ind_g18 ) { zeta_l_19  ~ normal(0, 1); zeta_l_sd_19 ~ cauchy(0, 5); }
    if(rand_ind_g19 ) { zeta_l_20  ~ normal(0, 1); zeta_l_sd_20 ~ cauchy(0, 5); }
    if(rand_ind_g20 ) { zeta_l_21  ~ normal(0, 1); zeta_l_sd_21 ~ cauchy(0, 5); }
    if(rand_ind_g21 ) { zeta_l_22  ~ normal(0, 1); zeta_l_sd_22 ~ cauchy(0, 5); }
    if(rand_ind_g22 ) { zeta_l_23  ~ normal(0, 1); zeta_l_sd_23 ~ cauchy(0, 5); }
    if(rand_ind_g23 ) { zeta_l_24  ~ normal(0, 1); zeta_l_sd_24 ~ cauchy(0, 5); }
    if(rand_ind_g24 ) { zeta_l_25  ~ normal(0, 1); zeta_l_sd_25 ~ cauchy(0, 5); }
    if(rand_ind_g25 ) { zeta_l_26  ~ normal(0, 1); zeta_l_sd_26 ~ cauchy(0, 5); }
    if(rand_ind_g26 ) { zeta_l_27  ~ normal(0, 1); zeta_l_sd_27 ~ cauchy(0, 5); }
    if(rand_ind_g27 ) { zeta_l_28  ~ normal(0, 1); zeta_l_sd_28 ~ cauchy(0, 5); }
    if(rand_ind_g28 ) { zeta_l_29  ~ normal(0, 1); zeta_l_sd_29 ~ cauchy(0, 5); }
    if(rand_ind_g29 ) { zeta_l_30  ~ normal(0, 1); zeta_l_sd_30 ~ cauchy(0, 5); }
    if(rand_ind_g30 ) { zeta_l_31  ~ normal(0, 1); zeta_l_sd_31 ~ cauchy(0, 5); }
    if(rand_ind_g31 ) { zeta_l_32  ~ normal(0, 1); zeta_l_sd_32 ~ cauchy(0, 5); }
  }
  // Correlated random effects for theta regression. Matrix of uncertainties.
  // Multivariate normal prior.
  if(any_rand_cor) {
    tau ~ cauchy(0, 2.5);
    Omega ~ lkj_corr(1);
    for(i in 1:u_Lzeta_cor) {
      zeta_c[i] ~ multi_normal(zeros_Lzeta_cor, quad_form_diag(Omega, tau));
    }
    if(rand_cor_g1) {
      tau_2 ~ cauchy(0, 2.5);
      Omega_2 ~ lkj_corr(1);
      for(i in 1:u_Lzeta_cor_2) {
        zeta_c_2[i] ~ multi_normal(zeros_Lzeta_cor_2, quad_form_diag(Omega_2, tau_2));
      }
    }
    if(rand_cor_g2) {
      tau_3 ~ cauchy(0, 2.5);
      Omega_3 ~ lkj_corr(1);
      for(i in 1:u_Lzeta_cor_3) {
        zeta_c_3[i] ~ multi_normal(zeros_Lzeta_cor_3, quad_form_diag(Omega_3, tau_3));
      }
    }
  }
  // independent random effects for alpha regression model. Estimate vector
  // of random effects and standard deviation in each random effect distribution.
  if(any_rand_ind_a) {
    aeta_l  ~ normal(0, 1);
    aeta_l_sd ~ cauchy(0, 2); // may need to hammer down on this?
  }
  // Correlated random effects for alpha regression. Matrix of uncertainties.
  // Multivariate normal prior.
  if(any_rand_cor_a) {
    tau_a ~ normal(0, 2);
    Omega_a ~ lkj_corr(1);
    for(i in 1:u_Laeta_cor) {
      aeta_c[i] ~ multi_normal(zeros_Laeta_cor, quad_form_diag(Omega_a, tau_a));
    }
  }
  // independent random effects for delta regression model. Estimate vector
  // of random effects and standard deviation in each random effect distribution.
  if(any_rand_ind_d) {
    deta_l  ~ normal(0, 1);
    deta_l_sd ~ cauchy(0, 5);
  }
  // Correlated random effects for delta regression. Matrix of uncertainties.
  // Multivariate normal prior.
  if(any_rand_cor_d) {
    tau_d ~ cauchy(0, 2.5);
    Omega_d ~ lkj_corr(1);
    for(i in 1:u_Ldeta_cor) {
      deta_c[i] ~ multi_normal(zeros_Ldeta_cor, quad_form_diag(Omega_d, tau_d));
    }
  }
  // increment the log likelihood (ordered logistic distribution customized for IRT)
  target += ordered_logistic_log_irt_vec(y, nu, c, eta3pl, Ncategi_jj, N_long, itype);
}
