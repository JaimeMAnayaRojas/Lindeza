// Autoregressive model for the Anas  larval data
data {
  int<lower=0> Nc;
  int<lower=0> Nm;
  int<lower=0> Nu;
  
  vector[Nc] Y_c;
  vector[Nc] t_c;
  int L_c[Nc];
  vector[Nc] Xa_c;

  vector[Nm] Y_m;
  vector[Nm] t_m;
  int L_m[Nm];
  vector[Nm] Xa_m;
  
  vector[Nu] Y_u;
  vector[Nu] t_u;
  int L_u[Nu];
  vector[Nu] Xa_u;
  
}
parameters {
  real alpha_c;
  real beta_c;
  real beta_Ac;
    
  real alpha_m;
  real beta_m;
  real beta_Am;
  
  real alpha_u;
  real beta_u;
  real beta_Au;
  
  real<lower=0> sigma;
  real<lower=0> sigma_L;
  
 vector[9] v_alpha;
  
  
  
}
model {

  v_alpha ~ normal( 0 , sigma_L );
  
  
  
  Y_c[2:Nc] ~ normal(alpha_c + beta_c * Y_c[1:(Nc - 1)] + beta_Ac * Xa_c[1:(Nc - 1)] + v_alpha[L_c[1:(Nc - 1)]], sigma);
  Y_m[2:Nm] ~ normal(alpha_m + beta_m * Y_m[1:(Nm - 1)] + beta_Am * Xa_c[1:(Nm - 1)] + v_alpha[L_m[1:(Nm - 1)]], sigma);
  Y_u[2:Nu] ~ normal(alpha_u + beta_u * Y_u[1:(Nu - 1)] + beta_Au * Xa_c[1:(Nu - 1)] + v_alpha[L_u[1:(Nu - 1)]], sigma);
  
  
}
