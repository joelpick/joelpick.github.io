
data {
  int<lower=0> N;                // number of observations
  int<lower=1> J;                // number of predictors
  int<lower=0> Nped;             // number of individuals in pedigree
  int<lower=0> N_NoParents;
  int<lower=0> N_ParentsOffspring;
  int<lower=0> N_ParentsNoOffspring;

  vector[N] y;                    // response  
  matrix[N, J] X;                 // matrix of predictors
  vector[Nped] MSV;               // Mendelian sampling variance

  int<lower=0, upper=Nped> animal[N];     // animal id relating y to animal in pedigree (refers to row number of dam and sire)
  int<lower=0, upper=Nped> dam[Nped];       // dam id in order of animal id
  int<lower=0, upper=Nped> sire[Nped];      // sire id  
  int<lower=0, upper=Nped> NoParents[N_NoParents];
  int<lower=0, upper=Nped> ParentsOffspring[N_ParentsOffspring];
  int<lower=0, upper=Nped> ParentsNoOffspring[N_ParentsNoOffspring];
}

transformed data{
  vector[Nped] MSsd = sqrt(MSV);
  // Compute, thin, and then scale QR decomposition
  matrix[N, J] Q = qr_Q(X)[, 1:J] * sqrt(N - 1.0);
  matrix[J, J] R = qr_R(X)[1:J, ] / sqrt(N - 1.0);
  matrix[J, J] R_inv = inverse(R);
}

parameters {
  vector[J] beta_tilde;              // fixed effects
  vector [N_NoParents] A_NoParents;  // scaled breeding value
  vector [N_ParentsOffspring] A_ParentsOffspring;
  vector [N_ParentsNoOffspring] A_ParentsNoOffspring;
  real <lower=0> sigma_A;            // sd of A
  real <lower=0> sigma_E;            // residual SD
}

model {
  vector [N] mu;
  vector [Nped] A_scaled;     // scaled breeding values
  vector [Nped] A;            // breeding values
  real A_mu;
  vector [N_ParentsNoOffspring] A_ParentsNoOffspring_mean;
 
  int PO_i;
  
  A_scaled[1] = 0;

  A_NoParents ~ normal( 0 , 1 );
  A_scaled[NoParents] = A_NoParents;

  for (i in 1:N_ParentsOffspring)
  {
    PO_i = ParentsOffspring[i];
 
    A_mu = (A_scaled[dam[PO_i]] + A_scaled[sire[PO_i]])*0.5;
    A_ParentsOffspring[i] ~ normal( A_mu, MSsd[PO_i]);
    A_scaled[PO_i] = A_ParentsOffspring[i];
  }
  
   A_ParentsNoOffspring_mean = (A_scaled[dam[ParentsNoOffspring]] + A_scaled[sire[ParentsNoOffspring]])*0.5;
  A_ParentsNoOffspring ~ normal( A_ParentsNoOffspring_mean , MSsd[ParentsNoOffspring]);
  A_scaled[ParentsNoOffspring] = A_ParentsNoOffspring ;
  
  A = A_scaled * sigma_A;
  mu = Q * beta_tilde + A[animal];
 
  y ~ normal(mu, sigma_E) ;
  
  sigma_A ~ cauchy(0,5);
  sigma_E ~ cauchy(0,5);
}

generated quantities {
  vector[J] beta = R_inv * beta_tilde; // coefficients on X
}
