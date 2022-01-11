
data {
  int<lower=0> N;                     // number of observations
  int<lower=1> J; // number of predictors
  int<lower=0> Nped;                     // number of individuals in pedigree
  int<lower=0> N_NoParentsOffspring;
  int<lower=0> N_ParentsOffspring;
  int<lower=0> N_NoOffspring;

  vector[N] y;                        // response  
  matrix[N, J] X; // matrix of predictors
  vector[Nped] MSV;                      // Mendelian sampling variance

  int<lower=0, upper=Nped> animal[N];     // animal id relating y to animal in pedigree (refers to row number of dam and sire)
  int<lower=0, upper=Nped> dam[Nped];       // dam id in order of animal id
  int<lower=0, upper=Nped> sire[Nped];      // sire id  
  int<lower=0, upper=Nped> NoParentsOffspring[N_NoParentsOffspring];
  int<lower=0, upper=Nped> ParentsOffspring[N_ParentsOffspring];
  int<lower=0, upper=Nped> NoOffspring[N_NoOffspring];
}

transformed data{
  // Compute, thin, and then scale QR decomposition
  matrix[N, J] Q = qr_Q(X)[, 1:J] * sqrt(N - 1.0);
  matrix[J, J] R = qr_R(X)[1:J, ] / sqrt(N - 1.0);
  matrix[J, J] R_inv = inverse(R);

  vector[Nped] MSsd = sqrt(MSV);
  vector[Nped] MSV_red;
  MSV_red[1] = 0; 
  MSV_red[NoParentsOffspring] = rep_vector(0.0,N_NoParentsOffspring);
  MSV_red[ParentsOffspring] = rep_vector(0.0,N_ParentsOffspring);
  MSV_red[NoOffspring] = MSV[NoOffspring];
}


parameters {
  vector[J] beta_tilde; // fixed effects
  vector [N_NoParentsOffspring] A_NoParentsOffspring;          // scaled breeding value
  vector [N_ParentsOffspring] A_ParentsOffspring;
  // vector [N_NoOffspring] A_NoOffspring;
  real <lower=1e-16> sigma_A;          // sd of A
  real <lower=1e-16> sigma_E;          // residual SD
}

model {
  vector [N] mu;
  vector [Nped] A_scaled;          // scaled breeding value
  vector [Nped] sigma_E_mix;
  int PO_i;
  real A_mu;
 
  A_scaled[1] = 0;

  A_NoParentsOffspring ~ normal( 0 , 1 );
  A_scaled[NoParentsOffspring] = A_NoParentsOffspring;

  for (i in 1:N_ParentsOffspring)
  {
    PO_i = ParentsOffspring[i];
    A_mu = (A_scaled[dam[PO_i]] + A_scaled[sire[PO_i]])*0.5;
    A_ParentsOffspring[i] ~ normal( A_mu, MSsd[PO_i]);  
    A_scaled[PO_i] = A_ParentsOffspring[i];
  }

  A_scaled[NoOffspring] = (A_scaled[dam[NoOffspring]] + A_scaled[sire[NoOffspring]])*0.5;
  
  sigma_A ~ cauchy(0,5);
  sigma_E ~ cauchy(0,5);

  mu = Q * beta_tilde + A_scaled[animal] * sigma_A;
  sigma_E_mix = sqrt(MSV_red * sigma_A^2 + sigma_E^2);

  y ~ normal(mu, sigma_E_mix[animal]) ;
  
}

generated quantities {
      vector[J] beta = R_inv * beta_tilde; // coefficients on x
}
