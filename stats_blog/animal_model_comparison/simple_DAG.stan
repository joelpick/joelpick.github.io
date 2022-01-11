// Only works when phenotype is known for ALL individuals in pedigree, which is almost never true


data {
  int<lower=0> N;                     // number of observations
  int<lower=1> J;                // number of predictors
  int<lower=0> Nped;                     // number of individuals in pedigree
  
  vector[N] y;                        // response
  matrix[N, J] X;                 // matrix of predictors
  vector[Nped] MSV;               // Mendelian sampling variance

  int<lower=0, upper=Nped> animal[N];     // animal id relating y to animal in pedigree (refers to row number of dam and sire)
  int<lower=0, upper=Nped> dam[Nped];       // dam id in order of animal id
  int<lower=0, upper=Nped> sire[Nped];      // sire id  
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

  vector [Nped] A_scaled;          // scaled breeding value
  real <lower=0> sigma_A;          // sd of A
  real <lower=0> sigma_E;          // residual SD
}

// transformed parameters {
// }

model {
  real mu; // intercept
  
  A_scaled[1] = 0;

  for (i in 2:Nped)
  {  
    A_scaled[i] ~ normal( (A_scaled[dam[i]] + A_scaled[sire[i]])*0.5 , MSsd[i]);
  }

  mu = Q * beta_tilde + sigma_A * A_scaled[animal];
  y ~ normal(mu, sigma_E) ;
  sigma_A ~ cauchy(0,5);
  sigma_E ~ cauchy(0,5);
}

generated quantities {
  vector[J] beta = R_inv * beta_tilde; // coefficients on X
}
