//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  int<lower =0> nspc; // n is the number of species in the data set
  vector[N] y;
  vector[N] logb;
  vector[N] inv_temp;
  int spc[N];
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu_Ao;
  real logsigma_Ao;
  vector[nspc] Ao_raw;
  real Eo;
  real n;
  real<lower=0> sigma;
}

transformed parameters {
 real Ao[nspc];
 real sigma_Ao;
 
 sigma_Ao = exp(logsigma_Ao);
 
 for (i in 1:nspc) Ao[i] = mu_Ao + Ao_raw[i] * sigma_Ao;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  real mu[N];
  
  // Priors
  for (i in 1:nspc) Ao_raw[i] ~ normal(0,1);
  logsigma_Ao ~ cauchy(0, 2);
  Eo ~ cauchy(0, 10);
  n ~ cauchy(0, 10); 
  sigma ~ cauchy(0, 10);
  mu_Ao ~ cauchy(20, 50);
  
  for (i in 1:N) mu[i] = Ao[spc[i]] + Eo * inv_temp[i] + n * logb[i];
  
  y ~ normal(mu, sigma);
}

