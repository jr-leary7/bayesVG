data {
  int<lower=1> M;  // number of spots
  int<lower=1> N;  // number of gene-spot pairs in long dataframe
  int<lower=1> G;  // number of genes
  int<lower=1> k;  // number of basis functions used to approximate GP
  array[N] int<lower=1, upper=M> spot_id;  // unique ID for each spot
  array[N] int<lower=1, upper=G> gene_id;  // unique ID for each gene
  matrix[M, k] phi;  // matrix of basis functions used to approximate GP
  vector[N] y;  // vector of normalized, scaled gene expression used as response variable
}

parameters {
  real beta0;  // global intercept
  matrix[G, k] alpha;  // matrix of gene-specific coefficients for each basis function
  real<lower=0> sigma_y;  // observation noise of response variable
  vector<lower=0>[G] amplitude;  // vector of gene-specific amplitudes of the approximate GP
  real mu_amplitude;  // mean for the amplitude
  real<lower=0> sigma_amplitude;  // SD for the amplitude
  real mu_beta0;  // mean for the gene-specific intercepts
  real<lower=0> sigma_beta0;  // SD for the gene-specific intercepts
  vector[k] mu_alpha;  // vector of means for the basis function coefficients
  vector<lower=0>[k] sigma_alpha;  // vector of SDs for the basis function coefficients
}

model {
  mu_beta0 ~ normal(0, 3);
  sigma_beta0 ~ normal(0, 2);
  mu_alpha ~ normal(0, 2);
  sigma_alpha ~ normal(0, 2);
  sigma_y ~ normal(0, 2);
  mu_amplitude ~ normal(0, 2);
  sigma_amplitude ~ normal(0, 1);
  for (i in 1:G) {
    for (j in 1:k) {
      alpha[i, j] ~ normal(mu_alpha[j], sigma_alpha[j]);
    }
    amplitude[i] ~ lognormal(mu_amplitude, sigma_amplitude);
  }
  for (i in 1:N) {
    int g = gene_id[i];
    int p = spot_id[i];
    real w_i = dot_product(phi[p], alpha[g]);
    y[i] ~ normal(beta0 + square(amplitude[g]) * w_i, sigma_y);
  }
}
